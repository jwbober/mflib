#include "cuspforms_modp.h"
#include "classnumbers.h"
#include "mfformat.h"

#include "flint/nmod_poly.h"
#include "flint/fmpz_poly.h"

#include "flint/NTL-interface.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <thread>

#include <unordered_map>

#include "NTL/ZZX.h"
#include "NTL/ZZXFactoring.h"

#include "acb_poly.h"

#include "arb_fmpz_poly.h"

#include "arb-extras.h"
#include "flint-extras.h"

#include "ThreadPool/ThreadPool.h"

using namespace std;

static int nthreads = 1;

const char * usage =
" level weight mfdb polydb [overwrite] [nthreads]\n"
"\n"
"Given an input file with embeddings of newforms, compute the decomposition\n"
"of each character space into Galois orbits over Q and put this information\n"
"into the output file polydb.\n"
"\n"
"This program expects that if\n"
"there is data about one character in a Galois orbit in the database, then there\n"
"will be data for all of the characters. So it is a little fragile in a few ways.\n"
"\n"
"By default, overwrite is set to 0. If it is nonzero, the contents of the\n"
"hecke_polynomials table in the output database (if it exists) will be deleted.\n";

struct pair_hash {
    std::size_t operator()(const pair<int,int> &x) const {
        return x.first + ((long)x.second << 32);
    }
};

struct pair_equal {
    bool operator()(const pair<int,int> &x, const pair<int,int> &y) const {
        return x == y;
    }
};

int hecke_polynomial_modular_approximation(fmpz_poly_t out, int level, int weight, int chi_number, vector<int> hecke_operator, fmpz_t precision) {
    int verbose2 = 0;
    int verbose = 1;
    // compute the characteristic polynomial of the hecke operator sum c_n T_n,
    // where the vector hecke_operator is [c_2, c_3, c_4, ..., c_k], acting on
    // the space
    //
    // (DIRECT SUM) S_weight^new (level, chi)
    //
    // where the direct sum is taken over the Galois orbit of characters chi
    // which contains the character numbered chinumber, accurate modulo at
    // least the input precision. On output, the variable precision will be
    // modified to contain the actual precision of the computation.
    //
    // returns 1 on success, 0 if somethign went wrong (like bad input),
    // in which case the output variables may not or may not have changed.
    //

    if(hecke_operator.size() < 1) return 0;

    DirichletGroup G(level);
    if(GCD(level, chi_number) != 1) return 0;
    DirichletCharacter chi = G.character(chi_number);
    if(chi.is_even() && weight % 2 == 1) return 0;
    if(!chi.is_even() && weight % 2 == 0) return 0;

    set<long> _orbit = chi.galois_orbit();
    vector<long> orbit(_orbit.begin(), _orbit.end());
    int orbitsize = orbit.size();
    if(verbose) {
        cout << "in hecke_polynomial_modular_approximation: Galois orbit size is " << orbitsize << endl;
    }

    fmpz_poly_t hecke1;
    fmpz_poly_t hecke2;

    fmpz_poly_init(hecke1);
    fmpz_poly_init(hecke2);

    fmpz_poly_zero(hecke1);
    fmpz_poly_zero(hecke2);

    nmod_poly_t hecke_poly_modp;

    fmpz_t a;
    fmpz_t b;
    fmpz_t modulus;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(modulus);
    fmpz_set_ui(modulus, 1);


    // we'll start the computation at 2^50
    long p = 1125899906842624;
    do {
        fmpz_poly_set(hecke1, hecke2);
        // hecke1 now contains the polynomial we've computed so far, accurate modulo modulus
        p += 1;

        // This is just a convenient way of choosing p. We don't actually
        // do anything with the first space constructed, but we do get back
        // a prime p which satisfies the appropriate congruence conditions
        // and will remain stable as we vary the character in the orbit.
        cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
        p = S->p;
        if(verbose) cout << p << endl;
        if(verbose) {
            cout << fmpz_bits(modulus) << " " << fmpz_bits(precision) << endl;
        }
        nmod_poly_init(hecke_poly_modp, p);
        nmod_poly_one(hecke_poly_modp);
        nmod_poly_t * local_polys = new nmod_poly_t[orbitsize];
        for(int k = 0; k < orbitsize; k++) {
            nmod_poly_init(local_polys[k], p);
        }
        //nmod_poly_init(f, p);

        //std::thread * threads = new std::thread[orbitsize];
        //for(long m : orbit) {
        ThreadPool * pool = new ThreadPool(min(nthreads, orbitsize));
        for(int k = 0; k < orbitsize; k++) {
            long m = orbit[k];
            chi = G.character(m);
            cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
            p = S->p;


            pool->enqueue( [local_polys, k, S, p, &hecke_operator](){
            //threads[k] = std::thread( [local_polys, k, S, p, &hecke_operator](){
                nmod_mat_t hecke_mat;
                nmod_mat_init(hecke_mat, S->new_dimension(), S->new_dimension(), p);
                for(int l = 2; l < hecke_operator.size() + 2; l++) {
                    int cl = hecke_operator[l - 2];
                    if(cl == 0) continue;
                    nmod_mat_t Tl;
                    S->hecke_matrix(Tl, l);
                    nmod_mat_scalar_mul_add(hecke_mat, hecke_mat, cl, Tl);
                    nmod_mat_clear(Tl);
                }
                nmod_mat_charpoly(local_polys[k], hecke_mat);
                //nmod_poly_mul(hecke_poly_modp, hecke_poly_modp, f);
                nmod_mat_clear(hecke_mat);
                });
        }
        delete pool;
        for(int k = 0; k < orbitsize; k++) {
            // threads[k].join();
            nmod_poly_mul(hecke_poly_modp, hecke_poly_modp, local_polys[k]);
            nmod_poly_clear(local_polys[k]);
        }

        delete [] local_polys;
        //delete [] threads;
        int degree = nmod_poly_degree(hecke_poly_modp);
        for(int k = 0; k <= degree; k++) {
            fmpz_poly_get_coeff_fmpz(a, hecke1, k);
            unsigned long z = nmod_poly_get_coeff_ui(hecke_poly_modp, k);
            fmpz_CRT_ui(a, a, modulus, z, p, 1);
            fmpz_poly_set_coeff_fmpz(hecke2, k, a);
        }
        fmpz_mul_ui(modulus, modulus, p);

        //nmod_poly_clear(f);
        nmod_poly_clear(hecke_poly_modp);
        clear_cuspforms_modp();
    //} while(!fmpz_poly_equal(hecke1, hecke2));
    } while(fmpz_cmp(modulus, precision) < 0);

    fmpz_set(precision, modulus);
    fmpz_poly_set(out, hecke2);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(modulus);

    fmpz_poly_clear(hecke1);
    fmpz_poly_clear(hecke2);
    return 1;
}

static void twist_poly(arb_poly_t r1,arb_poly_t r2,const fmpz_poly_t f,slong a,slong q, slong prec) {
    static arb_t c,s,t;
    static int init;
    slong d,n;

    if (!init) {
        arb_init(c); arb_init(s);
        arb_init(t);
        init = 1;
    }
    d = fmpz_poly_degree(f);
    arb_poly_truncate(r1,d+1);
    arb_poly_truncate(r2,d+1);
    for (;d>=0;d--) {
        n = a*d % q;
        if (n < 0) n += q;
        arb_set_si(t,2*n);
        arb_div_si(t,t,q,prec);
        arb_sin_cos_pi(s,c,t,prec);

        arb_mul_fmpz(t,c,fmpz_poly_get_coeff_ptr(f,d),prec);
        if (4*n == q || 4*n == 3*q) arb_zero(t);
        arb_poly_set_coeff_arb(r1,d,t);

        arb_mul_fmpz(t,s,fmpz_poly_get_coeff_ptr(f,d),prec);
        if (!n || 2*n == q) arb_zero(t);
        arb_poly_set_coeff_arb(r2,d,t);
    }
}

static int phase(acb_srcptr a,int order, slong prec) {
    static arb_t t,pi;
    static fmpz_t z;
    static int init;
    int res;

    if (order == 1) return 0;
    if (!init) {
        arb_init(t); arb_init(pi);
        fmpz_init(z);
        init = 1;
    }
    acb_arg(t,a,prec);
    arb_const_pi(pi,prec);
    arb_div(t,t,pi,prec);
    arb_mul_si(t,t,order,prec);
    if (!arb_get_unique_fmpz(z,t)) {
        printf("ehh cannot determine phase in arbgcd\n");
        cout << order << endl;
        acb_printd(a, 10);
        cout << endl;
        arb_printd(t, 10);
        cout << endl;
        return -1;
    }
    res = fmpz_get_si(z) % order;
    if (res < 0) res += order;
    return res;
}


int fmpz_vec_cmp(fmpz * x, fmpz * y, slong len) {
    for(int k = 0; k < len; k++) {
        int z = fmpz_cmp(x + k, y + k);
        if(z != 0) return z;
    }
    return 0;
}

bool match_eigenvalues_rotate_gcd(int * res,
                                 const fmpz_poly_struct * poly,
                                 int npolys,
                                 const acb_ptr a,
                                 int ncoeffs,
                                 int order,
                                 slong prec) {
    static arb_poly_t g;
    static arb_t t;
    static acb_t x,z;
    static int init;
    arb_poly_t *f;
    int i,j,u,phi,s;
    bool retval = true;

    if (!init) {
        arb_poly_init(g); arb_init(t);
        acb_init(x); acb_init(z);
        init = 1;
    }

    f = (arb_poly_t *)malloc(npolys*sizeof(f[0]));
    for (i=0;i<npolys;i++)
        arb_poly_init(f[i]);

    for (u=0,phi=0;u<order;u++)
        if (GCD(u,order) == 1) phi++;
    for (j=0;j<ncoeffs;j++)
        res[j] = npolys;
    for (u=0;u<order;u++)
        if (GCD(u,order) == 1) {
            acb_set_si(z,-u);
            acb_div_si(z,z,order,prec);
            acb_exp_pi_i(z,z,prec);
            for (i=0;i<npolys;i++) {
                //cout << order << endl;
                //fmpz_poly_print_pretty(poly + i, "x");
                //cout  << endl;
                twist_poly(f[i],g,poly + i,u,2*order, prec);
                //arb_poly_printd(f[i], 10);
                //cout << endl;
                arb_poly_printd(g, 10);
                //cout << endl;
                //cout << endl;
                // It may be the case that f[i] or g is zero, for example when all the
                // roots or already real, or in a more complicated exmaple, when the input
                // polynomial is something like x^4 + 100.
                //
                // In this second example the degree of the gcd will not be what we
                // expect generically, but we're ok because we know what the gcd is
                // because we know that one of the polynomials is exactly zero and
                // we only need to know the degree to compute the gcd/
                if(arb_poly_is_zero(f[i])) {
                    arb_poly_set(f[i], g);
                }
                else if(arb_poly_is_zero(g)) {
                    // do nothing
                }
                else {// neither is zero, so do the gcd
                    if (arb_poly_gcd(f[i],f[i],g,fmpz_poly_degree(poly + i)/phi, prec) < 0) {
                        cout << "ehh problem with arb_poly_gcd" << endl;
                        retval = false;
                        goto cleanup;
                    }
                }
            }
        for (j=0;j<ncoeffs;j++) {
            int thisphase = 0;
            if(acb_contains_zero(a + j)) {
                thisphase = u;
            }
            else {
                thisphase = phase(a + j, order, prec);
            }
            if(thisphase < 0) {
                retval = false;
                goto cleanup;
            }
            if (thisphase == u) {
                acb_mul(x,a + j,z,prec);
                for (i=0;i<npolys;i++) {
                    s = arb_poly_sign_change(f[i],acb_realref(x),prec);
                    if (!s || ((s < 0) && (res[j] < npolys))) {
                        arb_poly_printd(f[i], 10); cout << endl;
                        cout << x << " " << a + j << endl;
                        cout << "ehh found extra roots or couldn't find a sign change" << endl;
                        retval = false;
                        goto cleanup;
                    }
                    if (s < 0) res[j] = i;
                }
#if 0
printf("%d ",res[j]);
acb_printd(a + j,15);
printf("\n");
#endif
            }
        }
    }

    for (j=0;j<ncoeffs;j++) {
        if (res[j] == npolys) {
            cout << "ehh problem with root " << j << endl;
            retval = false;
        }
    }

cleanup:
    for (i=0;i<npolys;i++)
        arb_poly_clear(f[i]);
    free(f);
    return retval;
}

string format_hecke_operator(const vector<int> &hecke_operator) {
    stringstream ss;
    bool anything_so_far = false;
    for(int k = 0; k < hecke_operator.size(); k++) {
        int c = hecke_operator[k];
        if(c == 0) continue;
        if(anything_so_far) ss << " + ";
        if(c != 1) ss << c << ".";
        ss << "T" << k + 2;
        anything_so_far = true;
    }
    return ss.str();
}

int find_unique_eigenvalues(vector<int> &hecke_operator,
                            acb_ptr eigenvalues,
                            const int dimension,
                            const int full_dimension,
                            const unordered_map<pair<int,int>, acb_ptr, pair_hash, pair_equal> &coeff_table,
                            const vector<int> &chiorbit,
                            const int level,
                            const bool require_single_operator,
                            const int max_operator,
                            slong prec) {
    hecke_operator.clear();
    for(int k = 0; k < full_dimension; k++) {
        acb_zero(eigenvalues + k);
    }
    bool unique_eigenvalues = false;
    int k = 2;
    while(!unique_eigenvalues && (k <= max_operator)) {
        // we sometimes really prefer to compute the minimal polynomial of just
        // one hecke-eigenvalue, so we might try that first.

        if(require_single_operator) {
            while(GCD(k, level) != 1) {
                hecke_operator.push_back(0);
                k++;
            }
            hecke_operator.push_back(1);
        }
        else {
            hecke_operator.push_back(k - 1);
        }
        cout << "Trying Hecke operator " << format_hecke_operator(hecke_operator) << endl;
        for(int l = 0; l < chiorbit.size(); l++) {
            int chi_number = chiorbit[l];
            int chi_inverse = InvMod(chi_number, level);
            if(level == 1) chi_inverse = 1;
            for(int j = 0; j < dimension; j++) {
                if(chi_number <= chi_inverse) {
                    acb_ptr coeffs = coeff_table.at(make_pair(chi_number, j));
                    if(require_single_operator)
                        acb_set(eigenvalues + l * dimension + j, coeffs + k - 1);
                    else
                        acb_addmul_ui(eigenvalues + l * dimension + j, coeffs + k - 1, k - 1, prec);
                }
                else {
                    acb_ptr coeffs = coeff_table.at(make_pair(chi_inverse, j));
                    acb_conj(coeffs + k - 1, coeffs + k - 1);
                    if(require_single_operator) {
                        acb_set(eigenvalues + l * dimension + j, coeffs + k - 1);
                    }
                    else {
                        acb_addmul_ui(eigenvalues + l * dimension + j, coeffs + k - 1, k - 1, prec);
                    }
                    acb_conj(coeffs + k - 1, coeffs + k - 1);
                }
            }
        }
        unique_eigenvalues = true;
        for(int n = 0; n < dimension * chiorbit.size(); n++) {
            for(int m = n + 1; m < dimension * chiorbit.size(); m++) {
                if(acb_overlaps(eigenvalues + n, eigenvalues + m)) {
                    unique_eigenvalues = false;
                    break;
                }
            }
            if(!unique_eigenvalues && require_single_operator) {
                hecke_operator.pop_back();
                hecke_operator.push_back(0);
                break;
            }
        }
        k++;
    }
    if(unique_eigenvalues) {
        return k - 1;
    }
    hecke_operator.clear();
    return 0;
}

int get_integral_heckepoly(fmpz_poly_t zheckepoly, acb_poly_t heckepoly, int level, int weight, int chi_number, const vector<int> &hecke_operator, slong prec) {
    for(int k = 0; k <= acb_poly_degree(heckepoly); k++) {
        arb_zero(acb_imagref(acb_poly_get_coeff_ptr(heckepoly, k)));
    }
    if(!acb_poly_get_unique_fmpz_poly(zheckepoly, heckepoly)) {
        // we're going to have to up the accuracy with a mod p computation
        cout << "doing mod p computations to increase hecke polynomial precision" << endl;
        arb_t max_error, r;
        arb_init(max_error);
        arb_init(r);
        arb_zero(max_error);
        for(int k = 0; k < acb_poly_degree(heckepoly); k++) {
            arb_get_rad_arb(r, acb_realref(acb_poly_get_coeff_ptr(heckepoly, k)));
            if(arb_lt(max_error, r)) {
                arb_set(max_error, r);
            }
        }
        fmpz_t modular_precision;
        fmpz_init(modular_precision);

        arb_mul_2exp_si(max_error, max_error, 1);
        arb_get_ubound_fmpz(modular_precision, max_error, prec);
        fmpz_poly_t heckepoly_modular_approx;
        fmpz_poly_init(heckepoly_modular_approx);
        int result = hecke_polynomial_modular_approximation(heckepoly_modular_approx, level, weight, chi_number, hecke_operator, modular_precision);
        if(!result) return 0;
        result = acb_poly_get_unique_fmpz_poly_modular(zheckepoly, heckepoly, heckepoly_modular_approx, modular_precision);
        return result;
    }
    else {
        return 1;
    }
}

int main(int argc, char ** argv) {
    if(argc < 5) {
        cout << argv[0] << usage;
        return 0;
    }
    clock_t start_time = clock();
    long p = 1125899906842679l;
    string mfdbname = argv[3];
    string polydbname = argv[4];

    int level = atoi(argv[1]);
    int weight = atoi(argv[2]);
    set<long> chi_list;
    int * dimensions = NULL;
    init_classnumbers();
    load_factor_table();
    int prec = 0;

    sqlite3 * db;
    if(sqlite3_open(":memory:", &db) != SQLITE_OK) {
        cout << "ohno error opening in memory database." << endl;
        return 1;
    }
    int retval;

    retval = sqlite3_exec(db, "CREATE TABLE modforms (level INTEGER, weight INTEGER, chi INTEGER, orbit INTEGER, j INTEGER, prec INTEGER, exponent INTEGER, ncoeffs INTEGER, coefficients BLOB)", NULL, 0, NULL);
    retval = sqlite3_exec(db, "CREATE INDEX mf_level_weight_chi_j ON modforms (level, weight, chi, j);", NULL, 0, NULL);
    string attach_stmt = "ATTACH DATABASE '" + mfdbname + "' AS infiledb;";
    retval = sqlite3_exec(db, attach_stmt.c_str(), NULL, 0, NULL);
    if(retval != SQLITE_OK) {
        cout << "Error attaching input db." << endl;
        return 0;
    }
    string copy_stmt = "INSERT INTO modforms SELECT * from infiledb.modforms WHERE level=" + to_string(level) + " AND weight=" + to_string(weight);
    retval = sqlite3_exec(db, copy_stmt.c_str(), NULL, 0, NULL);
    if(retval != SQLITE_OK) {
        cout << "Error reading input db." << endl;
        return 0;
    }
    retval = sqlite3_exec(db, "DETACH DATABASE infiledb;", NULL, 0, NULL);

    cout << "input data loaded into memory." << endl;

    sqlite3 * polydb;
    if(sqlite3_open(polydbname.c_str(), &polydb) != SQLITE_OK) {
        cout << "ohno error opening output database." << endl;
        return 1;
    }
    polydb_init(polydb);

    int overwrite = 0;
    if(argc > 5) {
        overwrite = atoi(argv[5]);
    }
    if(overwrite) {
        cout << "deleting old database entries" << endl;
        sqlite3_exec(polydb, "DELETE FROM heckepolys;", NULL, 0, NULL);
    }
    if(argc > 6) {
        nthreads = atoi(argv[6]);
        cout << "setting nthreads = " << nthreads << endl;
    }

    mfheader * headers;
    int count = mfdb_contents(db, &headers);

    unordered_map<pair<int,int>, mfheader*, pair_hash, pair_equal> header_table;
    unordered_map<pair<int,int>, acb_ptr, pair_hash, pair_equal> coeff_table;

    for(int k = 0; k < count; k++) {
        level = headers[k].level;
        if(!dimensions) {
            if(level == 1)
                dimensions = new int[level + 1](); // stupid hack because the trivial character mod 1 is 1.1
                                                   // which is the only case where the character number might
                                                   // get as large as the level
            else
                dimensions = new int[level]();
        }
        header_table[make_pair((int)headers[k].chi, (int)headers[k].j)] = headers + k;
        weight = headers[k].weight;
        chi_list.insert(headers[k].chi);
        if(-headers[k].prec > prec) prec = -headers[k].prec;
        if(headers[k].j + 1 > dimensions[headers[k].chi]) dimensions[headers[k].chi] = headers[k].j+1;
        //cout << headers[k].level << " " << headers[k].weight << " " << headers[k].chi << " " << headers[k].j << endl;
    }
    if(prec < 0) prec = 200;
    prec = 30*prec;
    cout << "using working precision " << prec << endl;
    //free(headers);

    for(int chi = 1; chi <= level; chi++) {
        if(header_table.count(make_pair(chi, 0)) == 0) continue;
        for(int j = 0; j < dimensions[chi]; j++) {
            mfheader header;
            acb_ptr coeffs;
            int result = mfdb_get_entry(db, &header, &coeffs, level, weight, chi, j);
            if(!result) {
                cout << "ohno didn't find a database entry that we were expecting." << endl;
            }
            coeff_table[make_pair(chi, j)] = coeffs;
        }
    }

    DirichletGroup G(level);
    int * character_orbit_labels = new int[level];
    {
        int k = 0;
        for(auto S : G.ordered_galois_orbits()) {
            for(auto chi : S) {
                character_orbit_labels[chi] = k;
                cout << chi << " "  << k << endl;
            }
            k++;
        }
    }


    while(chi_list.size() > 0) {
        long chi_number = *chi_list.begin();
        DirichletCharacter chi = G.character(chi_number);
        set<long> _orbit = chi.galois_orbit();
        vector<int> orbit(_orbit.begin(), _orbit.end());
        vector<int> hecke_operator;
        for(long k : orbit) {
            //cout << k << " ";
            chi_list.erase(k);
        }

        mfheader header;
        acb_ptr coeffs;
        int dimension = dimensions[chi_number];
        int full_dimension = orbit.size() * dimension;
        //prec = full_dimension * 300;
        acb_ptr eigenvalues = _acb_vec_init(full_dimension);

        int unique_eigenvalues = find_unique_eigenvalues(hecke_operator,
                                        eigenvalues,
                                        dimension,
                                        full_dimension,
                                        coeff_table,
                                        orbit,
                                        level,
                                        false, // require single operator
                                        100, // max operator
                                        prec);


        cout << "for character " << chi_number << " trying hecke operator " << format_hecke_operator(hecke_operator) << endl;

        acb_poly_t heckepoly;
        acb_poly_init(heckepoly);
        cout << "forming product polynomial: " << clock() << endl;
        acb_poly_product_roots(heckepoly, eigenvalues, full_dimension, 300);
        cout << "finished forming product polynomial: " << clock() << endl;

        fmpz_poly_t zheckepoly;
        fmpz_poly_init(zheckepoly);

        if(!get_integral_heckepoly(zheckepoly, heckepoly, level, weight, chi_number, hecke_operator, prec)) {
            cout << "ohno something went wrong getting the hecke polynomial." << endl;
            continue;
        }

        NTL::ZZ content;
        NTL::ZZX ff;
        NTL::vec_pair_ZZX_long factors;
        fmpz_poly_get_ZZX(ff, zheckepoly);
        cout << "factoring polynomial: " << clock() << endl;
        NTL::factor(content, factors, ff);
        cout << "finished factoring polynomial: " << clock() << endl;
        fmpz_poly_t g;
        fmpz_poly_init(g);

        int nfactors = factors.length();

        acb_t t1;
        acb_init(t1);

        fmpz_poly_struct * fmpz_factors = new fmpz_poly_struct[nfactors];
        long * evaluation_bits = new long[nfactors];
        long max_coefficient_bits = 0;
        int * roots_found = new int[nfactors];

        for(int l = 0; l < nfactors; l++) {
            fmpz_poly_struct * g = fmpz_factors + l;
            fmpz_poly_init(g);
            fmpz_poly_set_ZZX(g, factors[l].a);
            long maxbits = 0;
            for(int k = 0; k <= fmpz_poly_degree(g); k++) {
                long bits = fmpz_bits(fmpz_poly_get_coeff_ptr(g, k));
                if(bits > maxbits) maxbits = bits;
            }
            if(maxbits > max_coefficient_bits) max_coefficient_bits = maxbits;
            evaluation_bits[l] = (maxbits + 10) * fmpz_poly_degree(g);
            if(evaluation_bits[l] < prec) evaluation_bits[l] = prec;
            roots_found[l] = 0;
        }

        cout << "factor degrees:";
        for(int l = 0; l < nfactors; l++) {
            cout << " " << fmpz_poly_degree(fmpz_factors + l);
        }
        cout << endl;

        bool matched_all_roots = false;
        int * matches = new int[full_dimension];
        if(nfactors == 1) {
            for(int k = 0; k < full_dimension; k++) matches[k] = 0;
            matched_all_roots = true;
        }

        if(!matched_all_roots) {
            for(int k = 0; k < full_dimension; k++) {
                matches[k] = -1;
            }
            // we first try to evaluate each of the polynomials as every possible
            // root. for each polynomial, if we find exactly the right number of
            // possible roots, then we know that we found them all, and we are
            // finished with that polynomial now.
            cout << "in root matching phase 1" << endl;
            matched_all_roots = true;

            bool * possible_root_table = new bool[full_dimension * nfactors]();
            for(int l = 0; l < nfactors; l++) {
                fmpz_poly_struct * g = fmpz_factors + l;
                if(roots_found[l] == fmpz_poly_degree(g)) continue;

                acb_ptr poly_values = _acb_vec_init(full_dimension);
                ThreadPool * pool = new ThreadPool(nthreads);
                for(int k = 0; k < full_dimension; k++) {
                    if(matches[k] != -1) continue;
                    pool->enqueue(arb_fmpz_poly_evaluate_acb, poly_values + k, g, eigenvalues + k, evaluation_bits[l]);

                }

                delete pool;

                vector<int> possible_roots;

                for(int k = 0; k < full_dimension; k++) {
                    if(matches[k] != -1) continue;
                    cout << "checking polynomial " << l << " (degree " << fmpz_poly_degree(g) << ") at root " << k;
                    if(acb_contains_zero(poly_values + k)) {
                        possible_root_table[k + full_dimension*l] = true;
                        possible_roots.push_back(k);
                        cout << " possible zero: " << possible_roots.size() << " out of " << fmpz_poly_degree(g) << " roots found." << endl;
                    }
                    else {
                        cout << " not a zero" << endl;
                    }
                }
                _acb_vec_clear(poly_values, full_dimension);

                if(possible_roots.size() == fmpz_poly_degree(g)) {
                    cout << "found all roots for polynomial " << l << endl;
                    for(int k : possible_roots) {
                        for(int l2 = 0; l2 < nfactors; l2++) {
                            if(l2 != l) possible_root_table[k + full_dimension*l2] = false;
                        }
                        matches[k] = l;
                    }
                    roots_found[l] = fmpz_poly_degree(g);
                }
                else {
                    cout << "ambiguity for polynomial " << l << endl;
                    matched_all_roots = false;
                }
            }

            bool progress_made = true;
            while(progress_made && (!matched_all_roots)) {
                progress_made = false;
                cout << "attempting phase 2 root matching..." << endl;
                // at this point it is possible that a later polynomial evaluation will have
                // removed the ambiguity in an earlier evaluation, and it is also possible
                // that we can certify some roots because they only evaluate to an interval
                // containing zero on a single polynomial, which is not something that we
                // have checked.
                //
                // (We still won't be able to handle the case where there are two
                //  eigenvalues which both could be roots of the same two polynomials.)

                // first we go through each eigenvalue and see if there is now
                // only one polynomial it might be a root of.

                matched_all_roots = true; // not really true, but assumed innocent
                                          // until found guilty...

                for(int k = 0; k < full_dimension; k++) {
                    if(matches[k] != -1) continue;
                    int full_possible_polys = 0;
                    int possible_poly = -1;
                    for(int l = 0; l < nfactors; l++) {
                        if(possible_root_table[k + full_dimension*l]) {
                            full_possible_polys += 1;
                            possible_poly = l;
                        }
                    }
                    if(full_possible_polys == 1) {
                        progress_made = true;
                        matches[k] = possible_poly;
                        for(int l = 0; l < nfactors; l++) {
                            if(l != possible_poly)
                                possible_root_table[k + full_dimension*l] = false;
                        }
                    }
                }

                // then we go though each polynomial and see if there are exactly the right
                // number of possible roots left
                for(int l = 0; l < nfactors; l++) {
                    if(roots_found[l] == fmpz_poly_degree(fmpz_factors + l)) continue;
                    vector<int> possible_roots;
                    for(int k = 0; k < full_dimension; k++) {
                        if(possible_root_table[k + full_dimension*l]) {
                            possible_roots.push_back(k);
                        }
                    }
                    if(possible_roots.size() == fmpz_poly_degree(fmpz_factors + l)) {
                        progress_made = true;
                        cout << "found all roots for polynomial " << l << endl;
                        for(int k : possible_roots) {
                            for(int l2 = 0; l2 < nfactors; l2++) {
                                if(l2 != l) possible_root_table[k + full_dimension*l2] = false;
                            }
                            matches[k] = l;
                        }
                        roots_found[l] = fmpz_poly_degree(g);
                    }
                    else {
                        cout << "ambiguity remains for polynomial " << l << endl;
                        matched_all_roots = false;
                    }
                }
            }
        }

        bool single_operator = false;
        if(unique_eigenvalues == 2) {
            single_operator = true;
        }
        else if(!matched_all_roots) {
            // trying to find single hecke operator with unique eigenvalue
            cout << "Trying harder to find a single hecke operator with squarefree characteristic polynomial." << endl;
            unique_eigenvalues = find_unique_eigenvalues(hecke_operator,
                                        eigenvalues,
                                        dimension,
                                        full_dimension,
                                        coeff_table,
                                        orbit,
                                        level,
                                        true, // require single operator
                                        100, // max operator
                                        prec);
            if(unique_eigenvalues != 0) {
                single_operator = true;
            }
        }

        long a = G.exponent(chi_number, unique_eigenvalues);
        int character_value_order = G.phi_q/GCD(a, G.phi_q);
        cout << character_value_order << endl;

        if(!matched_all_roots && single_operator) {
            cout << "single hecke operator, so attempting gcd trick." << endl;
            long gcd_prec = 4*(max_coefficient_bits + full_dimension + 100);
            while(!matched_all_roots && gcd_prec < 5000000) {
                cout << "trying gcd match with precsion " << gcd_prec << endl;
                matched_all_roots =  match_eigenvalues_rotate_gcd(matches,
                                         fmpz_factors,
                                         nfactors,
                                         eigenvalues,
                                         full_dimension,
                                         character_value_order,
                                         gcd_prec);
                gcd_prec = 2*gcd_prec;
            }

        }

        fmpz ** ztraces;
        vector<int> * mforbits;
        if(!matched_all_roots) {
            cout << "ohno we couldn't match all the eigenvalues to irreducible factors." << endl;
            goto cleanup;
        }

        cout << "checking traces to make sure that they are integers" << endl;
        ztraces = new fmpz*[nfactors];
        for(int l = 0; l < nfactors; l++) {
            //cout << l;
            arb_ptr traces = _arb_vec_init(100);
            for(int k = 0; k < full_dimension; k++) {
                if(matches[k] == l) {
                    int chi = orbit[k/dimension];
                    int chi_inverse = InvMod(chi, level);
                    if(level == 1) chi_inverse = 1;
                    int j = k % dimension;
                    //cout << " " << chi << "." << j;
                    acb_ptr coeffs;
                    if(chi <= chi_inverse) {
                        coeffs = coeff_table[make_pair(chi, j)];
                    }
                    else {
                        coeffs = coeff_table[make_pair(chi_inverse, j)];
                    }
                    for(int n = 0; n < 100; n++) {
                        arb_add(traces + n, traces + n, acb_realref(coeffs + n), prec);
                    }
                }
            }
            //fmpz_t t;
            //fmpz_init(t);
            ztraces[l] = _fmpz_vec_init(100);
            for(int n = 0; n < 100; n++) {
                if(!arb_get_unique_fmpz(ztraces[l] + n, traces + n)) {
                    cout << "ohno a trace was not an integer" << endl;
                    cout << "this is not supposed to happen (but maybe we are not using enough precision)" << endl;
                    return 1;
                }
                else {
                    cout << ztraces[l] + n << " ";
                }
            }
            cout << endl;
            _arb_vec_clear(traces, 100);
        }

        mforbits = new vector<int>[nfactors];
        for(int k = 0; k < full_dimension; k++) {
            if(matches[k] >= 0 && matches[k] < nfactors)
                mforbits[matches[k]].push_back(k);
        }

        // very simple sort on the traces and the factors...
        for(int l1 = 0; l1 < nfactors; l1++) {
            for(int l2 = l1 + 1; l2 < nfactors; l2++) {
                int z = fmpz_vec_cmp(ztraces[l1], ztraces[l2], 100);
                if(z == 0) {
                    cout << "ohno the first 100 traces match. continuing anyway" << endl;
                }
                else if (z > 0) {
                    //if(fmpz_vec_cmp(ztraces[l], ztraces
                    fmpz * tmp = ztraces[l2];
                    ztraces[l2] = ztraces[l1];
                    ztraces[l1] = tmp;
                    fmpz_poly_swap(fmpz_factors + l1, fmpz_factors + l2);
                    mforbits[l1].swap(mforbits[l2]);
                }
            }
        }

        for(int l = 0; l < nfactors; l++) {
            cout << l << ":";
            for(int k = 0; k < 100; k++) {
                cout << " " << ztraces[l] + k;
            }
            cout << endl;
        }

        // All went well.
        // Now to just record all the information we just computed...

        sqlite3_exec(polydb, "BEGIN TRANSACTION", NULL, NULL, NULL);
        for(int l = 0; l < nfactors; l++) {
            fmpz_poly_struct * g = fmpz_factors + l;
            //fmpz_poly_set_ZZX(g, factors[l].a);
            fmpz_poly_print_pretty(g, "x");
            vector<int> mforbit;
            for(int k : mforbits[l]) {
                mforbit.push_back(orbit[k/dimension]);
                mforbit.push_back(k % dimension);
            }
            //if(roots_found[l] == fmpz_poly_degree(g)) {
            //   for(int k = 0; k < full_dimension; k++) {
            //        if(matches[k] == l) {
                    //    cout << " " << orbit[k/dimension] << '.' << k % dimension;
            //            mforbit.push_back(orbit[k/dimension]);
            //            mforbit.push_back(k % dimension);
            //        }
            //    }
            //}
            polydb_insert(polydb, g, hecke_operator.data(), hecke_operator.size(), mforbit.data(), mforbit.size(),
                    level, weight, orbit[0], character_orbit_labels[orbit[0]], l, l, 100, ztraces[l]);
            cout << endl;
        }
        sqlite3_exec(polydb, "END TRANSACTION", NULL, NULL, NULL);

        for(int l = 0; l < nfactors; l++) {
            _fmpz_vec_clear(ztraces[l], 100);
        }
        delete [] ztraces;

cleanup:

        fmpz_poly_clear(g);

        fmpz_poly_clear(zheckepoly);
        acb_poly_clear(heckepoly);
        _acb_vec_clear(eigenvalues, orbit.size() * dimension);
        delete [] matches;
        for(int l = 0; l < nfactors; l++) {
            fmpz_poly_clear(fmpz_factors + l);
        }
        delete [] fmpz_factors;
        delete [] evaluation_bits;
        delete [] roots_found;
    }

    flint_cleanup();
    delete [] character_orbit_labels;
    cout << "main function exiting normally." << endl;
    return 0;
}
