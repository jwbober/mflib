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

#include "NTL/ZZX.h"
#include "NTL/ZZXFactoring.h"

#include "acb_poly.h"

#include "arb-extras.h"

using namespace std;

const char * usage =
" leve weight mfdb polydb [overwrite]\n"
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

ostream& operator << (ostream& out, fmpz_t in) {
    char * outstr = fmpz_get_str(NULL, 10, in);
    out.write(outstr, strlen(outstr));
    flint_free(outstr);
    return out;
}


void arb_get_ubound_fmpz(fmpz_t out, arb_t in, slong prec) {
    arf_t t;
    arf_init(t);
    arb_get_ubound_arf(t, in, prec);
    arf_get_fmpz(out, t, ARF_RND_UP);
    arf_clear(t);
}

int arb_get_unique_fmpz_modular(fmpz_t out, arb_t real_approx, fmpz_t modular_approx, fmpz_t mod) {
    // compute and return the unique integer in the interval represented by real_approx
    // and congruent to modular_approx modulo mod, if such an integer exists.
    //
    // Returns 1 on success, 0 if there is no such integer.

    fmpz_t a, b, e;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(e);

    arb_get_interval_fmpz_2exp(a, b, e, real_approx);

    // we really should be more careful with error checking here.
    slong exponent = fmpz_get_si(e);
    if(exponent > 0)
        fmpz_mul_2exp(a, a, exponent);
    else
        fmpz_fdiv_q_2exp(a, a, -exponent);
    fmpz_mod(b, a, mod);
    fmpz_sub(a, a, b);
    fmpz_add(a, a, modular_approx);

    // We are probably in the interval now, if there is such an integer
    // in the interval. If I were to think properly about what I'm doing,
    // then I could just check a + mod probably to make sure that it isn't
    // also in the interval. But there is rounding and stuff to worry about...

    int nfound = 0;
    fmpz_sub(a, a, mod);
    fmpz_sub(a, a, mod);

    if(arb_contains_fmpz(real_approx, a)) {
        nfound += 1;
        fmpz_set(out, a);
    }
    for(int k = 0; k < 4; k++) {
        fmpz_add(a, a, mod);
        if(arb_contains_fmpz(real_approx, a)) {
            nfound += 1;
            fmpz_set(out, a);
        }
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(e);

    if(nfound != 1) return 0;

    return 1;
}

int arb_poly_get_unique_fmpz_poly_modular(fmpz_poly_t out, arb_poly_t real_approx, fmpz_poly_t modular_approx, fmpz_t mod) {
    fmpz_t t;
    fmpz_init(t);
    slong degree = arb_poly_degree(real_approx);
    for(slong k = 0; k <= degree; k++) {
        int result = arb_get_unique_fmpz_modular(
                t,
                arb_poly_get_coeff_ptr(real_approx, k),
                fmpz_poly_get_coeff_ptr(modular_approx, k),
                mod);
        if(result == 0) {
            fmpz_clear(t);
            return 0;
        }
        fmpz_poly_set_coeff_fmpz(out, k, t);
    }
    fmpz_clear(t);
    return 1;
}

int acb_poly_get_unique_fmpz_poly_modular(fmpz_poly_t out, acb_poly_t real_approx, fmpz_poly_t modular_approx, fmpz_t mod) {
    fmpz_t t;
    fmpz_init(t);
    slong degree = acb_poly_degree(real_approx);
    for(slong k = 0; k <= degree; k++) {
        if(!arb_contains_zero(acb_imagref(acb_poly_get_coeff_ptr(real_approx, k)))) {
            fmpz_clear(t);
            return 0;
        }
        int result = arb_get_unique_fmpz_modular(
                t,
                acb_realref(acb_poly_get_coeff_ptr(real_approx, k)),
                fmpz_poly_get_coeff_ptr(modular_approx, k),
                mod);
        if(result == 0) {
            fmpz_clear(t);
            return 0;
        }
        fmpz_poly_set_coeff_fmpz(out, k, t);
    }
    fmpz_clear(t);
    return 1;
}

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
        nmod_poly_t f;
        nmod_poly_init(f, p);
        for(long m : orbit) {
            chi = G.character(m);
            cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
            p = S->p;
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
            nmod_mat_charpoly(f, hecke_mat);
            nmod_poly_mul(hecke_poly_modp, hecke_poly_modp, f);
            nmod_mat_clear(hecke_mat);
        }
        int degree = nmod_poly_degree(hecke_poly_modp);
        for(int k = 0; k <= degree; k++) {
            fmpz_poly_get_coeff_fmpz(a, hecke1, k);
            unsigned long z = nmod_poly_get_coeff_ui(hecke_poly_modp, k);
            fmpz_CRT_ui(a, a, modulus, z, p, 1);
            fmpz_poly_set_coeff_fmpz(hecke2, k, a);
        }
        fmpz_mul_ui(modulus, modulus, p);

        nmod_poly_clear(f);
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
    string copy_stmt = "INSERT INTO modforms SELECT * from infiledb.modforms WHERE level=" + to_string(level) + " AND weight=" + to_string(weight);
    retval = sqlite3_exec(db, copy_stmt.c_str(), NULL, 0, NULL);
    retval = sqlite3_exec(db, "DETACH DATABASE infiledb;", NULL, 0, NULL);

    cout << "input data loaded into memory." << endl;

    sqlite3 * polydb;
    if(sqlite3_open(polydbname.c_str(), &polydb) != SQLITE_OK) {
        cout << "ohno error opening output database." << endl;
        return 1;
    }
    polydb_init(polydb);

    int overwrite = 0;
    if(argc > 6) {
        overwrite = atoi(argv[6]);
    }
    if(overwrite) {
        sqlite3_exec(polydb, "DELETE FROM heckepolys;", NULL, 0, NULL);
    }

    mfheader * headers;
    int count = mfdb_contents(db, &headers);

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
        weight = headers[k].weight;
        chi_list.insert(headers[k].chi);
        if(-headers[k].prec > prec) prec = -headers[k].prec;
        if(headers[k].j + 1 > dimensions[headers[k].chi]) dimensions[headers[k].chi] = headers[k].j+1;
        //cout << headers[k].level << " " << headers[k].weight << " " << headers[k].chi << " " << headers[k].j << endl;
    }
    if(prec < 0) prec = 200;
    prec = 30*prec;
    cout << "using working precision " << prec << endl;
    free(headers);

    DirichletGroup G(level);
    while(chi_list.size() > 0) {
        long chi_number = *chi_list.begin();
        DirichletCharacter chi = G.character(chi_number);
        int realchi = false;
        if(order_mod(chi_number, level) < 3) realchi = true;
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
        prec = full_dimension * 300;
        acb_ptr eigenvalues = _acb_vec_init(full_dimension);
        bool unique_eigenvalues = false;
        int k = 2;

        while(!unique_eigenvalues && (k < 30)) {
            // we really prefer to compute the minimal polynomial of just one hecke-eigenvalue,
            // so let's try that first.
            cout << "Trying Hecke operator " << k << endl;
            hecke_operator.push_back(1);
            for(int l = 0; l < orbit.size(); l++) {
                long chi_number = orbit[l];
                long chi_inverse = InvMod(chi_number, level);
                if(level == 1) chi_inverse = 1;
                for(int j = 0; j < dimension; j++) {
                    if(chi_number <= chi_inverse) {
                        int result = mfdb_get_entry(db, &header, &coeffs, level, weight, chi_number, j);
                        acb_set(eigenvalues + l * dimension + j, coeffs + k - 1);
                    }
                    else {
                        int result = mfdb_get_entry(db, &header, &coeffs, level, weight, chi_inverse, j);
                        acb_conj(coeffs + k - 1, coeffs + k - 1);
                        acb_set(eigenvalues + l * dimension + j, coeffs + k - 1);
                    }
                    _acb_vec_clear(coeffs, header.ncoeffs);
                }
            }
            unique_eigenvalues = true;
            for(int n = 0; n < dimension * orbit.size(); n++) {
                for(int m = n + 1; m < dimension * orbit.size(); m++) {
                    if(acb_overlaps(eigenvalues + n, eigenvalues + m)) {
                        unique_eigenvalues = false;
                        break;
                    }
                }
                if(!unique_eigenvalues) {
                    hecke_operator.pop_back();
                    hecke_operator.push_back(0);
                    break;
                }
            }
            k++;

        }

        if(!unique_eigenvalues) {
            hecke_operator.clear();
            k = 2;
        }

        while(!unique_eigenvalues) {
            hecke_operator.push_back(k - 1);
            for(int l = 0; l < orbit.size(); l++) {
                long chi_number = orbit[l];
                long chi_inverse = InvMod(chi_number, level);
                if(level == 1) chi_inverse = 1;
                for(int j = 0; j < dimension; j++) {
                    if(chi_number <= chi_inverse) {
                        int result = mfdb_get_entry(db, &header, &coeffs, level, weight, chi_number, j);
                        acb_addmul_ui(eigenvalues + l * dimension + j, coeffs + k - 1, k - 1, prec);
                    }
                    else {
                        int result = mfdb_get_entry(db, &header, &coeffs, level, weight, chi_inverse, j);
                        acb_conj(coeffs + k - 1, coeffs + k - 1);
                        acb_addmul_ui(eigenvalues + l * dimension + j, coeffs + k - 1, k - 1, prec);
                    }
                    _acb_vec_clear(coeffs, header.ncoeffs);
                }
            }
            unique_eigenvalues = true;
            for(int n = 0; n < dimension * orbit.size(); n++) {
                for(int m = n + 1; m < dimension * orbit.size(); m++) {
                    if(acb_overlaps(eigenvalues + n, eigenvalues + m)) {
                        unique_eigenvalues = false;
                        break;
                    }
                }
                if(!unique_eigenvalues) break;
            }
            k++;
        }
        cout << "for character " << chi_number << " using hecke operator";
        for(int k : hecke_operator) cout << " " << k;
        cout << endl;

        acb_poly_t heckepoly;
        acb_poly_init(heckepoly);
        acb_poly_product_roots(heckepoly, eigenvalues, full_dimension, prec);
        //for(int l = 0; l < full_dimension; l++) {
        //    acb_poly_set_coeff_si(f, 1, 1);
        //    acb_mul_si(eigenvalues + l, eigenvalues + l, -1, prec);
        //    acb_poly_set_coeff_acb(f, 0, eigenvalues + l);
        //    acb_mul_si(eigenvalues + l, eigenvalues + l, -1, prec);
        //    acb_poly_mul(heckepoly, heckepoly, f, prec);
        //}
        //acb_poly_printd(heckepoly, 10);
        //cout << endl;

        fmpz_poly_t zheckepoly;
        fmpz_poly_init(zheckepoly);
        for(int k = 0; k <= full_dimension; k++) {
            arb_zero(acb_imagref(acb_poly_get_coeff_ptr(heckepoly, k)));
        }
        if(!acb_poly_get_unique_fmpz_poly(zheckepoly, heckepoly)) {
            // we're going to have to up the accuracy with a mod p computation
            cout << "doing mod p computations to increase hecke polynomial precision" << endl;
            arb_t max_error, r;
            arb_init(max_error);
            arb_init(r);
            arb_zero(max_error);
            for(int k = 0; k < full_dimension; k++) {
                arb_get_rad_arb(r, acb_realref(acb_poly_get_coeff_ptr(heckepoly, k)));
                //cout << r << " " << acb_realref(acb_poly_get_coeff_ptr(heckepoly, k)) << endl;
                //arb_printd(r, 10);
                //cout << endl;
                //arb_printd(acb_realref(acb_poly_get_coeff_ptr(heckepoly, k)), 10);
                //cout << endl;
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
            if(result == 0) {
                cout << "ohno" << endl;
                continue;
            }
            result = acb_poly_get_unique_fmpz_poly_modular(zheckepoly, heckepoly, heckepoly_modular_approx, modular_precision);
            if(result == 0) {
                cout << "ohno. something went wrong with modular approximation." << endl;
                continue;
            }
        }

        NTL::ZZ content;
        NTL::ZZX ff;
        NTL::vec_pair_ZZX_long factors;
        fmpz_poly_get_ZZX(ff, zheckepoly);
        NTL::factor(content, factors, ff);
        fmpz_poly_t g;
        fmpz_poly_init(g);

        int nfactors = factors.length();
        arb_poly_t arbfactor;
        arb_poly_init(arbfactor);
        int * matches = new int[full_dimension];
        for(int k = 0; k < full_dimension; k++) {
            matches[k] = -1;
        }

        acb_t t1;
        acb_init(t1);
        //for(int l = 0; l < factors.length(); l++) {
        //    cout << deg(factors[l].a);
        //    if(l < factors.length() - 1) cout << " ";
        //}
        //cout << endl;

        arb_poly_init(arbfactor);
        for(int l = 0; l < nfactors; l++) {
            fmpz_poly_set_ZZX(g, factors[l].a);
            long maxbits = 0;
            for(int k = 0; k <= fmpz_poly_degree(g); k++) {
                long bits = fmpz_bits(fmpz_poly_get_coeff_ptr(g, k));
                if(bits > maxbits) maxbits = bits;
            }
            arb_poly_set_fmpz_poly(arbfactor, g, maxbits + 10);
            long evaluation_bits = (maxbits + 10) * fmpz_poly_degree(g);
            if(evaluation_bits < prec) evaluation_bits = prec;
            for(int k = 0; k < full_dimension; k++) {
                arb_poly_evaluate_acb(t1, arbfactor, eigenvalues + k, evaluation_bits);
                if(acb_contains_zero(t1)) {
                    if(matches[k] != -1) {
                        cout << "ohno1" << endl;
                        continue;
                    }
                    matches[k] = l;
                }
            }
        }
        for(int k = 0; k < full_dimension; k++) {
            if(matches[k] == -1) cout << "ohno2" << endl;
        }

        // All went well.
        // Now to just record all the information we just computed...

        sqlite3_exec(polydb, "BEGIN TRANSACTION", NULL, NULL, NULL);
        for(int l = 0; l < nfactors; l++) {
            fmpz_poly_set_ZZX(g, factors[l].a);
            fmpz_poly_print_pretty(g, "x");
            vector<int> mforbit;
            for(int k = 0; k < full_dimension; k++) {
                if(matches[k] == l) {
                //    cout << " " << orbit[k/dimension] << '.' << k % dimension;
                    mforbit.push_back(orbit[k/dimension]);
                    mforbit.push_back(k % dimension);
                }
            }
            polydb_insert(polydb, g, hecke_operator.data(), hecke_operator.size(), mforbit.data(), mforbit.size(),
                    level, weight, orbit[0], l, 0);
            cout << endl;
        }
        sqlite3_exec(polydb, "END TRANSACTION", NULL, NULL, NULL);



        fmpz_poly_clear(g);
        //cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
        //int dim = S->new_dimension();
        //p = S->p;

        //cout << chi << " " << dimensions[chi] << endl;
        // now that we know what we are dealing with we do some actual computations...
        fmpz_poly_clear(zheckepoly);
        acb_poly_clear(heckepoly);
        _acb_vec_clear(eigenvalues, orbit.size() * dimension);
    }

    flint_cleanup();
    return 0;
}
