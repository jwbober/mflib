#include "cuspforms_modp.h"
#include "classnumbers.h"

#include "flint/nmod_poly.h"
#include "flint/fmpz_poly.h"

#include "flint/NTL-interface.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <thread>

#include "NTL/ZZX.h"
#include "NTL/ZZXFactoring.h"

#include "ThreadPool/ThreadPool.h"

using namespace std;

const char * usage =
" level weight chi p [verbose]\n"
"\n"
"Compute the splitting of S_weight^new (level, chi) into Galois orbits over QQ.\n"
"The first line of the output is the level, weight, and then a list of\n"
"characters in this Galois orbit. The second line is an integer n which indicates\n"
"that we are going to print the characteristic polynomials of\n"
"T_2 + T_3 + T_4 + ... + Tn and the third line is a list of integers giving the\n"
"dimensions of the Galois orbits. (The degrees of the fields of coefficients.)\n"
"After this there is a list of polynomials, one per line, which are the minimal\n"
"polynomials of a_2 + a_3 + ... + a_n.\n"
"\n"
"This program works by computing these polynomials mod p for enough small primes\n"
"to recover the exact polynomial over the integers.\n";


void coefficient_bound(fmpz_t bound, int level, int weight, int dimension, int X) {
    // compute an upper bound for the log of the size of the largest coefficient
    // of the characteristic polynomial of T2 + T3 + T4 + ... + TX acting
    // on a space of given weight and level.

    const int prec = 10000; // whatever

    arb_t zbound;
    arb_init(zbound);

    fmpz_t kk;

    arb_t c, t, s;
    arb_init(c);
    arb_init(t);
    arb_init(s);

    fmpz_init(kk);

    for(int n = 2; n <= X; n++) {
        arb_set_ui(c, n);
        arb_pow_ui(c, c, weight - 1, prec);
        arb_sqrt(c, c, prec);
        arb_mul_ui(c, c, ndivisors(n), prec);
        arb_add_ui(c, c, 1, prec);

        arb_set_ui(t, dimension + 1);
        arb_div(t, t, c, prec);
        arb_floor(t, t, prec);
        if(!arb_get_unique_fmpz(kk, t)) {
            cout << "ohno" << endl;
            exit(-1);
        }

        long k = fmpz_get_si(kk);

        arb_sub_ui(c, c, 1, prec);

        arb_set_ui(s, dimension);
        arb_sub_ui(s, s, k, prec);
        arb_pow(c, c, s, prec);
        fmpz_bin_uiui(kk, dimension, k);
        arb_mul_fmpz(c, c, kk, prec);

        arb_add(zbound, zbound, c, prec);
    }

    arf_t Z;
    arf_init(Z);
    arb_get_abs_ubound_arf(Z, zbound, prec);
    arf_get_fmpz(bound, Z, ARF_RND_CEIL);

    arf_clear(Z);

    fmpz_clear(kk);
    arb_clear(zbound);
    arb_clear(c);
    arb_clear(t);
    arb_clear(s);
}


int main(int argc, char ** argv) {
    clock_t start_time = clock();
    int level;
    int chi_number;
    int weight;
    long p;

    if(argc < 5) {
        cout << argv[0] << usage;
        return 0;
    }


    init_classnumbers();

    level = atoi(argv[1]);
    weight = atoi(argv[2]);
    chi_number = atoi(argv[3]);
    int nthreads = atoi(argv[4]);

    p = 1125899906842679;

    //cout << level << endl
    //     << weight << endl
    //     << chi_number << endl
    //     << p << endl;

    int verbose = 0;
    if(argc > 5) verbose = atoi(argv[5]);
    int verbose2 = 0;
    if(verbose > 0) verbose2 = verbose - 1;

    DirichletGroup G(level);
    if(GCD(level, chi_number) != 1) return 0;
    DirichletCharacter chi = G.character(chi_number);
    if(chi.is_even() && weight % 2 == 1) return 0;
    if(!chi.is_even() && weight % 2 == 0) return 0;

    set<long> _orbit = chi.galois_orbit();
    vector<long> orbit(_orbit.begin(), _orbit.end());
    cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
    int dim = S->new_dimension();
    int orbitsize = orbit.size();
    int total_dimension = dim * orbitsize;
    cout << "dimension over cyclotomic field = " << dim << endl;
    cout << "dimension over Q = " << total_dimension << endl;
    if(dim == 0) return 0;
    p = S->p;

    fmpz_poly_t hecke1;
    fmpz_poly_t hecke2;

    fmpz_poly_init(hecke1);
    fmpz_poly_init(hecke2);

    fmpz_poly_one(hecke1);
    fmpz_poly_one(hecke2);

    nmod_poly_t hecke_poly_modp;

    fmpz_t a;
    fmpz_t b;
    fmpz_t modulus;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(modulus);
    fmpz_set_ui(modulus, p);

    // start by trying to see if we can find some n such that the characteristic
    // polynomial of T(2) + T(3) + ... + T(n) is squarefree. I don't actually know
    // if this will work well mod p, but I guess it will for large p at least, and
    // I suppose I'll find out quickly if it doesn't.

    nmod_mat_t * hecke_mats = new nmod_mat_t[orbit.size()];
    for(unsigned int k = 0; k < orbit.size(); k++) {
        nmod_mat_init(hecke_mats[k], dim, dim, p);
        nmod_mat_zero(hecke_mats[k]);
    }

    clock_t loopstart = clock();
    nmod_poly_init(hecke_poly_modp, p);
    int n = 1; {
        nmod_poly_t f;
        nmod_poly_init(f, p);

        do {
            n += 1;

            if(verbose) {
                cout << n << endl;
            }

            ThreadPool * pool = new ThreadPool(min(nthreads, orbitsize));

            nmod_poly_t * local_polys = new nmod_poly_t[orbitsize];
            for(int k = 0; k < orbitsize; k++) {
                nmod_poly_init(local_polys[k], p);
            }

            for(int k = 0; k < orbitsize;  k++) {
                long m = orbit[k];
                chi = G.character(m);
                cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
                p = S->p;
                pool->enqueue( [hecke_mats, local_polys, k, S, p, &n]() {
                    nmod_mat_t Tn;
                    S->hecke_matrix(Tn, n);
                    nmod_mat_add(hecke_mats[k], hecke_mats[k], Tn);
                    nmod_mat_charpoly(local_polys[k], hecke_mats[k]);
                    //nmod_poly_mul(hecke_poly_modp, hecke_poly_modp, f);
                    nmod_mat_clear(Tn);
                });
            }
            delete pool;

            nmod_poly_one(hecke_poly_modp);
            for(int k = 0; k < orbitsize; k++) {
                nmod_poly_mul(hecke_poly_modp, hecke_poly_modp, local_polys[k]);
                nmod_poly_clear(local_polys[k]);
            }
            delete [] local_polys;

            //for(unsigned int k = 0; k < orbit.size(); k++) {
            //    chi = G.character(orbit[k]);
            //    cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
            //    nmod_mat_t Tn;
            //    S->hecke_matrix(Tn, n);
            //    nmod_mat_add(hecke_mats[k], hecke_mats[k], Tn);
            //    nmod_mat_charpoly(f, hecke_mats[k]);
            //    nmod_poly_mul(hecke_poly_modp, hecke_poly_modp, f);
            //    nmod_mat_clear(Tn);
            // }
        } while(!nmod_poly_is_squarefree(hecke_poly_modp));

        //nmod_poly_factor_t factorization;
        //nmod_poly_factor_init(factorization);
        //nmod_poly_factor(factorization, hecke_poly_modp);
        //cout << "factorization type: ";
        //for(int k = 0; k < factorization[0].num; k++) {
        //    cout << " (" << nmod_poly_degree(factorization[0].p + k)
        //         << "," << factorization[0].exp[k] << ")";
        //}
        //cout << endl;
        //nmod_poly_factor_clear(factorization);

        int degree = nmod_poly_degree(hecke_poly_modp);
        for(int k = 0; k <= degree; k++) {
            long z = nmod_poly_get_coeff_ui(hecke_poly_modp, k);
            if(z > p/2) {
                z -= p;
            }
            fmpz_poly_set_coeff_si(hecke2, k, z);
        }
        nmod_poly_clear(f);
    }
    for(unsigned int k = 0; k < orbit.size(); k++) {
        nmod_mat_clear(hecke_mats[k]);
    }

    delete [] hecke_mats;

    fmpz_t bound;
    fmpz_init(bound);

    coefficient_bound(bound, level, weight, total_dimension, n);

    nmod_poly_clear(hecke_poly_modp);
    if(verbose) {
        cout << "finished initial work." << endl;
    }

    do {
        //fmpz_poly_print_pretty(hecke2, "x"); cout << endl;
        fmpz_poly_set(hecke1, hecke2);
        p += 1;
        cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
        p = S->p;
        //if(verbose) cout << p << endl;
        if(verbose) {
            clock_t looptime = clock();
            double time_so_far = (looptime - loopstart)/(double)CLOCKS_PER_SEC;
            double percent_finished = fmpz_bits(modulus)/(double)fmpz_bits(bound);
            double total_estimate = time_so_far/percent_finished;
            cout << percent_finished << " " << time_so_far << " " << total_estimate << endl;
        }
        nmod_poly_init(hecke_poly_modp, p);
        nmod_poly_one(hecke_poly_modp);
        nmod_poly_t f;
        nmod_poly_init(f, p);
        nmod_poly_t * local_polys = new nmod_poly_t[orbitsize];
        for(int k = 0; k < orbitsize; k++) {
            nmod_poly_init(local_polys[k], p);
        }

        ThreadPool * pool = new ThreadPool(min(nthreads, orbitsize));

        //for(long m : orbit) {
        for(int k = 0; k < orbitsize;  k++) {
            long m = orbit[k];
            chi = G.character(m);
            cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
            p = S->p;
            pool->enqueue( [local_polys, k, S, p, &n]() {
                nmod_mat_t hecke_mat;
                S->hecke_matrix(hecke_mat, 2);
                for(int l = 3; l <= n; l++) {
                    nmod_mat_t Tl;
                    S->hecke_matrix(Tl, l);
                    nmod_mat_add(hecke_mat, hecke_mat, Tl);
                    nmod_mat_clear(Tl);
                }
                nmod_mat_charpoly(local_polys[k], hecke_mat);
                //nmod_poly_mul(hecke_poly_modp, hecke_poly_modp, f);
                nmod_mat_clear(hecke_mat);
            });
        }
        delete pool;
        for(int k = 0; k < orbitsize; k++) {
            nmod_poly_mul(hecke_poly_modp, hecke_poly_modp, local_polys[k]);
            nmod_poly_clear(local_polys[k]);
        }

        delete [] local_polys;

        //nmod_poly_factor_t factorization;
        //nmod_poly_factor_init(factorization);
        //nmod_poly_factor(factorization, hecke_poly_modp);
        //cout << "factorization type: ";
        //for(int k = 0; k < factorization[0].num; k++) {
        //    cout << " (" << nmod_poly_degree(factorization[0].p + k)
        //         << "," << factorization[0].exp[k] << ")";
        //}
        //cout << endl;
        //nmod_poly_factor_clear(factorization);


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
    } while(fmpz_cmp(modulus, bound) < 0);

    //while(!fmpz_poly_equal(hecke1, hecke2));

    {
        fmpz_t x;
        fmpz_t largest;
        fmpz_init(x); fmpz_init(largest);

        for(int k = 0; k < fmpz_poly_degree(hecke1); k++) {
            fmpz_poly_get_coeff_fmpz(x, hecke1, k);
            fmpz_abs(x, x);
            if(fmpz_cmp(largest, x) < 0) {
                fmpz_set(largest, x);
            }
        }

        cout << "largest coefficient: ";
        fmpz_print(largest);
        cout << endl;

        cout << "upper bound: ";
        fmpz_print(bound);
        cout << endl;

        fmpz_clear(x); fmpz_clear(largest);
    }

    if(!fmpz_poly_is_squarefree(hecke1)) {
        cout << "something went wrong." << endl;
        cout << n << " ";
        fmpz_poly_print_pretty(hecke1, "x");
        cout << endl;
        return 1;
    }
    cout << level << " " << weight;
    for(long m : orbit) {
        cout << " " << m;
    }
    cout << endl;
    cout << n << endl;
    NTL::ZZ content;
    NTL::ZZX ff;
    NTL::vec_pair_ZZX_long factors;
    fmpz_poly_get_ZZX(ff, hecke1);
    NTL::factor(content, factors, ff);
    fmpz_poly_t g;
    fmpz_poly_init(g);
    for(int l = 0; l < factors.length(); l++) {
        cout << deg(factors[l].a);
        if(l < factors.length() - 1) cout << " ";
    }
    cout << endl;
    for(int l = 0; l < factors.length(); l++) {
        fmpz_poly_set_ZZX(g, factors[l].a);
        fmpz_poly_print_pretty(g, "x");
        cout << endl;
    }

    clock_t end_time = clock();
    //cout << (end_time - start_time)/(double)CLOCKS_PER_SEC << endl;

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(modulus);

    fmpz_poly_clear(hecke1);
    fmpz_poly_clear(hecke2);
    flint_cleanup();
    return 0;
}
