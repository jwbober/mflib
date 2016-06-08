#include "modform_modp.h"
#include "classnumbers.h"

#include "flint/nmod_poly.h"
#include "flint/fmpz_poly.h"

#include "flint/NTL-interface.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

#include "NTL/ZZX.h"
#include "NTL/ZZXFactoring.h"


using namespace std;

int main(int argc, char ** argv) {
    clock_t start_time = clock();
    int level;
    int chi_number;
    int weight;
    long p;

    init_classnumbers();

    level = atoi(argv[1]);
    weight = atoi(argv[2]);
    chi_number = atoi(argv[3]);
    p = atol(argv[4]);

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

    nmod_mat_t hecke_mats[orbit.size()];
    for(unsigned int k = 0; k < orbit.size(); k++) {
        nmod_mat_init(hecke_mats[k], dim, dim, p);
        nmod_mat_zero(hecke_mats[k]);
    }

    nmod_poly_init(hecke_poly_modp, p);
    int n = 1; {
        nmod_poly_t f;
        nmod_poly_init(f, p);

        do {
            n += 1;
            if(verbose) {
                cout << n << endl;
            }
            nmod_poly_one(hecke_poly_modp);
            for(unsigned int k = 0; k < orbit.size(); k++) {
                chi = G.character(orbit[k]);
                cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
                nmod_mat_t Tn;
                S->hecke_matrix(Tn, n);
                nmod_mat_add(hecke_mats[k], hecke_mats[k], Tn);
                nmod_mat_charpoly(f, hecke_mats[k]);
                nmod_poly_mul(hecke_poly_modp, hecke_poly_modp, f);
                nmod_mat_clear(Tn);
            }
        } while(!nmod_poly_is_squarefree(hecke_poly_modp));

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

    nmod_poly_clear(hecke_poly_modp);
    if(verbose) {
        cout << "finished initial work." << endl;
    }

    do {
        fmpz_poly_set(hecke1, hecke2);
        p += 1;
        cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
        p = S->p;
        if(verbose) cout << p << endl;
        nmod_poly_init(hecke_poly_modp, p);
        nmod_poly_one(hecke_poly_modp);
        nmod_poly_t f;
        nmod_poly_init(f, p);
        for(long m : orbit) {
            chi = G.character(m);
            cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose2);
            p = S->p;
            nmod_mat_t hecke_mat;
            S->hecke_matrix(hecke_mat, 2);
            for(int l = 3; l <= n; l++) {
                nmod_mat_t Tl;
                S->hecke_matrix(Tl, l);
                nmod_mat_add(hecke_mat, hecke_mat, Tl);
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
    } while(!fmpz_poly_equal(hecke1, hecke2));

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
    cout << (end_time - start_time)/(double)CLOCKS_PER_SEC << endl;

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(modulus);

    fmpz_poly_clear(hecke1);
    fmpz_poly_clear(hecke2);
    flint_cleanup();
    return 0;
}
