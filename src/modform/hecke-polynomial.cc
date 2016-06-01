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

    DirichletGroup G(level);
    if(GCD(level, chi_number) != 1) return 0;
    DirichletCharacter chi = G.character(chi_number);
    if(chi.is_even() && weight % 2 == 1) return 0;
    if(!chi.is_even() && weight % 2 == 0) return 0;

    set<long> orbit = chi.galois_orbit();
    cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose);
    int dim = S->new_dimension();
    if(dim == 0) return 0;
    p = S->p;

    fmpz_poly_t hecke1;
    fmpz_poly_t hecke2;

    fmpz_poly_init(hecke1);
    fmpz_poly_init(hecke2);

    fmpz_poly_one(hecke1);
    fmpz_poly_one(hecke2);

    nmod_poly_t f;

    fmpz_t a;
    fmpz_t b;
    fmpz_t modulus;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(modulus);
    fmpz_zero(modulus);

    // start by trying to see if we can find some n such that the characteristic
    // polynomial of T(2) + T(3) + ... + T(n) is squarefree. I don't actually know
    // if this will work well mod p, but I guess it will for large p at least, and
    // I suppose I'll find out quickly if it doesn't.

    int n = 2; {
        nmod_mat_t hecke_mat;
        nmod_poly_t hecke_poly;
        nmod_poly_init(hecke_poly, p);
        S->hecke_matrix(hecke_mat, 2);
        nmod_mat_charpoly(hecke_poly, hecke_mat);
        while(!nmod_poly_is_squarefree(hecke_poly)) {
            n += 1;
            nmod_mat_t Tn;
            S->hecke_matrix(Tn, n);
            nmod_mat_add(hecke_mat, hecke_mat, Tn);
            nmod_mat_charpoly(hecke_poly, hecke_mat);
            nmod_mat_clear(Tn);
        }
        nmod_mat_clear(hecke_mat);
        nmod_poly_clear(hecke_poly);
    }
    bool finished = false;
    while(!finished) {
        nmod_poly_init(f, p);
        nmod_poly_one(f);
        for(long m : orbit) {
            chi = G.character(m);
            //cout << "computing space S_" << weight << "(" << level << "," << m << ")" << endl;
            cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose);
            p = S->p;
            nmod_mat_t hecke_mat;
            S->hecke_matrix(hecke_mat, 2);
            for(int l = 3; l <= n; l++) {
                nmod_mat_t Tl;
                S->hecke_matrix(Tl, l);
                nmod_mat_add(hecke_mat, hecke_mat, Tl);
                nmod_mat_clear(Tl);
            }
            //cout << "computed the space." << endl;
            //nmod_mat_print_pretty(Tl);
            nmod_poly_t hecke_poly;
            nmod_poly_init(hecke_poly, p);
            nmod_mat_charpoly(hecke_poly, hecke_mat);
            //cout << p << " " << m << " ";
            //nmod_poly_print_pretty(hecke_poly, "x");
            //cout << endl;
            nmod_poly_mul(f, f, hecke_poly);
            nmod_poly_clear(hecke_poly);
            nmod_mat_clear(hecke_mat);
        }
        if(fmpz_cmp_ui(modulus, 0) != 0) {
            int degree = nmod_poly_degree(f);
            for(int k = 0; k <= degree; k++) {
                fmpz_poly_get_coeff_fmpz(a, hecke1, k);
                unsigned long z = nmod_poly_get_coeff_ui(f, k);
                fmpz_CRT_ui(a, a, modulus, z, p, 1);
                fmpz_poly_set_coeff_fmpz(hecke2, k, a);
            }

            if(fmpz_poly_equal(hecke1, hecke2)) {
                finished = true;
            }
            else {
                fmpz_mul_ui(modulus, modulus, p);
                fmpz_poly_set(hecke1, hecke2);
            }
            fmpz_poly_get_coeff_fmpz(a, hecke1, 0);
        }
        else {
            fmpz_set_ui(modulus, p);
            int degree = nmod_poly_degree(f);
            for(int k = 0; k <= degree; k++) {
                unsigned long z = nmod_poly_get_coeff_ui(f, k);
                fmpz_poly_set_coeff_ui(hecke1, k, z);
            }
        }
        p += 1;
        cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose);
        p = S->p;
        nmod_poly_clear(f);
    }

    if(!fmpz_poly_is_squarefree(hecke1)) {
        cout << "something went wrong." << endl;
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

    return 0;
}
