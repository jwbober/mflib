#include "modform_modp.h"
#include "classnumbers.h"

#include "flint/nmod_poly.h"
#include "flint/fmpz_poly.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

int main(int argc, char ** argv) {
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

    bool finished = false;
    while(!finished) {
        nmod_poly_init(f, p);
        nmod_poly_one(f);
        for(long m : orbit) {
            chi = G.character(m);
            //cout << "computing space S_" << weight << "(" << level << "," << m << ")" << endl;
            cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose);
            p = S->p;
            nmod_mat_t Tl;
            S->hecke_matrix(Tl, 2);
            //cout << "computed the space." << endl;
            //nmod_mat_print_pretty(Tl);
            nmod_poly_t hecke_poly;
            nmod_poly_init(hecke_poly, p);
            nmod_mat_charpoly(hecke_poly, Tl);
            //cout << p << " " << m << " ";
            //nmod_poly_print_pretty(hecke_poly, "x");
            //cout << endl;
            nmod_poly_mul(f, f, hecke_poly);
            nmod_poly_clear(hecke_poly);
            nmod_mat_clear(Tl);
        }
        if(fmpz_cmp_ui(modulus, 0) != 0) {
            int degree = nmod_poly_degree(f);
            cout << degree << endl;
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
            cout << "a = ";
            fmpz_print(a);
            cout << endl;
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

    fmpz_poly_print_pretty(hecke1, "x");
    cout << endl;

    return 0;
}
