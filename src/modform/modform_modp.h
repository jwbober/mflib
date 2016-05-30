#ifndef __MODFORM_MODP_H__
#define __MODFORM_MODP_H__

#include <vector>

#include "flint/flint.h"
#include "flint/nmod_mat.h"

#include "slint.h"
#include "characters.h"

static long psi(int N) {
    int_factorization_t f;
    factor(N, f);

    long PSI = 1;
    for(int k = 0; k < f.nfactors; k++) {
        PSI *= f.factors[k].f/f.factors[k].p * (long)(f.factors[k].p + 1);
    }
    return PSI;
}

class cuspforms_modp;
cuspforms_modp * get_cuspforms_modp(DirichletCharacter &chi, int weight, long p, int verbose = 0);

class cuspforms_modp {
public:
    nmod_t modp;
    long p;
    int level;
    int weight;
    int chi;
    int conductor;
    int verbose;

    std::vector<int> sublevels;
    std::vector<int> divisors_of_level;
    std::vector<cuspforms_modp*> subspaces;
    std::vector<DirichletGroup> sublevel_dirichlet_groups;

    int * divisor_counts;
    long * psi_table;
    int * phi_table;
    int ** gcd_tables;
    long * chi_values;
    long * chip_values;

    std::vector<int> basis_rows;

    std::vector<long> traces;

    long trace(int n);
    long trace_TnTm(int n, int m);
    long trace_TpTnTm(int p, int n, int m);
            // p is getting confusing.
            // I suppose I should be computing cuspforms mod l, because
            // I really don't want to refer to the lth Fourier coeffcient...

    void compute_traces(int end);

    int dimension();
    int new_dimension();

    void newspace_basis(nmod_mat_t basis, int ncoeffs);
    void newforms(nmod_mat_t forms, int ncoeffs); // I'm not sure what I was
                                                  // planning for this function to do.

    void hecke_matrix(nmod_mat_t Tp, int n);

    const std::vector<int>& newspace_basis_data();
    void basis(nmod_mat_t * basis, int ncoeffs);

    long evalpoly(long t, long n);


    cuspforms_modp(DirichletCharacter &_chi, const int w, long _p, int _verbose = 0) {
        //if(w != 2) {
        //    std::cerr << "only weight 2 is supported at the moment..." << std::endl;
        //    exit(1);
        //}
        verbose = _verbose;
        chi = _chi.m;
        level = _chi.parent->q;
        long primitive_index;
        conductor = _chi.conductor(&primitive_index);
        weight = w;
        p = _p;

        int order = order_mod(_chi.m, level);
        if(p == 0) {
            p = order + 1;
            while(p < 10000) p += order;
        }
        else {
            if( (p - 1) % order != 0) {
                p += order - (p - 1) % order;
            }
        }
        while(!is_prime(p)) {p += order;}

        std::cout << level << " " << chi << " " << p << std::endl;
        nmod_init(&modp, p);

        psi_table = new long[level + 1];
        phi_table = new int[level + 1];
        gcd_tables = new int*[level + 1];

        for(int k = 1; k < level + 1; k++) {
            if(level % k == 0) {
                psi_table[k] = psi(k);
                gcd_tables[k] = new int[k];
                for(int j = 0; j < k; j++) {
                    gcd_tables[k][j] = GCD(k,j);
                }
                phi_table[k] = euler_phi(k);
            }
        }

        chi_values = new long[level + 1];
        divisor_counts = new int[level + 1];
        _chi.values_mod_p(p, chi_values);

        DirichletGroup G(conductor);
        DirichletCharacter chip = G.character(primitive_index);
        chip_values = new long[conductor + 1];
        chip.values_mod_p(p, chip_values);

        divisors_of_level = divisors(level);
        for(int M : divisors_of_level) {
            divisor_counts[M] = ndivisors(M);
            if(M % conductor == 0 && M != level) {
                sublevels.push_back(M);
                DirichletGroup G(M);
                DirichletCharacter psi = G.character(G.index_from_primitive_character(conductor, primitive_index));
                cuspforms_modp * S2 = get_cuspforms_modp(psi, weight, p, verbose);
                subspaces.push_back(S2);
            }
        }
    }

    ~cuspforms_modp() {
        delete [] psi_table;
        for(int k = 1; k < level + 1; k++) {
            if(level % k == 0) {
                delete [] gcd_tables[k];
            }
        }
        delete [] gcd_tables;
        delete [] phi_table;
        delete [] chi_values;
        delete [] divisor_counts;
        delete [] chip_values;
    }
};

#endif
