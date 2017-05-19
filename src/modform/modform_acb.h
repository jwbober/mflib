#ifndef __MODFORM_ACB_H__
#define __MODFORM_ACB_H__

// to suppress warnings about redefinition.
// (USE_ARB could be better named, like, _CHARACTERS_H_USE_ARB_...)
#ifdef USE_ARB
#undef USE_ARB
#endif

#define USE_ARB

#include <vector>
#include <complex>

#include "slint.h"

#include "characters.h"
#include "modform_modp.h"
#include "acb_mat.h"

//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>

//typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cmatrix_t;
//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rmatrix_t;

class cuspforms_acb;
cuspforms_acb * get_cuspforms_acb(DirichletCharacter &chi, int weight, int verbose = 0);
void clear_cuspform_cache();

class cuspforms_acb {
public:
    int level;
    int weight;
    int chi;
    int conductor;
    int verbose;

    std::vector<int> sublevels;
    std::vector<int> divisors_of_level;
    std::vector<cuspforms_acb*> subspaces;

    int * divisor_counts;
    int * psi_table;
    int * phi_table;
    int ** gcd_tables;
    int prec;

    acb_t * chi_values;
    acb_t * chip_values;

    std::vector<int> basis_cols;

    acb_t * traces;
    int traces_size = 0;
    int traces_computed = 0;

    //std::vector<std::complex<double> > traces;

    void trace(acb_t out, int n);
    void trace_TnTm(acb_t out, int n, int m);
    void trace_TpTnTm(acb_t out, int p, int n, int m);
    void compute_traces(int end);

    int dimension();
    long new_dimension();

    void newspace_basis(acb_mat_t B, int ncoeffs);
    const std::vector<int>& newspace_basis_data();
    //cmatrix_t basis(int ncoeffs);

    cuspforms_modp * modp_space;
    void newforms(acb_mat_t out, int ncoeffs);

    void evalpoly(fmpz_t out, fmpz_t t, fmpz_t n);

    cuspforms_acb(DirichletCharacter &_chi, const int w, int _verbose = 0) {
        verbose = _verbose;
        chi = _chi.m;
        level = _chi.parent->q;
        prec = _chi.parent->prec;
        long primitive_index;
        conductor = _chi.conductor(&primitive_index);
        weight = w;

        psi_table = new int[level + 1];
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

        modp_space = get_cuspforms_modp(_chi, weight, 10000000l, verbose);

        chi_values = new acb_t[level + 1];
        for(int k = 0; k <= level; k++) {
            acb_init(chi_values[k]);
            _chi.value(chi_values[k], k);
        }
        divisor_counts = new int[level + 1];

        DirichletGroup G(conductor, prec);
        DirichletCharacter chip = G.character(primitive_index);
        chip_values = new acb_t[conductor + 1];
        for(int k = 0; k <= conductor; k++) {
            acb_init(chip_values[k]);
            chip.value(chip_values[k], k);
        }

        divisors_of_level = divisors(level);
        for(int M : divisors_of_level) {
            divisor_counts[M] = ndivisors(M);
            if(M % conductor == 0 && M != level) {
                sublevels.push_back(M);
                DirichletGroup G(M, prec);
                DirichletCharacter psi = G.character(G.index_from_primitive_character(conductor, primitive_index));
                cuspforms_acb * S2 = get_cuspforms_acb(psi, weight, verbose);
                subspaces.push_back(S2);
            }
        }
    }

    ~cuspforms_acb() {
        delete [] psi_table;
        for(int k = 1; k < level + 1; k++) {
            if(level % k == 0) {
                delete [] gcd_tables[k];
            }
        }
        delete [] gcd_tables;
        delete [] phi_table;
        for(int k = 0; k <= level; k++) {
            acb_clear(chi_values[k]);
        }
        delete [] chi_values;
        delete [] divisor_counts;
        for(int k = 0; k <= conductor; k++) {
            acb_clear(chip_values[k]);
        }
        delete [] chip_values;

        if(traces_size > 0) {
            for(int k = 0; k < traces_size; k++) {
                acb_clear(traces[k]);
            }
            delete [] traces;
        }
    }
};

#endif
