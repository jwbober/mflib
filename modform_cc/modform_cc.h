#ifndef __MODFORM_CC_H__
#define __MODFORM_CC_H__

#include <vector>
#include <complex>

#include "slint.h"
#include "characters.h"
#include "modform_modp.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> cmatrix_t;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> rmatrix_t;

class cuspforms_cc;
cuspforms_cc * get_cuspforms_cc(DirichletCharacter &chi, int weight, int verbose = 0);

class cuspforms_cc {
public:
    int level;
    int weight;
    int chi;
    int conductor;
    int verbose;

    std::vector<int> sublevels;
    std::vector<int> divisors_of_level;
    std::vector<cuspforms_cc*> subspaces;
    std::vector<DirichletGroup> sublevel_dirichlet_groups;

    int * divisor_counts;
    int * psi_table;
    int * phi_table;
    int ** gcd_tables;
    std::complex<double> * chi_values;
    std::complex<double> * chip_values;

    std::vector<int> basis_rows;

    std::vector<std::complex<double> > traces;

    std::complex<double> trace(int n);
    std::complex<double> trace_TnTm(int n, int m);
    void compute_traces(int end);

    int dimension();
    int new_dimension();

    cmatrix_t newspace_basis(int ncoeffs);
    const std::vector<int>& newspace_basis_data();
    cmatrix_t basis(int ncoeffs);

    cuspforms_modp * modp_space;
    cmatrix_t newforms(int ncoeffs);

    double evalpoly(long t, long n);

    cuspforms_cc(DirichletCharacter &_chi, const int w, int _verbose = 0) {
        verbose = _verbose;
        chi = _chi.m;
        level = _chi.parent->q;
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

        modp_space = get_cuspforms_modp(_chi, weight, 10000000l);

        chi_values = new std::complex<double>[level + 1];
        for(int k = 0; k <= level; k++) {
            chi_values[k] = _chi.value(k);
        }
        divisor_counts = new int[level + 1];

        DirichletGroup G(conductor);
        DirichletCharacter chip = G.character(primitive_index);
        chip_values = new std::complex<double>[conductor + 1];
        for(int k = 0; k <= conductor; k++) {
            chip_values[k] = chip.value(k);
        }

        divisors_of_level = divisors(level);
        for(int M : divisors_of_level) {
            divisor_counts[M] = ndivisors(M);
            if(M % conductor == 0 && M != level) {
                sublevels.push_back(M);
                DirichletGroup G(M);
                DirichletCharacter psi = G.character(G.index_from_primitive_character(conductor, primitive_index));
                cuspforms_cc * S2 = get_cuspforms_cc(psi, weight, verbose);
                subspaces.push_back(S2);
            }
        }
    }

    ~cuspforms_cc() {
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
