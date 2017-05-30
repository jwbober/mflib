#ifndef __MODFORM_MODP_H__
#define __MODFORM_MODP_H__

#include <vector>

#include "flint/flint.h"
#include "flint/nmod_mat.h"

#include "slint.h"
#include "characters.h"

class cuspforms_modp;
cuspforms_modp * get_cuspforms_modp(DirichletCharacter &chi, int weight, long p, int verbose = 0);
void clear_cuspforms_modp();

class cuspforms_modp {
public:
    nmod_t modp;
    long p;
    int level;
    int weight;
    int chi;
    int conductor;
    int verbose;
    int _dimension = -1;

    std::vector<int> sublevels;
    std::vector<int> divisors_of_level;
    std::vector<cuspforms_modp*> subspaces;
    std::vector<DirichletGroup> sublevel_dirichlet_groups;

    std::vector<int> basis_rows;    // list of integers to specify a basis of traces prime to the level
    std::vector<int> basis_rows2;   // list of integers to speficy a basis with no coprimality condition

    std::vector<long> traces;       // (sieved) trace of Hecke operators

    int * divisor_counts;
    long * psi_table;
    int * phi_table;
    int ** gcd_tables;
    long * chi_values;
    long * chip_values;

    long trace(int n);              // Get Tn, computing it first if necessary
    long trace_TnTm(int n, int m);  // Trace of the product of two Hecke operators
    long trace_TpTnTm(int p, int n, int m); // Trace of the product of three Hecke operators

    void compute_traces(int end);

    int dimension();
    int new_dimension();

    void newspace_basis(nmod_mat_t basis, int ncoeffs);
    void newforms(nmod_mat_t forms, int ncoeffs); // I'm not sure what I was
                                                  // planning for this function to do.

    void hecke_matrix(nmod_mat_t Tp, int n);

    const std::vector<int>& newspace_basis_data(bool coprime_only = true);
    void basis(nmod_mat_t basis, int ncoeffs);

    long evalpoly(long t, long n);
    cuspforms_modp(DirichletCharacter &_chi, const int w, long _p, int _verbose = 0);
    ~cuspforms_modp();
};

#endif
