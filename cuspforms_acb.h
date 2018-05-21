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
#include "cuspforms_modp.h"
#include "acb_mat.h"

class cuspforms_acb;
cuspforms_acb * get_cuspforms_acb(DirichletCharacter &chi, int weight, int nthreads = 1, int verbose = 0);
void clear_cuspform_cache();

class cuspforms_acb {
public:
    int level;
    int weight;
    int chi;
    int conductor;
    int nthreads;
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
    long _new_dimension = -1;
    long new_dimension();

    void newspace_basis(acb_mat_t B, int ncoeffs);
    const std::vector<int>& newspace_basis_data();
    //cmatrix_t basis(int ncoeffs);

    cuspforms_modp * modp_space;
    void newforms(acb_mat_t out, int ncoeffs);

    void evalpoly(fmpz_t out, fmpz_t t, fmpz_t n);

    cuspforms_acb(DirichletCharacter &_chi, const int w, int _nthreads = 1, int _verbose = 0);
    ~cuspforms_acb();
};

#endif
