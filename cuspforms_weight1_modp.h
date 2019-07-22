#ifndef __CUSPFORMS_WEIGHT1_MODP_H__
#define __CUSPFORMS_WEIGHT1_MODP_H__

#include <vector>

#include "flint/flint.h"
#include "flint/nmod_mat.h"
#include "flint/nmod_poly.h"

#include "cuspforms_modp.h"

#include "slint.h"
#include "characters.h"



class cuspforms_weight1_modp;
cuspforms_weight1_modp * get_cuspforms_weight1_modp(DirichletCharacter &chi, long p, int verbose = 0);
void clear_cuspforms_weight1_modp();

void eisenstein3_weight1_modp(nmod_poly_t E3, unsigned int degree);
void eisenstein4_weight1_modp(nmod_poly_t E4, unsigned int degree);

class cuspforms_weight1_modp {
public:
    nmod_t modp;
    long p;
    int level;
    const int weight = 1;
    int chi;
    int conductor;
    int verbose;
    int _dimension = -1;
    int zero_for_trivial_reasons = 0;

    std::vector<int> sublevels;
    std::vector<int> divisors_of_level;
    std::vector<cuspforms_weight1_modp*> subspaces;
    std::vector<DirichletGroup> sublevel_dirichlet_groups;

    int * divisor_counts;
    int * phi_table;
    int ** gcd_tables;
    long * chi_values;
    long * chip_values;

    long q3;
    long q4;

    DirichletCharacter chi3;
    DirichletCharacter chi4;

    cuspforms_modp * S3;
    cuspforms_modp * S4;

    nmod_mat_t basis_transformation;
    int basis_data_computed = false;

    void compute_basis_data(); // Fill the matrix basis_transformation

    //long trace(int n);              // Get Tn, computing it first if necessary
    //long trace_TnTm(int n, int m);  // Trace of the product of two Hecke operators
    //long trace_TpTnTm(int p, int n, int m); // Trace of the product of three Hecke operators

    //void compute_traces(int end);

    int dimension();
    int new_dimension();

    int newspace_basis(nmod_mat_t basis, int ncoeffs);
    void newforms(nmod_mat_t forms, int ncoeffs);

    // hecke matrix with respect to the default basis of the full space.
    // Note that our hecke matrix acts on the LEFT.
    void hecke_matrix(nmod_mat_t Tp, int n);
    void hecke_matrix2(nmod_mat_t Tp, int n, nmod_mat_t basis_mat);

    const std::vector<int>& newspace_basis_data(bool coprime_only = true);

    // There are two options for computing a basis of the full space of
    // cusp forms. The "simple" option computes a basis of the space as cusp forms
    // in higher level and weight divided by the Eisenstein series E3.
    //
    // The default option (which will use the simple basis internally) is to
    // compute a basis which consists consists as much as possible of forms
    // lifted from lower levels.

    void basis(nmod_mat_t basis, int ncoeffs);
    void simplebasis(nmod_mat_t basis, int ncoeffs);

    cuspforms_weight1_modp(DirichletCharacter &_chi, long _p, int _verbose = 0);
    ~cuspforms_weight1_modp();
};

#endif
