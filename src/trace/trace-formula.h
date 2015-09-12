#ifndef __TRACE_FORMULA_H__
#define __TRACE_FORMULA_H__

#include <vector>
#include <complex>

#include <characters.h>
#include <flint/nmod_mat.h>


void trace_Tn_unsieved_weight2(std::complex<double> * traces, int start, int end, int level, DirichletCharacter& chi);
void trace_Tn_modp_unsieved_weight2(int * traces, int start, int end, int level, int p, int * chi_values, DirichletCharacter& chi);
void sieve_trace_Tn_modp_on_weight2_for_newspaces(std::vector<int> * traces, int start, int end, int level, int p, int ** chi_values, DirichletCharacter& chi);
int newspace_bases_weight2_modp(nmod_mat_t * bases, int& ncoeffs, int level, int& p, DirichletCharacter& chi, int extra_rows = 0);
void cuspform_basis_weight2_modp(nmod_mat_t basis, int ncoeffs, int level, int& p, DirichletCharacter& chi);
int trace_TmTn_mod_p(int * traces,
                     int * chi_values,
                     int k,
                     int level,
                     int m,
                     int n,
                     int p);

#endif
