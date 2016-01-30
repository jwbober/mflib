#ifndef __TRACE_FORMULA_H__
#define __TRACE_FORMULA_H__

#include <vector>
#include <complex>
#include <Eigen/Dense>

#include <characters.h>
#include <flint/nmod_mat.h>
#include "acb_mat.h"


void trace_Tn_unsieved_weight2(
        std::complex<double> * traces,
        int start,
        int end,
        int level,
        std::complex<double> * chi_values,
        DirichletCharacter& chi,
        int verbose = 0);

void trace_Tn_modp_unsieved_weight2(
        long * traces,
        int start,
        int end,
        int level,
        long p,
        long * chi_values,
        DirichletCharacter& chi,
        int verbose = 0);

void sieve_trace_Tn_modp_on_weight2_for_newspaces(
        std::vector<long> * traces,
        int start,
        int end,
        int level,
        long p,
        long ** chi_values,
        DirichletCharacter& chi,
        int verbose = 0);

void sieve_trace_Tn_on_weight2_for_newspaces(
        std::vector< std::complex<double> > * traces,
        int start,
        int end,
        int level,
        std::complex<double> ** chi_values,
        DirichletCharacter& chi,
        int verbose = 0);



int newspace_bases_weight2_modp(
        nmod_mat_t * bases,
        int& ncoeffs,
        int level,
        long& p,
        DirichletCharacter& chi,
        int extra_rows = 0,
        int verbose = 0);

void cuspform_basis_weight2_modp(
        nmod_mat_t basis,
        int ncoeffs,
        int level,
        long& p,
        DirichletCharacter& chi,
        int verbose = 0);

long trace_TmTn_mod_p(
        long * traces,
        long * chi_values,
        int k,
        int level,
        int m,
        int n,
        long p);

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> newspace_basis_weight2(
        std::vector<int>& rows,
        int& ncoeffs,
        int level,
        DirichletCharacter& chi,
        int verbose = 0);

#endif
