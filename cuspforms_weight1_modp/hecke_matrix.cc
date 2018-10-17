#include "cuspforms_weight1_modp.h"

#include <iostream>
using namespace std;

void cuspforms_weight1_modp::hecke_matrix(nmod_mat_t Tp, int n) {
    int dim = dimension();
    nmod_mat_init(Tp, dim, dim, p);
    if(dim == 0) {
        return;
    }

    int ncoeffs = (5 + 2*dim)*n;  // I don't know what we want here,
                                  // but this is a starting point.
    nmod_mat_t basis_mat;
    basis(basis_mat, ncoeffs);

    nmod_mat_t square_basis_mat;  // we are going to need a full rank
                                  // square matrix for the flint functions
                                  // below, so...
                                  // (maybe this is not true. Possibly the
                                  //  functions I was calling weren't working
                                  //  because I forgot to transpose the matrices.)

    nmod_mat_init(square_basis_mat, dim, dim, p);
    vector<int> basis_coeff_indices;
    for(int j = 1; j < dim + 1; j++) {
        basis_coeff_indices.push_back(j);
    }
    for(int k = 0; k < dim; k++) {
        for(int j = 1; j < dim + 1; j++) {
            nmod_mat_entry(square_basis_mat, k, j - 1) = nmod_mat_entry(basis_mat, k, j);
        }
    }
    if(nmod_mat_rank(square_basis_mat) != dim) {
        basis_coeff_indices.clear();
        int j = 1;
        int r = 0;
        nmod_mat_zero(square_basis_mat);
        while(r < dim) {
            for(int k = 0; k < dim; k++) {
                nmod_mat_entry(square_basis_mat, k, r) = nmod_mat_entry(basis_mat, k, j);
            }
            int r2 = nmod_mat_rank(square_basis_mat);
            if(r2 > r) basis_coeff_indices.push_back(j);
            r = r2;
            j++;
        }
    }

    if(basis_coeff_indices[dim - 1] * n > ncoeffs) {
        cout << "not implemented yet." << endl;
        exit(-1);
    }

    nmod_mat_t Tp_basis;
    nmod_mat_init(Tp_basis, dim, dim, p);

    for(int k = 0; k < dim; k++) {
        for(int j = 0; j < dim; j++) {
            int m = basis_coeff_indices[j];
            long S = 0;
            for(int d : divisors(int(GCD(n, m)))) {
                if(GCD(d, level) > 1) continue;
                long z = chi_values[d];
                S = S + nmod_mul(z, nmod_mat_entry(basis_mat, k, m*n/(d*d)), modp);
                S = S % p;
            }
            nmod_mat_entry(Tp_basis, k, j) = S;
        }
    }

    nmod_mat_transpose(square_basis_mat, square_basis_mat);
    nmod_mat_transpose(Tp_basis, Tp_basis);
    nmod_mat_init(Tp, dim, dim, p);
    int x = nmod_mat_solve(Tp, square_basis_mat, Tp_basis);
    nmod_mat_transpose(Tp, Tp);

    nmod_mat_clear(basis_mat);
    nmod_mat_clear(Tp_basis);
}
