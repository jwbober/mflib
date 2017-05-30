//
// I'm not really sure what I was intending for this function to do. Right now
// it doesn't do anything.
//


#include <iostream>
#include <vector>

#include "cuspforms_modp.h"

using std::cerr;
using std::cout;
using std::endl;
using std::vector;

void cuspforms_modp::newforms(nmod_mat_t forms, int ncoeffs) {
    int dim = new_dimension();
    vector<int> basis_data = newspace_basis_data();
    int coefficients_needed_for_full_rank = basis_data[dim - 1] + 1;
    ncoeffs = std::max(ncoeffs, coefficients_needed_for_full_rank);
    nmod_mat_t basis;
    newspace_basis(basis, ncoeffs);
    nmod_mat_t smallbasis; // flint wants all the matrices to be square
                           // for linear equation solving, so we're going
                           // to pick out a piece of the basis that we know
                           // gives full rank. At the same time, we'll transpose
                           // it.
    nmod_mat_t Tn;
    nmod_mat_t transformed_basis;
    nmod_mat_init(smallbasis, dim, dim, p);
    nmod_mat_init(transformed_basis, dim, dim, p);
    nmod_mat_init(Tn, dim, dim, p);


    for(int k = 0; k < dim; k++) {
        for(int j = 0; j < dim; j++) {
            nmod_mat_entry(smallbasis, j, k) = nmod_mat_entry(basis, k, basis_data[j]);
        }
    }

    int n = 2;
    nmod_mat_t X, Y;
    nmod_mat_init(X, dim, dim, p);
    nmod_mat_init(Y, dim, dim, p);
    while(n < 20) {
        nmod_mat_zero(transformed_basis);
        int max_coeff_needed = coefficients_needed_for_full_rank * n;
        if(max_coeff_needed > ncoeffs) {
            cerr << "need to compute more coefficients." << endl;
            exit(0);
        }
        for(int j = 0; j < dim; j++) {
            for(int k = 0; k < dim; k++) {
                int m = basis_data[k];
                for(auto d : divisors(GCD(m, n))) {
                    d = d % p;
                    long z = chi_values[d % level];
                    z = nmod_mul(z, nmod_pow_ui(d, weight - 1, modp), modp);
                    z = nmod_mul(z, nmod_mat_entry(basis, j, m*n/(d*d)), modp);
                    nmod_mat_entry(transformed_basis, k, j) =
                        nmod_add(nmod_mat_entry(transformed_basis, k, j), z, modp);
                    //transformed_basis(j, m) += chi_values[d % level] * pow((double)d, weight - 1) * basis(j, m*n/(d*d));
                }
            }
        }
        nmod_mat_solve(Tn, smallbasis, transformed_basis);
        //nmod_mat_print_pretty(Tn);
        for(int a = 0; a <= 2*n; a++) {
            for(int k = 0; k < dim; k++) {
                nmod_mat_entry(X, k, k) = a;
            }
            nmod_mat_sub(Y, Tn, X);
            if(nmod_mat_det(Y) == 0) cout << a << " ";
            if(a == 0) continue;
            a = p - a;
            for(int k = 0; k < dim; k++) {
                nmod_mat_entry(X, k, k) = a;
            }
            nmod_mat_sub(Y, Tn, X);
            if(nmod_mat_det(Y) == 0) cout << a - p << " ";
            a = p - a;
        }
        n++;
    }

    nmod_mat_clear(Tn);
    nmod_mat_clear(basis);
    nmod_mat_clear(smallbasis);
    nmod_mat_clear(transformed_basis);
}
