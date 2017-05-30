#include <iostream>
#include <vector>

#include "cuspforms_modp.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

void cuspforms_modp::hecke_matrix(nmod_mat_t Tp, int l) {
    int dim = new_dimension();
    vector<int> basis_data = newspace_basis_data(false);
    if(verbose) {
        cout << "using basis rows";
        for(int k : basis_data) {
            cout << " " << k;
        }
        cout << endl;
    }
    compute_traces(l * basis_data[dim - 1] * basis_data[dim - 1] + 5);

    nmod_mat_t basis;
    nmod_mat_init(basis, dim, dim, p);
    for(int k = 0; k < dim; k++) {
        for(int j = k; j < dim; j++) {
            int n = basis_data[k];
            int m = basis_data[j];
            long z = trace_TnTm(n, m);
            nmod_mat_entry(basis, j, k) = z;
            nmod_mat_entry(basis, k, j) = z;
        }
    }

    nmod_mat_t Tp_basis;
    nmod_mat_init(Tp_basis, dim, dim, p);

    for(int k = 0; k < dim; k++) {
        for(int j = k; j < dim; j++) {
            int n = basis_data[k];
            int m = basis_data[j];
            long z1 = trace_TpTnTm(l, n, m);
            nmod_mat_entry(Tp_basis, j, k) = z1;
            nmod_mat_entry(Tp_basis, k, j) = z1;
        }
    }

    nmod_mat_init(Tp, dim, dim, p);
    nmod_mat_solve(Tp, basis, Tp_basis);

    nmod_mat_clear(basis);
    nmod_mat_clear(Tp_basis);
}
