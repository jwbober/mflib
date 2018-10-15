#include <vector>
#include <iostream>

#include "cuspforms_modp.h"

using std::vector;

void cuspforms_modp::newspace_basis(nmod_mat_t basis, int ncoeffs) {
    newspace_basis_data();
    int d = new_dimension();
    if(d == 0) {
        nmod_mat_init(basis, 0, 0, p);
        return;
    }

    compute_traces(ncoeffs * basis_rows[d-1]);
    nmod_mat_init(basis, d, ncoeffs, p);
    for(int k = 0; k < d; k++) {
        int n = basis_rows[k];
        for(int m = 0; m < ncoeffs; m++) {
            nmod_mat_entry(basis, k, m) = trace_TnTm(n, m);
        }
    }
}

void cuspforms_modp::basis(nmod_mat_t basis_mat, int ncoeffs) {
    nmod_mat_init(basis_mat, dimension(), ncoeffs + 1, p);

    int r = 0;

    for(int k = 0; k <= subspaces.size(); k++) {
        cuspforms_modp * subspace;
        if(k < subspaces.size()) {
            subspace = subspaces[k];
        }
        else {
            subspace = this;
        }
        nmod_mat_t subspace_basis;
        subspace->newspace_basis(subspace_basis, ncoeffs);
        int M1 = subspace->level;
        for(auto M2 : divisors(level/M1)) {
            for(int l = 0; l < nmod_mat_nrows(subspace_basis); l++) {
                for(int j = 1; j*M2 < ncoeffs; j++) {
                    nmod_mat_entry(basis_mat, r, j * M2) = nmod_mat_entry(subspace_basis, l, j);
                }
                r++;
            }
        }
        nmod_mat_clear(subspace_basis);
    }

}

const vector<int>& cuspforms_modp::newspace_basis_data(bool coprime_only) {
    if(basis_rows.size() > 0 && coprime_only) return basis_rows;
    if(basis_rows2.size() > 0 && !coprime_only) return basis_rows2;
    int d = new_dimension();
    if(d == 0) return basis_rows;
    nmod_mat_t basis;
    nmod_mat_init(basis, d, d, p);
    compute_traces(d*d);
    int n = 1;
    int row = 0;
    vector<int> *basis_rows_ptr;
    if(coprime_only) {
        basis_rows_ptr = &basis_rows;
    }
    else {
        basis_rows_ptr = &basis_rows2;
    }
    while(row < d) {
        if(coprime_only) {while(GCD(n, level) > 1) n++;}
        basis_rows_ptr->push_back(n);
        int m = 1;
        int col = 0;
        while(col < d) {
            if(coprime_only) {while(GCD(m, level) > 1) n++;}
            //for(int m = 1; m < d + 1; m++) {
            nmod_mat_entry(basis, row, col) = trace_TnTm(n, m);
            col++;
        }
        row++;
        n++;
    }

    if(nmod_mat_det(basis) != 0) {
        nmod_mat_clear(basis);
        return *basis_rows_ptr;
    }

    int k = 0;
    n = 1;
    basis_rows_ptr->clear();

    while(k < d) {
        if(coprime_only) {while(GCD(n, level) != 1) n++;}
        for(int j = 0; j < k; j++) {
            //while(GCD(m, level) != 1) m++;
            int m = basis_rows_ptr->at(j);
            long t = trace_TnTm(n,m);
            nmod_mat_entry(basis, k, j) = t;
            nmod_mat_entry(basis, j, k) = t;
            m++;
        }
        nmod_mat_entry(basis, k, k) = trace_TnTm(n, n);
        nmod_mat_t topleft;
        nmod_mat_window_init(topleft, basis, 0, 0, k + 1, k + 1);
        if(nmod_mat_det(topleft) != 0) {
            k++;
            basis_rows_ptr->push_back(n);
        }
        n++;
        nmod_mat_window_clear(topleft);
    }
    nmod_mat_clear(basis);
    return *basis_rows_ptr;
}
