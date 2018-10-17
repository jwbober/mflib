#include "cuspforms_weight1_modp.h"

void cuspforms_weight1_modp::basis(nmod_mat_t basis_mat, int ncoeffs) {
    if(_dimension == -1) {
        compute_basis_data();
    }
    if(_dimension == 0) {
        nmod_mat_init(basis_mat, 0, 0, p);
        return;
    }
    nmod_mat_init(basis_mat, _dimension, ncoeffs + 1, p);

    int r = 0;

    for(int k = 0; k < subspaces.size(); k++) {
        cuspforms_weight1_modp * subspace;
        if(k < subspaces.size()) {
            subspace = subspaces[k];
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
    nmod_mat_t simplebasis_mat;
    simplebasis(simplebasis_mat, ncoeffs);

    int nrows_needed = _dimension - r;
    for(int k = 0; k < nrows_needed; k++) {
        for(int j = 0; j < ncoeffs; j++)
            nmod_mat_entry(basis_mat, _dimension - (k+1), j) = nmod_mat_entry(simplebasis_mat, k, j);
    }

    if(nmod_mat_rank(basis_mat) == _dimension) {
        // The basis entries are pretty random, so it seems like
        // we should typically be in this case. If something goes wrong,
        // we'll add one row at a time and check for increasing rank.
        return;
    }

    for(int k = 0; k < nrows_needed; k++) {
        for(int j = 0; j < ncoeffs; j++)
            nmod_mat_entry(basis_mat, _dimension - (k+1), j) = 0;
    }
    int k = 0;
    for(int i = 0; i < nrows_needed; i++) {
        do {
            for(int j = 0; j < ncoeffs; j++)
                nmod_mat_entry(basis_mat, r, j) = nmod_mat_entry(simplebasis_mat, k, j);
            k++;
        } while(nmod_mat_rank(basis_mat) == r);
        r++;
    }
}
