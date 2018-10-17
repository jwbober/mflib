#include "cuspforms_weight1_modp.h"

#include <iostream>

void cuspforms_weight1_modp::newspace_basis(nmod_mat_t basis_mat, int ncoeffs) {
    if(new_dimension() == 0) {
        nmod_mat_init(basis_mat, 0, 0, p);
        return;
    }
    if(subspaces.size() != 0) {
        std::cout << "newspace basis not implemented yet." << std::endl;
        exit(-1);
    }

    simplebasis(basis_mat, ncoeffs);
}
