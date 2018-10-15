#include "cuspforms_weight1_modp.h"

void cuspforms_weight1_modp::simplebasis(nmod_mat_t basis, int ncoeffs) {
    if(_dimension == -1) {
        compute_basis_data();
    }
    if(_dimension == 0) {
        nmod_mat_init(basis, 0, 0, p);
        return;
    }
    nmod_mat_t basis3;
    S3->basis(basis3, ncoeffs);

    nmod_mat_init(basis, nmod_mat_nrows(basis_transformation), nmod_mat_ncols(basis3), p);
    nmod_mat_mul(basis, basis_transformation, basis3);

    nmod_poly_t E3, f, g;
    nmod_poly_t E3_inverse;

    nmod_poly_init(E3, p);
    nmod_poly_init(f, p);
    nmod_poly_init(g, p);
    nmod_poly_init(E3_inverse, p);

    eisenstein3_weight1_modp(E3, ncoeffs);
    nmod_poly_inv_series(E3_inverse, E3, ncoeffs);

    for(int k = 0; k < _dimension; k++) {
        for(int j = 0; j < ncoeffs; j++) {
            nmod_poly_set_coeff_ui(g, j, nmod_mat_entry(basis, k, j));
        }
        nmod_poly_mul(f, g, E3_inverse);
        unsigned long z = 0;
        int j = 1;
        while(z == 0) {
            z = nmod_poly_get_coeff_ui(f, j);
            j++;
        }
        z = InvMod(z, p);
        nmod_poly_scalar_mul_nmod(f, f, z);
        for(int j = 0; j < ncoeffs; j++) {
            z = nmod_poly_get_coeff_ui(f, j);
            nmod_mat_entry(basis, k, j) = z;
        }
    }

    nmod_poly_clear(E3);
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(E3_inverse);
    nmod_mat_clear(basis3);
}
