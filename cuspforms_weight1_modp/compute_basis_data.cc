#include "cuspforms_weight1_modp.h"

#include <iostream>

using namespace std;

void cuspforms_weight1_modp::compute_basis_data() {
    int d3 = S3->dimension();
    int d4 = S4->dimension();
    int Q = d3 + d4 + 10;

    nmod_mat_t basis3;
    nmod_mat_t basis4;
    S3->basis(basis3, Q);
    S4->basis(basis4, Q);


    if(d3 != nmod_mat_rank(basis3) || d4 != nmod_mat_rank(basis4)) {
        cout << "ohno didn't get a matrix of full rank." << endl;
        cout << "d3 = " << d3 << endl;
        cout << "for chi3 = (" << q3 << ", " << chi3.m << ") rank(basis3) = " << nmod_mat_rank(basis3) << endl;
        cout << "d4 = " << d4 << endl;
        cout << "for chi4 = (" << q3 << ", " << chi4.m << ") rank(basis4) = " << nmod_mat_rank(basis4) << endl;
        cout << "I think maybe this shouldn't happen? If it does, see extra_coefficients larger." << endl;
        nmod_mat_clear(basis3);
        nmod_mat_clear(basis4);
        exit(-1);
    }

    nmod_poly_t E3;  nmod_poly_init(E3, p);
    nmod_poly_t E4;  nmod_poly_init(E4, p);
    nmod_poly_t f;   nmod_poly_init(f, p);
    nmod_poly_t g;   nmod_poly_init(g, p);

    eisenstein3_weight1_modp(E3, Q);
    eisenstein4_weight1_modp(E4, Q);

    nmod_mat_t M;
    nmod_mat_init(M, Q, d3 + d4, p);

       for(int j = 0; j < d3; j++) {
        nmod_poly_zero(g);
        for(int k = 0; k < Q; k++) {
            nmod_poly_set_coeff_ui(g, k, nmod_mat_entry(basis3, j, k));
        }
        nmod_poly_mul(f, g, E4);
        for(int k = 0; k < Q; k++) {
            //nmod_mat_entry(M, j, k) = nmod_poly_get_coeff_ui(f, k);
            nmod_mat_entry(M, k, j) = nmod_poly_get_coeff_ui(f, k);
        }
    }

    for(int j = d3; j < d3 + d4; j++) {
        nmod_poly_zero(g);
        for(int k = 0; k < Q; k++) {
            nmod_poly_set_coeff_ui(g, k, nmod_mat_entry(basis4, j - d3, k));
        }
        nmod_poly_mul(f, g, E3);
        for(int k = 0; k < Q; k++) {
            //nmod_mat_entry(M, j, k) = nmod_poly_get_coeff_ui(f, k);
            nmod_mat_entry(M, k, j) = nmod_poly_get_coeff_ui(f, k);
        }
    }

    _dimension = d3 + d4 - nmod_mat_rank(M);
    nmod_mat_t nullspace;
    nmod_mat_init(nullspace, d3 + d4, _dimension, p);
    nmod_mat_nullspace(nullspace, M);

    nmod_mat_t half_nullspace;

    nmod_mat_window_init(half_nullspace, nullspace, 0, 0, d3, _dimension);
    nmod_mat_init(basis_transformation, _dimension, d3, p);
    nmod_mat_transpose(basis_transformation, half_nullspace);

    nmod_mat_window_clear(half_nullspace);

    nmod_mat_clear(M);
    nmod_mat_clear(basis3);
    nmod_mat_clear(basis4);
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(E3);
    nmod_poly_clear(E4);
}
