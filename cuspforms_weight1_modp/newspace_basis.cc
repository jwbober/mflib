#include <vector>

#include "cuspforms_weight1_modp.h"

#include "flint-extras.h"

using std::vector;

#include <iostream>
using std::cout;
using std::endl;


int full_rank_square_basis(nmod_mat_t square_basis_mat, nmod_mat_t basis_mat) {
    // We're assuming that the first column of basis_mat is all zeros, since
    // these are cusp forms. This is silly but we do that everywhere.
    int dim = nmod_mat_nrows(square_basis_mat);
    int ncoeffs = nmod_mat_ncols(basis_mat);
    if(dim == 0) return 1;
    if(ncoeffs < dim + 1) return 0;

    for(int k = 0; k < dim; k++) {
        for(int j = 1; j < dim + 1; j++) {
            nmod_mat_entry(square_basis_mat, k, j - 1) = nmod_mat_entry(basis_mat, k, j);
        }
    }
    if(nmod_mat_det(square_basis_mat) != 0) {
        return 1;
    }
    int j = 1;
    int r = 0;
    nmod_mat_zero(square_basis_mat);
    while(r < dim) {
        if(j > ncoeffs) return 0;
        for(int k = 0; k < dim; k++) {
            nmod_mat_entry(square_basis_mat, k, r) = nmod_mat_entry(basis_mat, k, j);
        }
        r = nmod_mat_rank(square_basis_mat);
        j++;
    }
    return 1;
}

void cuspforms_weight1_modp::hecke_matrix2(nmod_mat_t Tp, int n, nmod_mat_t basis_mat) {
    // compute a hecke matrix with respect to the given basis.
    // For testing purposes.
    int dim = nmod_mat_nrows(basis_mat);
    nmod_mat_init(Tp, dim, dim, p);
    if(dim == 0) {
        return;
    }

    int ncoeffs = nmod_mat_ncols(basis_mat);

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

    nmod_mat_clear(Tp_basis);
}


int cuspforms_weight1_modp::newspace_basis(nmod_mat_t basis_mat, int ncoeffs) {
    // return the number of newforms (normalized eigenforms) in the basis.
    // these newforms will be the last rows of the returned basis.
    int d = dimension();
    int newd = new_dimension();
    if(newd == 0) {
        nmod_mat_init(basis_mat, 0, 0, p);
        return 0;
    }
    //if(subspaces.size() == 0) {
    //    return 0;
        //std::cout << "newspace basis not implemented yet." << std::endl;
        //exit(-1);
    //}

    // Once we find a Hecke operator that acts invertibly
    // on the new subspace, we should be able to put it into
    // block diagonal form, isolating the new subspace.

    nmod_mat_t hecke_mat;
    nmod_mat_init(hecke_mat, d, d, p);

    nmod_mat_t Tn;
    int n = 2;
    hecke_matrix(Tn, n);

    nmod_mat_add(hecke_mat, hecke_mat, Tn);
    nmod_mat_clear(Tn);

    nmod_mat_t new_window;
    nmod_mat_t old_window;
    nmod_poly_t new_poly; nmod_poly_init(new_poly, p);
    nmod_poly_t old_poly; nmod_poly_init(old_poly, p);

    nmod_poly_t gcd; nmod_poly_init(gcd, p);

    nmod_mat_window_init(old_window, hecke_mat, 0, 0, d - newd, d - newd);
    nmod_mat_window_init(new_window, hecke_mat, d - newd, d - newd, d, d);
    nmod_mat_charpoly(new_poly, new_window);
    nmod_mat_charpoly(old_poly, old_window);
    nmod_poly_gcd(gcd, new_poly, old_poly);
    while(nmod_poly_evaluate_nmod(new_poly, 0) == 0 ||
          nmod_poly_evaluate_nmod(old_poly, 0) == 0 ||
          !nmod_poly_is_squarefree(new_poly) ||
          !nmod_poly_is_one(gcd)) {
        n++;
        hecke_matrix(Tn, n);
        nmod_mat_add(hecke_mat, hecke_mat, Tn);
        nmod_mat_charpoly(new_poly, new_window);
        nmod_mat_charpoly(old_poly, old_window);
        nmod_poly_gcd(gcd, new_poly, old_poly);
        nmod_mat_clear(Tn);
    }

    nmod_mat_clear(new_window);
    nmod_mat_clear(old_window);

    nmod_mat_t P, Pold, Pinv, X;
    nmod_mat_t Pnew;
    nmod_mat_init(P, d, d, p);
    nmod_mat_init(Pinv, d, d, p);
    nmod_mat_init(X, d, d, p);
    nmod_mat_window_init(Pold, P, 0, 0, d, d - newd);

    nmod_mat_t Z;
    nmod_mat_init(Z, d, d, p);

    nmod_poly_factor_t factorization;
    nmod_poly_factor_init(factorization);
    nmod_poly_factor(factorization, new_poly);

    // We want to include the eigenforms that we find at the end of
    // the basis, so we insert the linear factors at r2, which we keep
    // at the end, and the nonlinear ones at r1.
    int r1 = d - newd;
    int r2 = d;
    int eigenform_count = 0;
    for(int k = 0; k < factorization->num; k++) {
        int degree = nmod_poly_degree(factorization->p + k);

        if(degree == 1) {
            eigenform_count++;
            r2--;
            nmod_mat_window_init(Pnew, P, 0, r2, d, r2 + 1);
        }
        else {
            nmod_mat_window_init(Pnew, P, 0, r1, d, r1 + degree);
            r1 = r1 + degree;
        }
        nmod_poly_evaluate_mat(Z, factorization->p + k, hecke_mat);
        nmod_mat_transpose(Z, Z);
        int j = nmod_mat_nullspace(Pnew, Z);
        if(j != degree) {
            cout << "ohno nullspace size was not as expected." << endl;
            exit(-1);
        }

        nmod_mat_window_clear(Pnew);
    }
    nmod_poly_factor_clear(factorization);


    nmod_poly_evaluate_mat(Z, old_poly, hecke_mat);
    nmod_mat_transpose(Z, Z);
    nmod_mat_nullspace(Pold, Z);

    nmod_mat_clear(Pold);
    nmod_poly_clear(old_poly);
    nmod_poly_clear(new_poly);

    nmod_mat_transpose(P, P);
    nmod_mat_inv(Pinv, P);

    nmod_mat_mul(X, P, hecke_mat);
    nmod_mat_mul(Z, X, Pinv);

    nmod_mat_t full_basis_mat;

    basis(full_basis_mat, ncoeffs);

    nmod_mat_t basis_transformation;
    nmod_mat_window_init(basis_transformation, P, d - newd, 0, d, d);

    //nmod_mat_t full_basis_mat2;
    //nmod_mat_init_set(full_basis_mat2, full_basis_mat);
    //nmod_mat_mul(full_basis_mat2, P, full_basis_mat);

    //for(int n = 2; n < 20; n++) {
    //    nmod_mat_clear(hecke_mat);
    //    hecke_matrix2(hecke_mat, n, full_basis_mat2);
    //    //nmod_mat_mul(hecke_mat, hecke_mat, basis_transformation);
    //    cout << "T" << n << ": ";
    //    nmod_mat_print_pretty(hecke_mat);
    //}

    nmod_mat_init(basis_mat, newd, ncoeffs + 1, p);
    nmod_mat_mul(basis_mat, basis_transformation, full_basis_mat);
    for(int k = 0; k < newd; k++) {
        if(nmod_mat_entry(basis_mat, k, 1) != 0) {
            long x = nmod_inv(nmod_mat_entry(basis_mat, k, 1), modp);
            for(int j = 0; j < ncoeffs; j++) {
                nmod_mat_entry(basis_mat, k, j) = nmod_mul(x, nmod_mat_entry(basis_mat, k, j), modp);
            }
        }
    }

    nmod_mat_clear(basis_transformation);
    nmod_mat_clear(hecke_mat);
    nmod_mat_clear(full_basis_mat);
    nmod_mat_clear(P);
    nmod_mat_clear(Pinv);
    nmod_mat_clear(X);
    nmod_mat_clear(Z);

    return eigenform_count;
}
