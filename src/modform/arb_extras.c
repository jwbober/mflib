#include "arb_mat.h"

//int arb_cmp_for_sort(arb_struct a, arb_struct b) {
//    return arf_cmp(arb_midref(&a), arb_midref(&b));
//}

//void arb_sort(arb_ptr x, size_t count) {
//    qsort(x, count, sizeof(arb_struct), arb_cmp_for_sort);
//}

void fprintarf(FILE *fp,const arf_t x) {
    static int init;
    static fmpz_t m,e;

    if (!init) {
        fmpz_init(m); fmpz_init(e);
        init = 1;
    }
   arf_get_fmpz_2exp(m,e,x);
   fmpz_fprint(fp,m); fprintf(fp," "); fmpz_fprint(fp,e);
}

void fprintarb(FILE *fp,const arb_t x) {
    static int init;
    static arf_t a;

    if (!init) {
        arf_init(a);
        init = 1;
    }

    fprintarf(fp,arb_midref(x));
    fprintf(fp," ");
    arf_set_mag(a,arb_radref(x));
    fprintarf(fp,a);
}

void fprintacb(FILE *fp,const acb_t x) {
    fprintarb(fp,acb_realref(x));
    fprintf(fp," ");
    fprintarb(fp,acb_imagref(x));
}

void arb_square(arb_t out, const arb_t in, slong prec) {
    mag_mul(arb_radref(out), arb_radref(in), arb_radref(in));
    int inexact = arf_mul(arb_midref(out), arb_midref(in), arb_midref(in), ARB_RND, prec);
    if (inexact)
        arf_mag_add_ulp(arb_radref(out), arb_radref(out), arb_midref(out), prec);

}


void arb_mat_L2norm(arb_t out, const arb_mat_t in, slong prec) {
    int nrows = arb_mat_nrows(in);
    int ncols = arb_mat_ncols(in);
    arb_zero(out);
    for(int i = 0; i < nrows; i++) {
        for(int j = 0; j < ncols; j++) {
            arb_addmul(out, arb_mat_entry(in, i, j), arb_mat_entry(in, i, j), prec);
        }
    }
    arb_sqrtpos(out, out, prec);
}

void arb_set_exact(arb_t x) {
    mag_zero(arb_radref(x));
}
void arb_mat_set_exact(arb_mat_t A) {
    int nrows = arb_mat_nrows(A);
    int ncols = arb_mat_nrows(A);
    for(int i = 0; i < nrows; i++) {
        for(int j = 0; j < ncols; j++) {
            arb_set_exact(arb_mat_entry(A, i, j));
        }
    }
}

void arb_mat_print_sage_float(const arb_mat_t A) {
    int nrows = arb_mat_nrows(A);
    int ncols = arb_mat_ncols(A);
    printf("[");
    for(int j = 0; j < nrows; j++) {
        printf("[");
        for(int k = 0; k < ncols; k++) {
            double x = arf_get_d(arb_midref(arb_mat_entry(A, j, k)), ARF_RND_NEAR);
            printf("%e", x);
            if(k < nrows - 1)
                printf(", ");
        }
        printf("],\n");
    }
    printf("]\n");
}

void arb_mat_cholesky(arb_mat_t out, const arb_mat_t in, slong prec) {
    int nrows = arb_mat_nrows(in);
    for(int j = 0; j < nrows; j++) {
        for(int i = j; i < nrows; i++) {
            arb_set(arb_mat_entry(out, i, j), arb_mat_entry(in, i, j));
            for(int k = 0; k < j; k++) {
                arb_submul(arb_mat_entry(out, i, j), arb_mat_entry(out, i, k), arb_mat_entry(out, j, k), prec);
            }
            if(i == j) {
                arb_sqrt(arb_mat_entry(out, i, j), arb_mat_entry(out, i, j), prec);
            }
            else {
                arb_div(arb_mat_entry(out, i, j), arb_mat_entry(out, i, j), arb_mat_entry(out, j, j), prec);
            }
        }
    }
}

void arb_twobytwo_diag(arb_t u1, arb_t u2, const arb_t a, const arb_t b, const arb_t d, slong prec) {
    // Compute the orthogonal matrix that diagonalizes
    //
    //    A = [a b]
    //        [b d]
    //
    // This matrix will have the form
    //
    //    U = [cos x , -sin x]
    //        [sin x, cos x]
    //
    // where the diagonal matrix is U^t A U.
    // We set u1 = cos x, u2 = -sin x.

    if(arb_contains_zero(b)) {
        // this is not quite right (doesn't set error intervals)
        arb_set_ui(u1, 1);
        arb_set_ui(u2, 0);
        return;
    }
    arb_t x; arb_init(x);

    arb_mul(u1, b, b, prec);            // u1 = b^2
    arb_sub(u2, a, d, prec);            // u2 = a - d
    arb_mul_2exp_si(u2, u2, -1);        // u2 = (a - d)/2
    arb_mul(u2, u2, u2, prec);          // u2 = ( (a - d)/2 )^2
    arb_add(u1, u1, u2, prec);          // u1 = b^2 + ( (a-d)/2 )^2
    arb_sqrt(u1, u1, prec);             // u1 = sqrt(above)

    arb_mul_2exp_si(u1, u1, 1);         // u1 = 2 (sqrt (above) )
    arb_add(u1, u1, d, prec);           // u1 += d
    arb_sub(u1, u1, a, prec);           // u1 -= a
    arb_mul_2exp_si(u1, u1, -1);        // u1 = (d - a)/2 + sqrt(b^2 + ( (a-d)/2 )^2)

    arb_mul(x, u1, u1, prec);
    arb_addmul(x, b, b, prec);          // x = u1^2 + b^2
    arb_sqrt(x, x, prec);               // x = sqrt(u1^2 + b^2)
    arb_div(u2, u1, x, prec);
    arb_div(u1, b, x, prec);
    arb_neg(u1, u1);

    arb_clear(x);
}

void arf_twobytwo_diag(arf_t u1, arf_t u2, const arf_t a, const arf_t b, const arf_t d, slong prec) {
    // Compute the orthogonal matrix that diagonalizes
    //
    //    A = [a b]
    //        [b d]
    //
    // This matrix will have the form
    //
    //    U = [cos x , -sin x]
    //        [sin x, cos x]
    //
    // where the diagonal matrix is U^t A U.
    // We set u1 = cos x, u2 = -sin x.

    if(arf_is_zero(b)) {
        arf_set_ui(u1, 1);
        arf_set_ui(u2, 0);
        return;
    }
    arf_t x; arf_init(x);

    arf_mul(u1, b, b, prec, ARF_RND_NEAR);            // u1 = b^2
    arf_sub(u2, a, d, prec, ARF_RND_NEAR);            // u2 = a - d
    arf_mul_2exp_si(u2, u2, -1);                      // u2 = (a - d)/2
    arf_mul(u2, u2, u2, prec, ARF_RND_NEAR);          // u2 = ( (a - d)/2 )^2
    arf_add(u1, u1, u2, prec, ARF_RND_NEAR);          // u1 = b^2 + ( (a-d)/2 )^2
    arf_sqrt(u1, u1, prec, ARF_RND_NEAR);             // u1 = sqrt(above)

    arf_mul_2exp_si(u1, u1, 1);                       // u1 = 2 (sqrt (above) )
    arf_add(u1, u1, d, prec, ARF_RND_NEAR);           // u1 += d
    arf_sub(u1, u1, a, prec, ARF_RND_NEAR);           // u1 -= a
    arf_mul_2exp_si(u1, u1, -1);                      // u1 = (d - a)/2 + sqrt(b^2 + ( (a-d)/2 )^2)

    arf_mul(x, u1, u1, prec, ARF_RND_NEAR);
    arf_addmul(x, b, b, prec, ARF_RND_NEAR);          // x = u1^2 + b^2
    arf_sqrt(x, x, prec, ARF_RND_NEAR);               // x = sqrt(u1^2 + b^2)
    arf_div(u2, u1, x, prec, ARF_RND_NEAR);
    arf_div(u1, b, x, prec, ARF_RND_NEAR);
    arf_neg(u1, u1);

    arf_clear(x);
}


int arb_mat_jacobi(arb_mat_t D, arb_mat_t P, const arb_mat_t A, slong prec) {
    //
    // Given a d x d real symmetric matrix A, compute an orthogonal matrix
    // P and a diagonal D such that A = P D P^t = P D P^(-1).
    //
    // D should have already been initialized as a d x 1 matrix, and Pp
    // should have already been initialized as a d x d matrix.
    //
    // If the eigenvalues can be certified as unique, then a nonzero int is
    // returned, and the eigenvectors should have reasonable error bounds. If
    // the eigenvalues cannot be certified as unique, then some of the
    // eigenvectors will have infinite error radius.

#define B(i,j) arb_mat_entry(B, i, j)
#define D(i) arb_mat_entry(D, i, 0)
#define P(i,j) arb_mat_entry(P, i, j)
    int dim = arb_mat_nrows(A);
    if(dim == 1) {
        arb_mat_set(D, A);
        arb_mat_one(P);
        return 0;
    }
    arb_mat_t B;
    arb_mat_init(B, dim, dim);

    arf_t * B1 = (arf_t*)malloc(dim * sizeof(arf_t));
    arf_t * B2 = (arf_t*)malloc(dim * sizeof(arf_t));
    arf_t * row_max = (arf_t*)malloc((dim - 1) * sizeof(arf_t));
    int * row_max_indices = (int*)malloc((dim - 1) * sizeof(int));

    for(int k = 0; k < dim; k++) {
        arf_init(B1[k]);
        arf_init(B2[k]);
    }
    for(int k = 0; k < dim - 1; k++) {
        arf_init(row_max[k]);
    }

    arf_t x1, x2;
    arf_init(x1);
    arf_init(x2);

    arf_t Gii, Gij, Gji, Gjj;
    arf_init(Gii);
    arf_init(Gij);
    arf_init(Gji);
    arf_init(Gjj);

    arb_mat_set(B, A);
    arb_mat_one(P);

    for(int i = 0; i < dim - 1; i++) {
        for(int j = i + 1; j < dim; j++) {
            arf_abs(x1, arb_midref(B(i,j)));
            if(arf_cmp(row_max[i], x1) < 0) {
                arf_set(row_max[i], x1);
                row_max_indices[i] = j;
            }
        }
    }


    int finished = 0;

    while(!finished) {
        arf_zero(x1);
        int i = 0;
        int j = 0;
        for(int k = 0; k < dim - 1; k++) {
            if(arf_cmp(x1, row_max[k]) < 0) {
                arf_set(x1, row_max[k]);
                i = k;
            }
        }
        j = row_max_indices[i];

        slong bound = arf_abs_bound_lt_2exp_si(x1);
        if(bound < -prec * .9) {
            finished = 1;
            break;
        }
        else {
            //printf("%ld\n", arf_abs_bound_lt_2exp_si(x1));
            //arb_mat_printd(B, 10);
            //printf("\n");
        }

        arf_twobytwo_diag(Gii, Gij, arb_midref(B(i,i)), arb_midref(B(i,j)), arb_midref(B(j,j)), 2*prec);
        arf_neg(Gji, Gij);
        arf_set(Gjj, Gii);

        //printf("%d %d\n", i, j);
        //arf_printd(Gii, 100);
        //printf(" ");
        //arf_printd(Gij, 100);
        //printf("\n");
        if(arf_is_zero(Gij)) {  // If this happens, we're
            finished = 1;       // not going to do any better
            break;              // without increasing the precision.
        }

        for(int k = 0; k < dim; k++) {
            arf_mul(B1[k], Gii, arb_midref(B(i,k)), prec, ARF_RND_NEAR);
            arf_addmul(B1[k], Gji, arb_midref(B(j,k)), prec, ARF_RND_NEAR);

            arf_mul(B2[k], Gij, arb_midref(B(i,k)), prec, ARF_RND_NEAR);
            arf_addmul(B2[k], Gjj, arb_midref(B(j,k)), prec, ARF_RND_NEAR);
        }
        for(int k = 0; k < dim; k++) {
            arf_set(arb_midref(B(i,k)), B1[k]);
            arf_set(arb_midref(B(j,k)), B2[k]);
        }
        for(int k = 0; k < dim; k++) {
            arf_mul(B1[k], Gii, arb_midref(B(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B1[k], Gji, arb_midref(B(k,j)), prec, ARF_RND_NEAR);

            arf_mul(B2[k], Gij, arb_midref(B(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B2[k], Gjj, arb_midref(B(k,j)), prec, ARF_RND_NEAR);
        }
        for(int k = 0; k < dim; k++) {
            arf_set(arb_midref(B(k,i)), B1[k]);
            arf_set(arb_midref(B(k,j)), B2[k]);
        }

        for(int k = 0; k < dim; k++) {
            arf_mul(B1[k], Gii, arb_midref(P(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B1[k], Gji, arb_midref(P(k,j)), prec, ARF_RND_NEAR);

            arf_mul(B2[k], Gij, arb_midref(P(k,i)), prec, ARF_RND_NEAR);
            arf_addmul(B2[k], Gjj, arb_midref(P(k,j)), prec, ARF_RND_NEAR);
        }
        for(int k = 0; k < dim; k++) {
            arf_set(arb_midref(P(k,i)), B1[k]);
            arf_set(arb_midref(P(k,j)), B2[k]);
        }

        if(i < dim - 1)
            arf_set_ui(row_max[i], 0);
        if(j < dim - 1)
            arf_set_ui(row_max[j], 0);

        // Update the max in any row where the maximum
        // was in a column that changed.
        for(int k = 0; k < dim - 1; k++) {
            if(row_max_indices[k] == j || row_max_indices[k] == i) {
                arf_abs(row_max[k], arb_midref(B(k,k+1)));
                row_max_indices[k] = k+1;
                for(int l = k+2; l < dim; l++) {
                    arf_abs(x1, arb_midref(B(k,l)));
                    if(arf_cmp(row_max[k], x1) < 0) {
                        arf_set(row_max[k], x1);
                        row_max_indices[k] = l;
                    }
                }
            }
        }

        // Update the max in the ith row.
        for(int k = i + 1; k < dim; k++) {
            arf_abs(x1, arb_midref(B(i, k)));
            if(arf_cmp(row_max[i], x1) < 0) {
                arf_set(row_max[i], x1);
                row_max_indices[i] = k;
            }
        }

        // Update the max in the jth row.
        for(int k = j + 1; k < dim; k++) {
            arf_abs(x1, arb_midref(B(j, k)));
            if(arf_cmp(row_max[j], x1) < 0) {
                arf_set(row_max[j], x1);
                row_max_indices[j] = k;
            }
        }

        // Go through column i to see if any of
        // the new entries are larger than the
        // max of their row.
        for(int k = 0; k < i; k++) {
            if(k == dim) continue;
            arf_abs(x1, arb_midref(B(k, i)));
            if(arf_cmp(row_max[k], x1) < 0) {
                arf_set(row_max[k], x1);
                row_max_indices[k] = i;
            }
        }

        // And then column j.
        for(int k = 0; k < j; k++) {
            if(k == dim) continue;
            arf_abs(x1, arb_midref(B(k, j)));
            if(arf_cmp(row_max[k], x1) < 0) {
                arf_set(row_max[k], x1);
                row_max_indices[k] = j;
            }
        }
    }

    for(int k = 0; k < dim; k++) {
        arb_set(D(k), B(k,k));
        arb_set_exact(D(k));
    }

    // At this point we've done that diagonalization and all that remains is
    // to certify the correctness and compute error bounds.

    arb_mat_t e;

    arb_t error_norms[dim];
    for(int k = 0; k < dim; k++) arb_init(error_norms[k]);

    arb_mat_init(e, dim, 1);

    arb_t z1, z2;
    arb_init(z1);
    arb_init(z2);
    for(int j = 0; j < dim; j++) {
        arb_mat_set(B, A);
        for(int k = 0; k < dim; k++) {
            arb_sub(B(k, k), B(k, k), D(j), prec);
        }
        for(int k = 0; k < dim; k++) {
            arb_set(arb_mat_entry(e, k, 0), P(k, j));
        }
        arb_mat_L2norm(z2, e, prec);
        arb_mat_mul(e, B, e, prec);
        arb_mat_L2norm(error_norms[j], e, prec);

        arb_div(z2, error_norms[j], z2, prec); // and now z1 is an upper bound for the
                                               // error in the eigenvalue
        arb_add_error(D(j), z2);
    }

    int unique_eigenvalues = 1;
    for(int j = 0; j < dim; j++) {
        if(j == 0) {
            arb_sub(z1, D(j), D(1), prec);
        }
        else {
            arb_sub(z1, D(j), D(0), prec);
        }
        arb_get_abs_lbound_arf(x1, z1, prec);
        for(int k = 1; k < dim; k++) {
            if(k == j) continue;
            arb_sub(z1, D(j), D(k), prec);
            arb_get_abs_lbound_arf(x2, z1, prec);
            if(arf_cmp(x2, x1) < 0) {
                arf_set(x1, x2);
            }
        }
        if(arf_is_zero(x1)) {
            unique_eigenvalues = 0;
        }
        arb_div_arf(z1, error_norms[j], x1, prec);
        for(int k = 0; k < dim; k++) {
            arb_add_error(P(k, j), z1);
        }
    }

    arb_mat_clear(e);
    arb_clear(z1);
    arb_clear(z2);
    for(int k = 0; k < dim; k++) arb_clear(error_norms[k]);

    arf_clear(x1);
    arf_clear(x2);
    arb_mat_clear(B);
    for(int k = 0; k < dim; k++) {
        arf_clear(B1[k]);
        arf_clear(B2[k]);
    }
    for(int k = 0; k < dim - 1; k++) {
        arf_clear(row_max[k]);
    }
    arf_clear(Gii);
    arf_clear(Gij);
    arf_clear(Gji);
    arf_clear(Gjj);
    free(B1);
    free(B2);
    free(row_max);
    free(row_max_indices);

    if(unique_eigenvalues) return 0;
    else return 1;
#undef B
#undef D
#undef P
}

int arb_mat_generalized_eigenproblem_symmetric_positive_definite(arb_mat_t D, arb_mat_t P, const arb_mat_t A, const arb_mat_t B, slong prec) {
    // solve the generalized eigenvalue problem Ax = lambda * Bx, where B
    // is symmetric positive definite and A is symmetric.
    //
    // Returns 0 on success.
    //
    // If there is a problem inverting the Cholesky factorization of B, returns -1.
    // If the eigenvalues are not distinct (or cannot be certified as distinct),
    // returns 1. In this case the returned eigenvalues will be correct, but
    // some eigenvectors will be nonsense, or will be accurate eigenvectors but with
    // infinite radius, or something like that.

    int dim = arb_mat_nrows(B);
    arb_mat_t L;
    arb_mat_t X;
    arb_mat_init(L, dim, dim);
    arb_mat_init(X, dim, dim);

    arb_mat_cholesky(L, B, prec);
    int result = arb_mat_inv(L, L, prec);
    if(!result) {
        arb_mat_clear(L);
        arb_mat_clear(X);
        return -1;
    }
    arb_mat_mul(X, L, A, prec);
    arb_mat_transpose(L, L);
    arb_mat_mul(X, X, L, prec);
    result = arb_mat_jacobi(D, P, X, prec);
    arb_mat_mul(P, L, P, prec);
    //arb_mat_transpose(P, P);

    arb_mat_clear(L);
    arb_mat_clear(X);
    return result;
}
