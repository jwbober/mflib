#include <stdlib.h>
#include <time.h>
#include "arb_extras.h"

int test_cholesky(flint_rand_t state, int dim, int prec) {
    arb_mat_t mat;
    arb_mat_init(mat, dim, dim);
    for(int j = 0; j < dim; j++) {
        for(int k = 0; k <= j; k++) {
            arb_set_ui(arb_mat_entry(mat, j, k), 0);
            arf_randtest(arb_midref(arb_mat_entry(mat, j, k)), state, prec, 3);
            arf_abs(arb_midref(arb_mat_entry(mat, j, k)), arb_midref(arb_mat_entry(mat, j, k)));
            arb_set(arb_mat_entry(mat, k,j), arb_mat_entry(mat, j, k));
        }
    }
    arb_mat_t L, U;
    arb_mat_init(L, dim, dim);
    arb_mat_init(U, dim, dim);

    arb_mat_transpose(L, mat);
    arb_mat_mul(mat, L, mat, prec);
    arb_mat_zero(L);

    arb_mat_cholesky(L, mat, prec);
    arb_mat_transpose(U, L);
    arb_mat_mul(L, L, U, prec);

    if(arb_mat_contains(L, mat)) {
        return 1;
    }
    else {
        return 0;
    }

    arb_mat_clear(L);
    arb_mat_clear(U);
    arb_mat_clear(mat);
}

int test_jacobi(flint_rand_t state, int dim, int prec) {
    arb_mat_t A;
    arb_mat_t D;
    arb_mat_t P;
    arb_mat_t Pt;
    arb_mat_init(A, dim, dim);
    arb_mat_init(D, dim, dim);
    arb_mat_init(P, dim, dim);
    arb_mat_init(Pt, dim, dim);
    for(int j = 0; j < dim; j++) {
        for(int k = 0; k <= j; k++) {
            arf_randtest(arb_midref(arb_mat_entry(A, j, k)), state, prec, 3);
            arf_abs(arb_midref(arb_mat_entry(A, j, k)), arb_midref(arb_mat_entry(A, j, k)));
            arb_set(arb_mat_entry(A, k,j), arb_mat_entry(A, j, k));
        }
    }
    arb_mat_print_sage_float(A);
    printf("\n");
    int result = arb_mat_jacobi(D, P, A, prec);
    arb_mat_printd(D, 10);
    printf("\n");
    arb_mat_printd(P, 10);

    arb_mat_transpose(Pt, P);
    arb_mat_mul(A, Pt, A, prec);
    arb_mat_mul(A, A, P, prec);
    printf("\n");
    arb_mat_printd(A, 10);
    printf("%d\n", result);
    return 1;
}

int test_jacobi2() {
#define A(i,j) arb_mat_entry(A, i, j)
    arb_mat_t A;
    arb_mat_init(A, 5, 5);
    arb_set_si(A(0,0), 1);
    arb_set_si(A(0,1), -1);

    arb_set_si(A(1,0), -1);
    arb_set_si(A(1,1),  2);
    arb_set_si(A(1,2), -1);

    arb_set_si(A(2,1), -1);
    arb_set_si(A(2,2),  2);
    arb_set_si(A(2,3), -1);

    arb_set_si(A(3,2), -1);
    arb_set_si(A(3,3),  2);
    arb_set_si(A(3,4), -1);

    arb_set_si(A(4,3), -1);
    arb_set_si(A(4,4),  1);
#undef A

    arb_mat_scalar_mul_2exp_si(A, A, -8);

    arb_mat_t P;
    arb_mat_t D;

    arb_mat_init(P, 5, 5);
    arb_mat_init(D, 5, 1);

    arb_mat_jacobi(D, P, A, 4);
    arb_mat_printd(A, 15);
    printf("\n");
    arb_mat_printd(P, 15);
    printf("\n");
    arb_mat_printd(D, 15);
}

int main(int argc, char ** argv) {
    int dim = atoi(argv[1]);
    int prec = atoi(argv[2]);
    long seed = time(NULL);
    if(argc > 3)
        seed = atol(argv[3]);

    flint_rand_t state;
    flint_randinit(state);
    flint_randseed(state, seed, seed);

    //test_jacobi(state, dim, prec);
    test_jacobi2();
    return 0;
}
