#include "flint/nmod_mat.h"

void nmod_mat_print_sage(const nmod_mat_t A, const char * name) {
    int nrows = nmod_mat_nrows(A);
    int ncols = nmod_mat_ncols(A);
    printf("%s = Matrix([", name);
    for(int j = 0; j < nrows; j++) {
        printf("[");
        for(int k = 0; k < ncols; k++) {
            printf("%ld", nmod_mat_entry(A, j, k));
            if(k < nrows - 1)
                printf(", ");
        }
        printf("],\n");
    }
    printf("], ring = Integers(%ld))\n", A[0].mod.n);
}
