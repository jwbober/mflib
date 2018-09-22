#include <stdio.h>
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"

void x_fmpz_vec_read_raw(long * len, fmpz ** vec, const void * data, size_t datasize) {
    FILE * pseudofile = fmemopen( (void*)data, datasize, "r");
    size_t result = fread((void*)len, sizeof(*len), 1, pseudofile);
    if(result == 0) {
        *len = -1;
        return;
    }
    *vec = _fmpz_vec_init(*len);
    for(long k = 0; k < *len; k++) {
        fmpz_inp_raw((*vec) + k, pseudofile);
    }
    fclose(pseudofile);
}
