#include <stdio.h>
#include "flint/fmpz.h"

void x_fmpz_vec_write_raw(char ** data, size_t * datasize, long len, fmpz * vec) {
    FILE * pseudofile = open_memstream(data, datasize);

    fwrite((void*)&len, sizeof(len), 1, pseudofile);
    for(long k = 0; k < len; k++) {
        fmpz_out_raw(pseudofile, vec + k);
    }
    fclose(pseudofile);
}
