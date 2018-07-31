#include <stdint.h>
#include "flint/fmpz_poly.h"

void x_fmpz_poly_write_raw(char ** data, size_t * datasize, fmpz_poly_t f) {
    // Allocate and fill memory for a portable binary representation of f.
    // On output, *data contains the representation, which takes *datasize bytes.
    // The calling function should free the allocated data (with a call to free())
    // when finished with it.
    FILE * pseudofile = open_memstream(data, datasize);

    int64_t degree = fmpz_poly_degree(f);
    fwrite((void*)&degree, sizeof(degree), 1, pseudofile);
    fmpz_t x;
    fmpz_init(x);
    for(int k = degree; k >= 0; k--) {
        fmpz_poly_get_coeff_fmpz(x, f, k);
        fmpz_out_raw(pseudofile, x);
    }
    fmpz_clear(x);
    fclose(pseudofile);
}

void x_fmpz_poly_read_raw(fmpz_poly_t f, const void * data, size_t datasize) {
    FILE * pseudofile = fmemopen( (void *)data, datasize, "r");
    int64_t degree;
    size_t result = fread((void*)&degree, sizeof(degree), 1, pseudofile);
    if(result == 0) return;
    fmpz_t x;
    fmpz_init(x);
    fmpz_poly_zero(f);
    for(slong k = degree; k >= 0; k--) {
        fmpz_inp_raw(x, pseudofile);
        fmpz_poly_set_coeff_fmpz(f, k, x);
    }
    fmpz_clear(x);
    fclose(pseudofile);
}


