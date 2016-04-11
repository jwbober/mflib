#include <stdlib.h>
#include "acb.h"

int main(int argc, char ** argv) {
    FILE * infile = fopen(argv[1], "r");

    int prec;
    fread( &prec, sizeof(prec), 1, infile);
    fmpz_t x;
    fmpz_t y;
    fmpz_init(x);
    fmpz_init(y);
    acb_t z;
    acb_init(z);

    acb_t a;
    acb_init(a);

    for(int k = 0; k < 4000; k++) {
        fmpz_inp_raw(x, infile);
        fmpz_inp_raw(y, infile);
        acb_set_fmpz_fmpz(z, x, y);
        mag_set_ui(arb_radref(acb_realref(z)), 1);
        mag_set_ui(arb_radref(acb_imagref(z)), 1);
        acb_mul_2exp_si(z, z, -prec);
        acb_printd(z, 10);
        printf("\n");
    }
    //acb_printd(a, 10);

    fmpz_clear(x);
    fmpz_clear(y);
}
