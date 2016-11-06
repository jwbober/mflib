#include <stdlib.h>
#include <stdint.h>
#include "acb.h"

#include "mfformat.h"

int main(int argc, char ** argv) {
    FILE * infile = fopen(argv[1], "r");

    struct mfheader header;
    if(!read_mfheader(infile, &header)) {
        printf("error reading header.\n");
        return -1;
    }

    printf("%d.%d.%d.%d.%d\n", header.level, header.weight, header.chi, header.orbit, header.j);
    fmpz_t x;
    fmpz_t y;
    fmpz_init(x);
    fmpz_init(y);
    acb_t z;
    acb_init(z);

    acb_t a;
    acb_init(a);

    for(unsigned int k = 0; k < header.ncoeffs; k++) {
        fmpz_inp_raw(x, infile);
        if(header.chi != 1) {
            fmpz_inp_raw(y, infile);
        }
        acb_set_mfcoeff(z, x, y, &header);
        acb_printd(z, 20);
        printf("\n");
    }

    fmpz_clear(x);
    fmpz_clear(y);
}
