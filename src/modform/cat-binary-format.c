#include <stdlib.h>
#include <stdint.h>
#include "acb.h"

struct mfheader {
    uint32_t level;
    uint32_t weight;
    uint32_t chi;
    uint32_t orbit;
    uint32_t j;
    uint32_t prec;
    uint32_t ncoeffs;
};

int write_mfheader(FILE * outfile, struct mfheader * header) {
    if(!fwrite((char*)&header->level, sizeof(header->level), 1, outfile)) return 0;
    if(!fwrite((char*)&header->weight, sizeof(header->weight), 1, outfile)) return 0;
    if(!fwrite((char*)&header->chi, sizeof(header->chi), 1, outfile)) return 0;
    if(!fwrite((char*)&header->orbit, sizeof(header->orbit), 1, outfile)) return 0;
    if(!fwrite((char*)&header->j, sizeof(header->j), 1, outfile)) return 0;
    if(!fwrite((char*)&header->prec, sizeof(header->prec), 1, outfile)) return 0;
    if(!fwrite((char*)&header->ncoeffs, sizeof(header->ncoeffs), 1, outfile)) return 0;
    return 1;
}

int read_mfheader(FILE * infile, struct mfheader * header) {
    if(!fread((char*)&header->level, sizeof(header->level), 1, infile)) return 0;
    if(!fread((char*)&header->weight, sizeof(header->weight), 1, infile)) return 0;
    if(!fread((char*)&header->chi, sizeof(header->chi), 1, infile)) return 0;
    if(!fread((char*)&header->orbit, sizeof(header->orbit), 1, infile)) return 0;
    if(!fread((char*)&header->j, sizeof(header->j), 1, infile)) return 0;
    if(!fread((char*)&header->prec, sizeof(header->prec), 1, infile)) return 0;
    if(!fread((char*)&header->ncoeffs, sizeof(header->ncoeffs), 1, infile)) return 0;
    return 1;
}

int main(int argc, char ** argv) {
    FILE * infile = fopen(argv[1], "r");

    unsigned int magic;
    fread( &magic, sizeof(magic), 1, infile);
    if(magic != 3229261) {
        printf("wrong file format or error reading file.\n");
        return -1;
    }

    struct mfheader header;
    read_mfheader(infile, &header);

    printf("%d.%d.%d.%d.%d\n", header.level, header.weight, header.chi, header.orbit, header.j);
    fmpz_t x;
    fmpz_t y;
    fmpz_init(x);
    fmpz_init(y);
    acb_t z;
    acb_init(z);

    acb_t a;
    acb_init(a);

    for(int k = 0; k < header.ncoeffs; k++) {
        fmpz_inp_raw(x, infile);
        if(header.chi != 1)
            fmpz_inp_raw(y, infile);
        acb_set_fmpz_fmpz(z, x, y);
        mag_set_ui(arb_radref(acb_realref(z)), 1);
        if(header.chi != 1)
            mag_set_ui(arb_radref(acb_imagref(z)), 1);
        acb_mul_2exp_si(z, z, -(long)header.prec);
        acb_printd(z, 20);
        printf("\n");
    }
    //acb_printd(a, 10);

    fmpz_clear(x);
    fmpz_clear(y);
}
