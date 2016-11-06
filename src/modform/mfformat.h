#ifndef _MFFORMAT_H_
#define _MFFORMAT_H_

#ifdef __cplusplus
#include <cstdint>
#include <cstdio>
using std::fwrite;
using std::fread;
using std::FILE;
#else
#include <stdint.h>
#include <stdio.h>
#endif

#include "acb.h"

const long MF_PREC_EXACT = 2147483647;

struct mfheader {
    uint32_t version; // == 3229261 (ascii string "MF1\0")
    uint32_t level;
    uint32_t weight;
    uint32_t chi;
    uint32_t orbit;
    uint32_t j;
    int32_t prec;
    int32_t exponent;
    uint32_t ncoeffs;
    char reserved[92];
}; // For a grand total of 128 bytes, since that is a nice round number.

static int write_mfheader(FILE * outfile, struct mfheader * header) {
    if(!fwrite((char*)&header->version, sizeof(header->version), 1, outfile)) return 0;
    if(!fwrite((char*)&header->level, sizeof(header->level), 1, outfile)) return 0;
    if(!fwrite((char*)&header->weight, sizeof(header->weight), 1, outfile)) return 0;
    if(!fwrite((char*)&header->chi, sizeof(header->chi), 1, outfile)) return 0;
    if(!fwrite((char*)&header->orbit, sizeof(header->orbit), 1, outfile)) return 0;
    if(!fwrite((char*)&header->j, sizeof(header->j), 1, outfile)) return 0;
    if(!fwrite((char*)&header->prec, sizeof(header->prec), 1, outfile)) return 0;
    if(!fwrite((char*)&header->exponent, sizeof(header->exponent), 1, outfile)) return 0;
    if(!fwrite((char*)&header->ncoeffs, sizeof(header->ncoeffs), 1, outfile)) return 0;
    if(!fwrite((char*)&header->reserved, sizeof(header->reserved), 1, outfile)) return 0;
    return 1;
}

static int read_mfheader(FILE * outfile, struct mfheader * header) {
    if(!fread((char*)&header->version, sizeof(header->version), 1, outfile)) return 0;
    if(header->version != 3229261) return 0;
    if(!fread((char*)&header->level, sizeof(header->level), 1, outfile)) return 0;
    if(!fread((char*)&header->weight, sizeof(header->weight), 1, outfile)) return 0;
    if(!fread((char*)&header->chi, sizeof(header->chi), 1, outfile)) return 0;
    if(!fread((char*)&header->orbit, sizeof(header->orbit), 1, outfile)) return 0;
    if(!fread((char*)&header->j, sizeof(header->j), 1, outfile)) return 0;
    if(!fread((char*)&header->prec, sizeof(header->prec), 1, outfile)) return 0;
    if(!fread((char*)&header->exponent, sizeof(header->exponent), 1, outfile)) return 0;
    if(!fread((char*)&header->ncoeffs, sizeof(header->ncoeffs), 1, outfile)) return 0;
    if(!fread((char*)&header->reserved, sizeof(header->reserved), 1, outfile)) return 0;
    return 1;
}

static void acb_set_mfcoeff(acb_t out, fmpz_t x, fmpz_t y, struct mfheader * header) {
    acb_set_fmpz_fmpz(out, x, y);
    acb_mul_2exp_si(out, out, header->exponent);
    if(header->prec != MF_PREC_EXACT) {
        mag_set_ui(arb_radref(acb_realref(out)), 1);
        mag_mul_2exp_si(arb_radref(acb_realref(out)), arb_radref(acb_realref(out)), header->prec);
        if(header->chi != 1) {
            mag_set_ui(arb_radref(acb_imagref(out)), 1);
            mag_mul_2exp_si(arb_radref(acb_imagref(out)), arb_radref(acb_imagref(out)), header->prec);
        }
    }
}

#endif
