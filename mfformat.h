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
#include "sqlite3.h"

#ifdef __cplusplus
#define MF_PREC_EXACT 2147483647
#define MFV2 3294797
#define MFV1 3229261
extern "C" {
#else
const int MF_PREC_EXACT = 2147483647;
const int MFV2 = 3294797;
const int MFV1 = 3229261;
#endif

struct mfheader {
    uint32_t version; // == MFV1 or MFV2
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

int write_mfheader(FILE * outfile, struct mfheader * header);       // read or write an mfheader.
int read_mfheader(FILE * outfile, struct mfheader * header);        // return 0 on failue, 1 on success.
                                                                    // NOTE: when 0 is returned, some bytes
                                                                    // might have been read from the file.

void acb_set_mfcoeff(acb_t out, fmpz_t x, fmpz_t y, struct mfheader * header);  // Set out to the acb_t represented
                                                                                // by x and y, according to the
                                                                                // precision and exponent in header.

size_t acb_write_mfcoeff(FILE * outfile, struct mfheader * header, acb_t coeff);
    // write one or two fmpz_ts to outfile (dependin on whether or not chi is trivial)
    // representing coeff, according to precision and exponent in header. Return number
    // of bytes written. (Can fail if coeff is not precise enough, or possibly because
    // of other flakiness. For example, it is probably not possible currently to call acb_write_mfcoeff
    // on a coefficient that was read with acb_read_mfcoeff.)
    //
    // Doesn't handle failure in the case where one cal to fmpz_out_raw succeeds and the next fails,
    // and maybe doesn't handle other failure well. (If fmpz_out_raw returns nonzero, does that mean
    // it succeeded?)

int acb_attempt_rounding_mfcoeff(struct mfheader * header, acb_t coeff);
    // attempts to round coeff to the precision spefified in header. If this
    // succeeds, then a call to acb_write_mfcoeff should succeed as well, unless
    // there is an IO error.

int read_mffile(FILE * infile, struct mfheader * header, acb_ptr * coeffs);
    // populate header and *coeffs with data from infile.
    // on return, *header will have header->ncoeffs entries and should subsequently
    // be freed with _acb_vec_free(*coeffs, header->ncoeffs);

int write_mffile(FILE * outfile, struct mfheader * header, acb_ptr coeffs);
    // write mffile to outfile. expects header->version == MFV2, and otherwise
    // will return 0.
    // Returns 0 on failure, 1 on success. (Like with other functions, failure does
    // not mean that 0 bytes have been written.)

int read_mfdata(FILE * infile, struct mfheader * header, acb_ptr * coeffs);
    // like read_mffile, but doesn't read the header.

int write_mfdata(FILE * outfile, struct mfheader * header, acb_ptr coeffs);
    // like write_mffile, but doesn't write the header

int insert_into_sqlite(sqlite3 * db, struct mfheader * header, acb_ptr coeffs);

void iterate_through_sqlitefile(
        const char * filename,
        int (*callback)(struct mfheader * header, int coeff_datasize, const void * coeff_data, acb_ptr coeffs),
        int populate_coefficients);
    // go through the entire sqlite database filename, and for each entry populate header and coeffs and then call callback()
    // the callback() should NOT free the coeffcients or the header
    // If callback() ever returns nonzero, we stop the iteration.
    // If populate_coefficients is nonzero, we convert the raw coefficient data and store it in coeffs, otherwise we just give the raw data (which is useful for writing out files, for example)

#ifdef __cplusplus
} // extern "C"
#endif

#endif
