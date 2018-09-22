#ifndef _FLINT_EXTRAS_H_
#define _FLINT_EXTRAS_H_

// some convenience functions for flint types.
//
// note: we don't include any flint headers here, so those need to be included
// first. maybe that's not good, but I don't want to worry about whether it is
// #include "flint/whatever" or #include "whatever"

#ifdef __cplusplus

#include <ostream>
std::ostream& operator << (std::ostream& out, fmpz_t in);
std::ostream& operator << (std::ostream& out, fmpz_poly_t in);

extern "C" {
#endif

void x_fmpz_poly_write_raw(char ** data, size_t * datasize, fmpz_poly_t f);
// Allocate and fill memory for a portable binary representation of f.
// On output, *data contains the representation, which takes *datasize bytes.
// The calling function should free the allocated data (with a call to free())
// when finished with it.
// (The point of this is that the resulting output can be parsed by x_fmpz_read_raw().
//
// The "format", such as it is, is an int64_t (the degree of the polynomial)
// followed by degree + 1 calls to fmpz_out_raw to write out the coefficients
// from the highest degree downwards.
//
// (Technically, the degree of an fmpz appears to be a flint slong, which is
// "guaranteed to be the same size as GMP's mp_limb_signed_t". One day if
// this is 128 bits and someone tries to write out a >64 bit polynomial, there
// are going to be problems.)

void x_fmpz_poly_read_raw(fmpz_poly_t f, const void * data, size_t datasize);
// Set the polynomial f to the polynomial represented by data, which is a blob
// of datasize bytes. (Which was probably written out by x_fmp_poly_write_raw().)
//
// If datasize = 0, will probably do nothing at all, I think.

// The following two functions are basically the same as the fmpz_poly_t versions
// (as in, a vector could be read as a polynomial and a polynomial could be read
// as a vector).
void x_fmpz_vec_write_raw(char ** data, size_t * datasize, long len, fmpz * vec);
void x_fmpz_vec_read_raw(long * len, fmpz ** vec, const void * data, size_t datasize);


#ifdef __cplusplus
}
#endif

#endif // _FLINT_EXTRAS_H_
