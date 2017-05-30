//
// File: mpfi_fft.h
//
// 8th October 2013
//
// Dave Platt, University of Bristol
//
//
// 'C' code to augment MPFI with complex intervals
// and to provide basic dft functionality
//
// Tries to keep to a sensible naming scheme
// all mpfi_c routines are mpfi_c_xxx where xxx is similar to MPFI
// 
// All routines allow one or more arguments to be the same,
// e.g. mpfi_c_add(z,z,z) should double z.
// but it is probably more efficient to avoid this
//
// Unpredicatble results will happen if the real and imaginary parts of
// an mpfi_c_t are actually the same mpfi_t (why would you do that?)
//
// To do:
// Get rid of remaining op(op1,op1,op2) usage
// 
#ifndef MPFI_FFT
#define MPFI_FFT
#include "mpfi.h"
#include "inttypes.h"
#include "stdlib.h"
#include "stdio.h"

#ifdef __cplusplus
extern "C" {
#endif

// First the mpfi_c stuff

/* define a complex version of an MPFI interval */
typedef struct{
    mpfi_t re;
    mpfi_t im;
} _mpfi_c_struct;

// structs to make parameter passing easier
typedef _mpfi_c_struct mpfi_c_t[1];
typedef _mpfi_c_struct *mpfi_c_ptr;

/* initialisation */
static inline void mpfi_c_init(mpfi_c_ptr z)
{
    mpfi_init(z->re);
    mpfi_init(z->im);
};

/* clearing */
static inline void mpfi_c_clear(mpfi_c_ptr z)
{
    mpfi_clear(z->re);
    mpfi_clear(z->im);
};

/* efficient swapping */
static inline void mpfi_c_swap(mpfi_c_ptr z1,mpfi_c_ptr z2)
{
    mpfi_swap(z1->re,z2->re);
    mpfi_swap(z1->im,z2->im);
};

/* simple printing */
void mpfi_c_print(mpfi_c_ptr z);
void mpfi_c_print_str(const char *str, mpfi_c_ptr x);
void mpfi_c_printn(mpfi_c_ptr z,int n);
void mpfi_print(mpfi_ptr x);
void mpfi_print_str(const char *str, mpfi_ptr x);
void mpfi_printn(mpfi_ptr x, int n);

// setting from various data types
static inline void mpfi_c_set(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_set(z1->re,z2->re);
    mpfi_set(z1->im,z2->im);
};

static inline void mpfi_c_set_d(mpfi_c_ptr z, double re, double im)
{
  mpfi_set_d(z->re,re);
  mpfi_set_d(z->im,im);
};

static inline void mpfi_c_set_ui(mpfi_c_ptr z, unsigned long int re, unsigned long int im)
{
  mpfi_set_ui(z->re,re);
  mpfi_set_ui(z->im,im);
};

// just set the real part
static inline void mpfi_c_set_re(mpfi_c_ptr z,mpfi_ptr x)
{
    mpfi_set(z->re,x);
};

static inline void mpfi_c_set_re_d(mpfi_c_ptr z,double x)
{
    mpfi_set_d(z->re,x);
};

// just set the imaginary part
static inline void mpfi_c_set_im(mpfi_c_ptr z, mpfi_ptr x)
{
    mpfi_set(z->im,x);
};

// z1=z2+z3
static inline void mpfi_c_add(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_add(z1->re,z2->re,z3->re);
    mpfi_add(z1->im,z2->im,z3->im);
};

// z1=z2-z3
static inline void mpfi_c_sub(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_sub(z1->re,z2->re,z3->re);
    mpfi_sub(z1->im,z2->im,z3->im);
};

// z1=-z2
static inline void mpfi_c_neg(mpfi_c_ptr z1,mpfi_c_ptr z2)
{
    mpfi_neg(z1->im,z2->im);
    mpfi_neg(z1->re,z2->re);
};

// res=|z|^2
void mpfi_c_norm(mpfi_ptr res, mpfi_c_ptr z);
// res=|z|
void mpfi_c_abs(mpfi_ptr res, mpfi_c_ptr z);

/* use the obvious method, anything cleverer loses precision */
void mpfi_c_mul(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3);
// better than mpfi_c_mul(z1,z2,z2)
void mpfi_c_sqr(mpfi_c_ptr z1, mpfi_c_ptr z2);
/* multiply (mpfi_c_t z2) by (mpfi_t x) into (mpfi_c_t z1) safely */

static inline void mpfi_c_mul_i(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_ptr x)
{
    mpfi_mul(z1->re,z2->re,x);
    mpfi_mul(z1->im,z2->im,x);
};

// scalar mult
static inline void mpfi_c_mul_d(mpfi_c_ptr z1, mpfi_c_ptr z2, double x)
{
    mpfi_mul_d(z1->re,z2->re,x);
    mpfi_mul_d(z1->im,z2->im,x);
}

// and again
static inline void mpfi_c_mul_ui(mpfi_c_ptr z1, mpfi_c_ptr z2, unsigned long int x)
{
    mpfi_mul_ui(z1->re,z2->re,x);
    mpfi_mul_ui(z1->im,z2->im,x);
}

/* multiply by i in situ */
static inline void mpfi_c_muli(mpfi_c_ptr z)   
{
    mpfi_swap(z->re,z->im);
    mpfi_neg(z->re,z->re);
};

/* multiply (mpfi_c_t z2) by (mpfr_t x) into (mpfi_c_t z1) safely */
static inline void mpfi_c_mul_fr(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfr_ptr x)
{
    mpfi_mul_fr(z1->re,z2->re,x);
    mpfi_mul_fr(z1->im,z2->im,x);
};

// scalar mult by mpz_t
static inline void mpfi_c_mul_z(mpfi_c_ptr z1,mpfi_c_ptr z2, mpz_ptr x)
{
    mpfi_mul_z(z1->re,z2->re,x);
    mpfi_mul_z(z1->im,z2->im,x);
};
    

/* z1<-conjugate(z2) safely */
static inline void mpfi_c_conj(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_c_set(z1,z2);
    mpfi_neg(z1->im,z1->im);
};

// z1=z2/z3
void mpfi_c_div(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3);

// scalar division by mpz
static inline void mpfi_c_div_z(mpfi_c_ptr z1,mpfi_c_ptr z2, mpz_ptr x)
{
    mpfi_div_z(z1->re,z2->re,x);
    mpfi_div_z(z1->im,z2->im,x);
}

// scalar division by mpfi_t
static inline void mpfi_c_div_i (mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_ptr x)
{
	mpfi_div(z1->re,z1->re,x);
	mpfi_div(z1->im,z1->im,x);

}

// scalar division by ui
static inline void mpfi_c_div_ui (mpfi_c_ptr z1, mpfi_c_ptr z2, unsigned long int i)
{
    mpfi_div_ui(z1->re,z2->re,i);
    mpfi_div_ui(z1->im,z2->im,i);
};

// scalar division by double
static inline void mpfi_c_div_d(mpfi_c_ptr z1, mpfi_c_ptr z2, double x)
{
  mpfi_div_d(z1->re,z2->re,x);
  mpfi_div_d(z1->im,z2->im,x);
}

// z1=exp(z2)
void mpfi_c_exp(mpfi_c_ptr z1, mpfi_c_ptr z2);

// res=x^(re_s+I im_s)
void mpfi_c_pow_double_to_doubles(mpfi_c_ptr res, double x, double re_s, double im_s);

// res=x^y in reals
void mpfi_pow(mpfi_ptr res, mpfi_ptr x, mpfi_ptr y);

// res=x^z Complex=Real^Complex
void mpfi_c_pow_i_c (mpfi_c_ptr res, mpfi_ptr x, mpfi_c_ptr z);

// this appear to have got defined in mpfi somewhere along the way
/*
inline void mpfi_atan2 (mpfi_ptr res, mpfi_ptr y, mpfi_ptr x)
{
  mpfi_div(res,y,x);
  mpfi_atan(res,res);
  if(mpfi_is_neg(x))
    {
      if(mpfi_is_neg(y))
	mpfi_sub(res,res,mpfi_pi);
      else
	mpfi_add(res,res,mpfi_pi);
    }
}
*/

 // res=log(z)
void mpfi_c_log (mpfi_c_ptr res, mpfi_c_ptr z);

// aka mpfi_c_abs
#define mpfi_c_mod(a,b) mpfi_c_abs(a,b)

// arg(z)
void mpfi_c_arg(mpfi_ptr res, mpfi_c_ptr z);

static inline int mpfi_contains_zero(mpfi_ptr x)
{
  return(mpfi_cmp_d(x,0.0)==0);
}

static inline int mpfi_c_contains_zero(mpfi_c_ptr z)
{
  return(mpfi_contains_zero(z->re)&&mpfi_contains_zero(z->im));
}

//z=0+0i
static inline void mpfi_c_zero(mpfi_c_ptr z)
{
   mpfi_set_ui(z->re,0);
   mpfi_set_ui(z->im,0);
}

//z=1+0i
static inline void mpfi_c_one(mpfi_c_ptr z)
{
  mpfi_set_ui(z->re,1);
  mpfi_set_ui(z->im,0);
}

// returns root in right half plane
void mpfi_c_sqrt(mpfi_c_ptr res, mpfi_c_ptr z);

//
// Now the DFT stuff
//

/* perform an in place FFT
   n a power of 2
   x contains n contiguous mpfi_c_t
   w contains n/2 values of exp(-2 pi a/n) a=0..n/2-1
*/
void fft(mpfi_c_t *x,unsigned int n, mpfi_c_t *w);

/* perform an in place inverse FFT */
void ifft(mpfi_c_t *x,unsigned int n, mpfi_c_t *w);

/* perform a non-normalised inverse FFT */
void nifft(mpfi_c_t *x,unsigned int n, mpfi_c_t *w);

/* circular convolution; f and g must be distinct */
void convolve(mpfi_c_t *result,mpfi_c_t *f,mpfi_c_t *g,unsigned int n, mpfi_c_t *w); 

// Do iDFT on a Hermitian vector of length 2N, i.e. one that produces real values
// N is the reduced length and the ws relate to that N
// ws[n]=e(n/N) for n=0..N/2-1
// x is of length 2N with 0..N input values and (N+1..2N-1) anything
// omega=exp(2 pi I/2/N)
// on exit x[n]->re are set, x[n]->im are garbage.
//
void hermidft(mpfi_c_t *x, int N, mpfi_c_t *ws, mpfi_c_ptr omega);

void dft_init(uint64_t prec);

// return smallest power of 2 > 2n-1 
uint64_t calc_conv_size(uint64_t n);

// use Bluestein's algorithm to do a length len DFT on x
// where each datum is stride apart
void bluesteinfft(mpfi_c_t *x, uint64_t len, uint64_t stride);

// do a simple dft by O(n^2) method
// input in x[0],x[stride].., output to x 
void simple_dft(mpfi_c_t *x, unsigned long int N, unsigned long int stride);
void dft(mpfi_c_t *x, unsigned long int len, unsigned long int stride);

//
// do an ndft
// total length of x is N made up of num_dims dimensions whose lengths are
// specified in dims
//
void ndft(mpfi_c_t *x, unsigned long int N, unsigned long int num_dims, unsigned long int *dims);

#ifdef __cplusplus
}
#endif

#endif
