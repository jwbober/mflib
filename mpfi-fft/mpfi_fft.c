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

#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "math.h"
#include "mpfi.h"
#include "mpfi_fft.h"
// First the mpfi_c stuff

/* simple printing */
void mpfi_c_print(mpfi_c_ptr z)
{
    mpfi_out_str(stdout,10,0,z->re);
    printf(" + ");
    mpfi_out_str(stdout,10,0,z->im);
    printf("i\n");
};

void mpfi_c_print_str(const char *str, mpfi_c_ptr x)
{
  printf("%s",str);
  mpfi_c_print(x);
}

void mpfi_c_printn(mpfi_c_ptr z,int n)
{
    mpfi_out_str(stdout,10,n,z->re);
    printf(" + ");
    mpfi_out_str(stdout,10,n,z->im);
    printf("i\n");
};

void mpfi_print(mpfi_ptr x)
{
    mpfi_out_str(stdout,10,0,x);
    printf("\n");
};

void mpfi_print_str(const char *str, mpfi_ptr x)
{
  printf("%s",str);
  mpfi_print(x);
}

void mpfi_printn(mpfi_ptr x, int n)
{
    mpfi_out_str(stdout,10,n,x);
    printf("\n");
};

// globally defined variables
// used to avoid repeated init/clear
//
mpfi_c_t d_tmp3,d_tmp5,c_sqrt_tmp;
mpfi_t m_tmp1,m_tmp2,m_tmp3,d_tmp1,d_tmp2,d_tmp4,e_tmp,p_tmp1,p_tmp2;   
mpfi_t norm_tmp,norm_tmp1,c_log_tmp,c_log_tmp1,c_arg_tmp;
mpfi_t mod_tmp,mpfi_pi,mpfi_2_pi,mpfi_ln_2_pi,mpfi_sqrt_pi,mpfi_ln_pi;
mpfi_c_t lng_z,lng_tmp,lng_tmp1,lng_tmp2,ln_gamma_err,ln_gamma_err1;
mpfi_t log_2_pi_over_2;
mpfi_c_t c_one,c_minus_one,c_i,c_minus_i,c_zero;

mpfi_t mpfi_pi_by_2,mpfi_pi_by_4;
mpfi_t mpfi_log_2;


mpfi_t pp_x_sigma,pp_tlnx,pow_tmp,pow_tmp1;


void mpfi_c_setup(int prec)
{
/* now initialise all global variables.
   we don't bother clearing them down afterwards, just exit. */
  mpfr_set_default_prec(prec);
  mpfi_init(m_tmp1);
  mpfi_init(m_tmp2);
  mpfi_init(m_tmp3);
  mpfi_init(d_tmp1);
  mpfi_init(d_tmp2);
  mpfi_init(d_tmp4);
  mpfi_c_init(d_tmp3);
  mpfi_c_init(d_tmp5);
  mpfi_init(e_tmp);
  mpfi_init(p_tmp1);
  mpfi_init(p_tmp2);
  mpfi_init(norm_tmp);
  mpfi_init(norm_tmp1);
  mpfi_init(mod_tmp);
  mpfi_init(mpfi_pi);
  mpfi_const_pi(mpfi_pi);
  mpfi_init(mpfi_sqrt_pi);
  mpfi_sqrt(mpfi_sqrt_pi,mpfi_pi);
  mpfi_init(mpfi_2_pi);
  mpfi_init(mpfi_ln_2_pi);
  mpfi_mul_ui(mpfi_2_pi,mpfi_pi,2);
  mpfi_log(mpfi_ln_2_pi,mpfi_2_pi);
  mpfi_init(log_2_pi_over_2);
  mpfi_set(log_2_pi_over_2,mpfi_ln_2_pi);
  mpfi_div_ui(log_2_pi_over_2,log_2_pi_over_2,2);
  mpfi_init(c_log_tmp);
  mpfi_init(c_log_tmp1);
  mpfi_init(c_arg_tmp);
  mpfi_c_init(c_sqrt_tmp);
  mpfi_init(pp_tlnx);
  mpfi_init(pp_x_sigma);
  mpfi_init(pow_tmp);
  mpfi_init(pow_tmp1);
  mpfi_c_init(c_one);
  mpfi_c_init(c_minus_one);
  mpfi_c_init(c_i);
  mpfi_c_init(c_minus_i);
  mpfi_c_init(c_zero);
  mpfi_c_set_ui(c_one,1,0);
  mpfi_c_neg(c_minus_one,c_one);
  mpfi_c_set_ui(c_i,0,1);
  mpfi_c_neg(c_minus_i,c_i);
  mpfi_c_set_ui(c_zero,0,0);
  mpfi_init(mpfi_pi_by_2);
  mpfi_div_2ui(mpfi_pi_by_2,mpfi_pi,1);
  mpfi_init(mpfi_pi_by_4);
  mpfi_div_2ui(mpfi_pi_by_4,mpfi_pi,2);
  mpfi_init(mpfi_log_2);
  mpfi_const_log2(mpfi_log_2);
  mpfi_init(mpfi_ln_pi);
  mpfi_log(mpfi_ln_pi,mpfi_pi);
}

// res=|z|^2
void mpfi_c_norm(mpfi_ptr res, mpfi_c_ptr z)
{
  mpfi_sqr(norm_tmp1,z->re);
  mpfi_sqr(norm_tmp,z->im);
  mpfi_add(res,norm_tmp1,norm_tmp);
}

// res=|z|
void mpfi_c_abs(mpfi_ptr res, mpfi_c_ptr z)
{
   mpfi_c_norm(res,z);
   mpfi_sqrt(res,res);
}


/* use the obvious method, anything cleverer loses precision */
void mpfi_c_mul(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_mul(m_tmp1,z2->re,z3->re);
    mpfi_mul(m_tmp2,z2->im,z3->im);
    mpfi_sub(m_tmp3,m_tmp1,m_tmp2);
    mpfi_mul(m_tmp1,z2->re,z3->im);
    mpfi_mul(m_tmp2,z2->im,z3->re);
    mpfi_add(z1->im,m_tmp1,m_tmp2);
    mpfi_swap(z1->re,m_tmp3); // for cache reasons, might actually be quicker to use mpfi_set
};

// better than mpfi_c_mul(z1,z2,z2)
void mpfi_c_sqr(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
  mpfi_sqr(m_tmp1,z2->re);
  mpfi_sqr(m_tmp2,z2->im);
  mpfi_mul(m_tmp3,z2->re,z2->im);
  mpfi_sub(z1->re,m_tmp1,m_tmp2);
  mpfi_mul_2ui(z1->im,m_tmp3,1);
}

// z1=z2/z3
void mpfi_c_div(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_sqr(d_tmp1,z3->re);
    mpfi_sqr(d_tmp2,z3->im);
    mpfi_add(d_tmp4,d_tmp1,d_tmp2);
    mpfi_c_conj(d_tmp3,z3);
    mpfi_c_mul(d_tmp5,d_tmp3,z2);
    mpfi_div(z1->re,d_tmp5->re,d_tmp4);
    mpfi_div(z1->im,d_tmp5->im,d_tmp4);
};

// z1=exp(z2)
void mpfi_c_exp(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_exp(e_tmp,z2->re);
    mpfi_cos(z1->re,z2->im);
    mpfi_sin(z1->im,z2->im);
    mpfi_c_mul_i(z1,z1,e_tmp);
};

// res=x^(re_s+I im_s)
void mpfi_c_pow_double_to_doubles(mpfi_c_ptr res, double x, double re_s, double im_s)
{
  mpfi_set_d(p_tmp1,x);     
  mpfi_log(p_tmp1,p_tmp1);        // log(x)
  mpfi_mul_d(p_tmp2,p_tmp1,re_s);  // sigma*log(x)
  mpfi_exp(p_tmp2,p_tmp2);        // x^sigma
  mpfi_mul_d(p_tmp1,p_tmp1,im_s);  // t*log(x)
  mpfi_cos(res->re,p_tmp1);       // Re(res)=cos(t*log(x))
  mpfi_sin(res->im,p_tmp1);       // Im(res)=sin(t*log(x))
  mpfi_c_mul_i(res,res,p_tmp2);   // res*=x^sigma
}

// res=x^y in reals
void mpfi_pow(mpfi_ptr res, mpfi_ptr x, mpfi_ptr y)
{
  mpfi_log(pow_tmp,x);
  mpfi_mul(pow_tmp1,pow_tmp,y);
  mpfi_exp(res,pow_tmp1);
}

// res=x^z Complex=Real^Complex
void mpfi_c_pow_i_c (mpfi_c_ptr res, mpfi_ptr x, mpfi_c_ptr z)
{
  mpfi_pow(pp_x_sigma,x,z->re);
  mpfi_log(pp_tlnx,x);
  mpfi_mul(pp_tlnx,pp_tlnx,z->im);
  mpfi_cos(res->re,pp_tlnx);
  mpfi_sin(res->im,pp_tlnx);
  mpfi_c_mul_i(res,res,pp_x_sigma);
}

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
void mpfi_c_log (mpfi_c_ptr res, mpfi_c_ptr z)
{
  mpfi_c_norm(c_log_tmp,z);
  mpfi_log(c_log_tmp1,c_log_tmp);
  mpfi_div_ui(c_log_tmp,c_log_tmp1,2);
  mpfi_atan2(res->im,z->im,z->re);
  mpfi_set(res->re,c_log_tmp);
}

// arg(z)
void mpfi_c_arg(mpfi_ptr res, mpfi_c_ptr z);

// returns root in right half plane
void mpfi_c_sqrt(mpfi_c_ptr res, mpfi_c_ptr z)
{
   mpfi_c_log(c_sqrt_tmp,z);
   mpfi_div_2ui(c_sqrt_tmp->re,c_sqrt_tmp->re,1);
   mpfi_div_2ui(c_sqrt_tmp->im,c_sqrt_tmp->im,1);
   mpfi_c_exp(res,c_sqrt_tmp);
   if(mpfi_is_neg(res->re))
     mpfi_c_neg(res,res);
}

//
// Now the DFT stuff
//

#define LOG2_MAX_FFT_LEN (20)
#define MAX_FFT_LEN (((uint64_t) 1)<<LOG2_MAX_FFT_LEN)
mpfi_c_t X[MAX_FFT_LEN];
#define MAX_SIMPLE_DFT (31)

#ifndef fatal_error
#define fatal_error(str) {fprintf(stderr,"%s Exiting.\n",str);exit(1);}
#endif

// a structure to hold the persistent data
typedef struct{
  uint64_t conv_size; // the power of two larger than 2n-1
  mpfi_c_t *b;        // the exp(2*pi*i*k^2/n) after the dft
  mpfi_c_t *b_conj;   // the conj of exp(2*pi*i*k^2/n)
} bluestein_struct;

bluestein_struct B[MAX_FFT_LEN];

mpfi_c_t ctemp,htemp,*ws[LOG2_MAX_FFT_LEN+1],A[MAX_FFT_LEN];

int64_t log_2(uint64_t n)
{
  if(__builtin_popcount(n)!=1)
    return(-1);

  int64_t res=0;
  while(n!=1)
    {
      res++;n>>=1;
    }
  return(res);
}

int power2(uint64_t n)
{
  return(__builtin_popcount(n)==1);
}

/* perform an in place FFT
   n a power of 2
   x contains n contiguous mpfi_c_t
   w contains n/2 values of exp(-2 pi a/n) a=0..n/2-1
*/
void fft(mpfi_c_t *x,unsigned int n, mpfi_c_t *w) {
	int i,j,k,l;
	mpfi_c_t *p,*xend=x+n;

	/* swap each element with one with bit-reversed index */
	for (i=0,l=n>>1;i<l;++i) {
		/* j = bit reversal of i */
		for (k=1,j=0;k<n;k<<=1) {
			j <<= 1;
			if (i & k) j |= 1;
		}
                 
		if (i < j)
		  mpfi_c_swap(x[i],x[j]);
		else if (i > j)
			mpfi_c_swap(x[n-1-i],x[n-1-j]);
		++i, j |= l;
		mpfi_c_swap(x[i],x[j]);
	}

	for (k=1,l=n/2;k<n;k<<=1,l>>=1)
		for (p=x;p<xend;p+=k)
			for (j=0;j<n/2;j+=l,p++) {
				mpfi_c_mul(ctemp,p[k],w[j]);
				mpfi_c_sub(p[k],p[0],ctemp);
				mpfi_c_add(p[0],p[0],ctemp);
			}
}

/* perform an in place inverse FFT */
void ifft(mpfi_c_t *x,unsigned int n, mpfi_c_t *w) {
	int i,l=n>>1;

	fft(x,n,w);
	mpfi_c_div_ui(x[0],x[0],n);
	mpfi_c_div_ui(x[l],x[l],n);
	for (i=1;i<l;i++) {
		mpfi_c_div_ui(x[i],x[i],n);
		mpfi_c_div_ui(x[n-i],x[n-i],n);
		mpfi_c_swap(x[i],x[n-i]);
	}
}

/* perform a non-normalised inverse FFT */
void nifft(mpfi_c_t *x,unsigned int n, mpfi_c_t *w) {
	int i,l=n>>1;

	fft(x,n,w);
	for (i=1;i<l;i++)
	  mpfi_c_swap(x[i],x[n-i]);
}


/* circular convolution; f and g must be distinct */
void convolve(mpfi_c_t *result,mpfi_c_t *f,mpfi_c_t *g,unsigned int n, mpfi_c_t *w) {
	uint64_t i;
	fft(f,n,w);
        
        fft(g,n,w);
	for (i=0;i<n;i++)
	  mpfi_c_mul(result[i],f[i],g[i]);
	ifft(result,n,w);
}

// Do iDFT on a Hermitian vector of length 2N, i.e. one that produces real values
// N is the reduced length and the ws relate to that N
// ws[n]=e(n/N) for n=0..N/2-1
// x is of length 2N with 0..N input values and (N+1..2N-1) anything
// omega=exp(2 pi I/2/N)
// on exit x[n]->re are set, x[n]->im are garbage.
//
void hermidft(mpfi_c_t *x, int N, mpfi_c_t *ws, mpfi_c_ptr omega)
{
  int n,m;
  mpfi_c_t *res=x+N; // use back "half" of vector for iDFT

  mpfi_sub(res[0]->im,x[0]->re,x[N]->re);
  mpfi_add(res[0]->re,x[0]->re,x[N]->re); // x[N] now destroyed
  for(n=1;n<N;n++)
    {
      // compute e(n/(2N))*i=e((n+N/2)/(2N))
      mpfi_c_conj(ctemp,x[N-n]);
      mpfi_c_set(res[n],ctemp);
      mpfi_c_add(ctemp,ctemp,x[n]);
      mpfi_c_sub(res[n],x[n],res[n]);
      mpfi_c_muli(res[n]);
      m=(n/2)%N;
      if(m>=N/2)
	{
	  mpfi_c_mul(res[n],res[n],ws[m-N/2]);
	  mpfi_c_neg(res[n],res[n]);
	}
      else
	mpfi_c_mul(res[n],res[n],ws[m]);
      if(n&1)
	mpfi_c_mul(res[n],res[n],omega);
      mpfi_c_add(res[n],res[n],ctemp);
    }
  fft(res,N,ws);
  for(n=0,m=0;n<N+N;n+=2,m++)
    {
      mpfi_set(x[n]->re,res[m]->re);
      mpfi_set(x[n+1]->re,res[m]->im);
    }
}

mpfi_c_t *simple_ws[MAX_SIMPLE_DFT];
mpfi_c_t dft_spare[MAX_SIMPLE_DFT+1];

void dft_init(uint64_t prec)
{
  uint64_t n;

  mpfi_c_setup(prec);

  for(n=0;n<=MAX_SIMPLE_DFT;n++)
    mpfi_c_init(dft_spare[n]);
  for(n=0;n<MAX_FFT_LEN;n++)
    {
      mpfi_c_init(X[n]);
      mpfi_c_init(A[n]);
    }
  mpfi_c_init(ctemp);
}

// return smallest power of 2 > 2n-1 
uint64_t calc_conv_size(uint64_t n)
{
  uint64_t res=1;
  while(res<(n<<1)-1)
    res<<=1;
  return(res);
}

// use Bluestein's algorithm to do a length len DFT on x
// where each datum is stride apart
void bluesteinfft(mpfi_c_t *x, uint64_t len, uint64_t stride)
{
  //printf("In Bluestein FFT with len=%lu\n",len);
  uint64_t i,i2,l,s;
  int64_t l2;
  mpfi_c_t *w;
  if(B[len].conv_size==0) // not seen this size before as a bluestein
    {
      //printf("New length for Bluestein.\n");
      l=calc_conv_size(len);
      l2=log_2(l);
      if(l2<0)
	fatal_error("Bad vector length in bluesteinfft.");
      //printf("Full length =%lu log2=%lu\n",l,l2);
      B[len].conv_size=l;
      B[len].b=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*l);
      B[len].b_conj=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*len);
      if((!B[len].b)||(!B[len].b_conj))
	fatal_error("Error allocating memory for Bluestein.");

      mpfi_c_init(B[len].b[0]);
      mpfi_c_set_ui(B[len].b[0],1,0);
      mpfi_div_ui(dft_spare[0]->re,mpfi_pi,len); // pi/n
      // this should be OK till 4*len overflows 64 bits!!
      for(i=1,i2=1;i<len;i++,i2=(i2+i+i-1)%(len<<1))
	{
	  mpfi_mul_ui(dft_spare[1]->re,dft_spare[0]->re,i2); // pi k^2 /n
	  mpfi_c_init(B[len].b[i]);
	  mpfi_cos(B[len].b[i]->re,dft_spare[1]->re);
	  mpfi_sin(B[len].b[i]->im,dft_spare[1]->re);
	  mpfi_c_init(B[len].b[l-i]);
	  mpfi_c_set(B[len].b[l-i],B[len].b[i]);
	}
      for(;i<=l-len;i++)
	{
	  mpfi_c_init(B[len].b[i]);
	  mpfi_c_set_ui(B[len].b[i],0,0);
	}
      for(i=0;i<len;i++)
	{
	  mpfi_c_init(B[len].b_conj[i]);
	  mpfi_c_conj(B[len].b_conj[i],B[len].b[i]);
	}

      if(!ws[l2]) // we've not seen this conv_size as a straight FFT
	          // so construct the ws
	{
	  ws[l2]=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*l/2);
	  if(!ws[l2])
	    fatal_error("Failed to allocate memory for ws.");
	  w=ws[l2];
	  mpfi_c_init(w[0]);
	  mpfi_c_set_ui(w[0],1,0);
	  mpfi_div_d(dft_spare[0]->re,mpfi_2_pi,-((double) l)); // -2Pi/N
	  for(i=1;i<l/2;i++)
	    {
	      mpfi_mul_ui(dft_spare[0]->im,dft_spare[0]->re,i);
	      mpfi_c_init(w[i]);
	      mpfi_cos(w[i]->re,dft_spare[0]->im);
	      mpfi_sin(w[i]->im,dft_spare[0]->im);
	    }
	}
      else // we have seen this length as an fft, so just point to ws
	w=ws[l2];
      // do the FFT of b once only
      fft(B[len].b,l,w);
    }
  else // seen this length as a bluestein
    {
      l=B[len].conv_size;
      l2=log_2(l);
      w=ws[l2];
    }
  
  mpfi_c_set(A[0],x[0]);
  for(i=1,s=stride;i<len;i++,s+=stride)
    mpfi_c_mul(A[i],x[s],B[len].b_conj[i]);
  for(;i<l;i++)
    mpfi_c_set_ui(A[i],0,0);

  // now do the convolution, but B[len].b has been dft'd already
  fft(A,l,w);
  for(i=0;i<l;i++)
    mpfi_c_mul(A[i],A[i],B[len].b[i]);
  ifft(A,l,w);
  mpfi_c_set(x[0],A[0]);
  for(i=1,s=stride;i<len;i++,s+=stride)
    mpfi_c_mul(x[s],A[i],B[len].b_conj[i]);
}


// do a simple dft by O(n^2) method
// input in x[0],x[stride].., output to x 
void simple_dft(mpfi_c_t *x, unsigned long int N, unsigned long int stride)
{
  long unsigned int n,m,ptr;
  if(!simple_ws[N-1]) // never seen this length before
    {
      simple_ws[N-1]=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*N);
      if(!simple_ws[N-1])
	fatal_error("Failed to allocate memory for simple_ws.");
      mpfi_c_init(simple_ws[N-1][0]);
      mpfi_c_one(simple_ws[N-1][0]);
      mpfi_div_ui(dft_spare[0]->im,mpfi_2_pi,N); // 2 Pi/N
      mpfi_neg(dft_spare[0]->re,dft_spare[0]->im); // -2 Pi/N
      for(n=1;n<N;n++)
	{
	  //printf("Doing %lu\n",n);
	  mpfi_mul_ui(dft_spare[0]->im,dft_spare[0]->re,n); // -2pi n/N
	  mpfi_c_init(simple_ws[N-1][n]);
	  mpfi_cos(simple_ws[N-1][n]->re,dft_spare[0]->im);
	  mpfi_sin(simple_ws[N-1][n]->im,dft_spare[0]->im);
	}
    }
  /*
  for(n=0;n<N;n++)
    {
      printf("ws[%lu]=",n);
      mpfi_c_print(simple_ws[N-1][n]);
    }
  */

  mpfi_c_set(dft_spare[0],x[0]);
  for(m=1;m<N;m++)
    mpfi_c_add(dft_spare[0],dft_spare[0],x[m*stride]);
  for(n=1;n<N;n++)
    {
      mpfi_c_set(dft_spare[n],x[0]); // x[0]*e(0)
      for(m=1;m<N;m++)
	{
	  ptr=(m*n)%N;
	  mpfi_c_mul(dft_spare[MAX_SIMPLE_DFT],x[m*stride],simple_ws[N-1][ptr]); // x[m]+e(-mn/N)
	  //printf("term for n=%lu m=%lu is",n,m);mpfi_c_print(dft_spare[MAX_SIMPLE_DFT]);
	  mpfi_c_add(dft_spare[n],dft_spare[n],dft_spare[MAX_SIMPLE_DFT]);
	}
    }
  for(n=0;n<N;n++)
    mpfi_c_set(x[n*stride],dft_spare[n]);
}

void dft(mpfi_c_t *x, unsigned long int len, unsigned long int stride)
{
  //printf("In dft with len=%lu and stride=%lu\n",len,stride);
  if(len==2) // easy
    {
      mpfi_c_add(dft_spare[0],x[0],x[stride]);
      mpfi_c_sub(dft_spare[1],x[0],x[stride]);
      mpfi_c_set(x[0],dft_spare[0]);
      mpfi_c_set(x[stride],dft_spare[1]);
      return;
    }
  if(len==4) // almost as easy
    {
      mpfi_c_add(dft_spare[0],x[0],x[stride<<1]);
      mpfi_c_add(dft_spare[1],x[stride],x[3*stride]);
      mpfi_c_sub(dft_spare[2],x[0],x[stride<<1]);
      mpfi_c_sub(dft_spare[3],x[stride],x[3*stride]);
      mpfi_c_add(x[0],dft_spare[0],dft_spare[1]);
      mpfi_c_sub(x[stride<<1],dft_spare[0],dft_spare[1]);
      mpfi_c_muli(dft_spare[3]);
      mpfi_c_sub(x[stride],dft_spare[2],dft_spare[3]);
      mpfi_c_add(x[3*stride],dft_spare[2],dft_spare[3]);
      return;
    }
  if(len<=MAX_SIMPLE_DFT) // use the O(n^2) algorithm
    {
      simple_dft(x,len,stride);
      return;
    }
  
  if(power2(len)) // use the Gauss FFT
    {
      int64_t l2=log_2(len);
      uint64_t i,s;
      if(l2>LOG2_MAX_FFT_LEN)
	fatal_error("FFT too long.");
      mpfi_c_t *w;
      if(!ws[l2])
	{
	  ws[l2]=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*len/2);
	  if(!ws[l2])
	    fatal_error("Failed to allocate memory for ws.");
	  w=ws[l2];
	  mpfi_c_init(w[0]);
	  mpfi_c_set_ui(w[0],1,0);
	  mpfi_div_d(dft_spare[0]->re,mpfi_2_pi,-((double) len)); // -2Pi/N
	  for(i=1;i<len/2;i++)
	    {
	      mpfi_mul_ui(dft_spare[0]->im,dft_spare[0]->re,i);
	      mpfi_c_init(w[i]);
	      mpfi_cos(w[i]->re,dft_spare[0]->im);
	      mpfi_sin(w[i]->im,dft_spare[0]->im);
	    }
	}
      else
	w=ws[l2];

      if(stride==1) // data is contiguous, so just do it
	fft(x,len,w);
      else // copy it to make it contiguous
	{
	  for(i=0,s=0;i<len;i++,s+=stride)
	    mpfi_c_swap(X[i],x[s]); // is this faster than mpfi-c_set ??
	  fft(X,len,w);
	  for(i=0,s=0;i<len;i++,s+=stride)
	    mpfi_c_swap(x[s],X[i]);
	}
      return;
    }
  // going to have to use the full Bluestein algorithm
  if(len>MAX_FFT_LEN/2)
    fatal_error("FFT too big for Bluestein.");
  bluesteinfft(x,len,stride);
}

//
// do an ndft
// total length of x is N made up of num_dims dimensions whose lengths are
// specified in dims
//
void ndft(mpfi_c_t *x, unsigned long int N, unsigned long int num_dims, unsigned long int *dims)
{
  unsigned long int i,k,l,stride=N,depth;

  for(i=0;i<num_dims;i++)
    {
      depth=N/stride;
      stride/=dims[i];
      for(k=0;k<depth*dims[i]*stride;k+=dims[i]*stride)
	for(l=0;l<stride;l++) {
	    dft(x+k+l,dims[i],stride);
        }
    }
}
