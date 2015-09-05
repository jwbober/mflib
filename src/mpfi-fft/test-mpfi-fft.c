//
// version of mpfi_fft.h includes all necessary mpfi_c stuff from mpfi_c.h
// 
//

#include "mpfi_fft.h"
#include "inttypes.h" // uint64_t etc.


int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("Usage:- %s <PREC>.\n",argv[0]);
      exit(0);
    }

  uint64_t prec; // how many bits do we want MPFI to use?
  prec=atoi(argv[1]);
  if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
   {
      printf("Usage:- %s <PREC>.\n",argv[0]);
      exit(0);
    }
  dft_init(prec);

#define Q (5*37*41)
#define PHI_Q (4*36*40)

  mpfi_c_t vec[PHI_Q];
  uint64_t dims[3]={4,36,40};

  uint64_t i;
  for(i=0;i<PHI_Q;i++)
    {
      mpfi_c_init(vec[i]);
      mpfi_c_set_ui(vec[i],rand(),0);
    }

  //
  // this will do 36*40 length 4 dft's of i, i+36*40, i+2*36*40 and i+3*36+40
  // then 4*40 length 36 dfts
  // then 4*36 length 40 dfts
  //

  ndft(vec,PHI_Q,3,dims);

  mpfi_c_print_str("res[100]: ",vec[100]);

  return(0);
}
