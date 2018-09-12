#ifndef _ARB_EXTRAS_H_
#define _ARB_EXTRAS_H_

#include "acb_poly.h"
#include "arb_mat.h"

#ifdef __cplusplus

#include <ostream>

std::ostream& operator << (std::ostream& out, mag_t in);
std::ostream& operator << (std::ostream& out, arf_t in);
std::ostream& operator << (std::ostream& out, arb_t in);
std::ostream& operator << (std::ostream& out, acb_t in);

extern "C" {
#endif

void twobytwo_diag(arb_t u1, arb_t u2, const arb_t a, const arb_t b, const arb_t d, slong prec);
void arb_mat_cholesky(arb_mat_t out, const arb_mat_t in, slong prec);
int arb_mat_jacobi(arb_mat_t e, arb_mat_t E, const arb_mat_t A, slong prec);
int arb_mat_generalized_eigenproblem_symmetric_positive_definite(arb_mat_t D, arb_mat_t P, const arb_mat_t A, const arb_mat_t B, slong prec);
void arb_mat_print_sage_float(const arb_mat_t A);

void fprintacb(FILE *fp, const acb_t x);
void fprintarb(FILE *fp, const arb_t x);
void fprintarf(FILE *fp, const arf_t x);

void acb_reasonable_sqrt(acb_t out, const acb_t in, slong prec);

int arb_poly_gcd(arb_poly_t res,const arb_poly_t a,const arb_poly_t b,slong d, slong prec);
int arb_poly_sign_change(arb_poly_t f,arb_srcptr x,slong prec); // sgn(f(a)*f(b)), where x=[a,b]
int arb_get_unique_fmpz_modular(fmpz_t out, arb_t real_approx, fmpz_t modular_approx, fmpz_t mod);
int arb_poly_get_unique_fmpz_poly_modular(fmpz_poly_t out, arb_poly_t real_approx, fmpz_poly_t modular_approx, fmpz_t mod);
int acb_poly_get_unique_fmpz_poly_modular(fmpz_poly_t out, acb_poly_t real_approx, fmpz_poly_t modular_approx, fmpz_t mod);
void arb_get_ubound_fmpz(fmpz_t out, arb_t in, slong prec);

#ifdef __cplusplus
}
#endif
#endif
