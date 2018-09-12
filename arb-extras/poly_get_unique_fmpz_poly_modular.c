#include "arb-extras.h"

int arb_poly_get_unique_fmpz_poly_modular(fmpz_poly_t out, arb_poly_t real_approx, fmpz_poly_t modular_approx, fmpz_t mod) {
    fmpz_t t;
    fmpz_init(t);
    slong degree = arb_poly_degree(real_approx);
    for(slong k = 0; k <= degree; k++) {
        int result = arb_get_unique_fmpz_modular(
                t,
                arb_poly_get_coeff_ptr(real_approx, k),
                fmpz_poly_get_coeff_ptr(modular_approx, k),
                mod);
        if(result == 0) {
            fmpz_clear(t);
            return 0;
        }
        fmpz_poly_set_coeff_fmpz(out, k, t);
    }
    fmpz_clear(t);
    return 1;
}

int acb_poly_get_unique_fmpz_poly_modular(fmpz_poly_t out, acb_poly_t real_approx, fmpz_poly_t modular_approx, fmpz_t mod) {
    fmpz_t t;
    fmpz_init(t);
    slong degree = acb_poly_degree(real_approx);
    for(slong k = 0; k <= degree; k++) {
        if(!arb_contains_zero(acb_imagref(acb_poly_get_coeff_ptr(real_approx, k)))) {
            fmpz_clear(t);
            return 0;
        }
        int result = arb_get_unique_fmpz_modular(
                t,
                acb_realref(acb_poly_get_coeff_ptr(real_approx, k)),
                fmpz_poly_get_coeff_ptr(modular_approx, k),
                mod);
        if(result == 0) {
            fmpz_clear(t);
            return 0;
        }
        fmpz_poly_set_coeff_fmpz(out, k, t);
    }
    fmpz_clear(t);
    return 1;
}
