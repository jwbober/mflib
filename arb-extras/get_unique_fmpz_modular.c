#include "arb.h"

int arb_get_unique_fmpz_modular(fmpz_t out, arb_t real_approx, fmpz_t modular_approx, fmpz_t mod) {
    // compute and return the unique integer in the interval represented by real_approx
    // and congruent to modular_approx modulo mod, if such an integer exists.
    //
    // Returns 1 on success, 0 if there is no such integer.

    fmpz_t a, b, e;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(e);

    arb_get_interval_fmpz_2exp(a, b, e, real_approx);

    // we really should be more careful with error checking here.
    slong exponent = fmpz_get_si(e);
    if(exponent > 0)
        fmpz_mul_2exp(a, a, exponent);
    else
        fmpz_fdiv_q_2exp(a, a, -exponent);
    fmpz_mod(b, a, mod);
    fmpz_sub(a, a, b);
    fmpz_add(a, a, modular_approx);

    // We are probably in the interval now, if there is such an integer
    // in the interval. If I were to think properly about what I'm doing,
    // then I could just check a + mod probably to make sure that it isn't
    // also in the interval. But there is rounding and stuff to worry about...

    int nfound = 0;
    fmpz_sub(a, a, mod);
    fmpz_sub(a, a, mod);

    if(arb_contains_fmpz(real_approx, a)) {
        nfound += 1;
        fmpz_set(out, a);
    }
    for(int k = 0; k < 4; k++) {
        fmpz_add(a, a, mod);
        if(arb_contains_fmpz(real_approx, a)) {
            nfound += 1;
            fmpz_set(out, a);
        }
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(e);

    if(nfound != 1) return 0;

    return 1;
}
