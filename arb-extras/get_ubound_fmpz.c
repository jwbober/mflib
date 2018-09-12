#include "arb.h"

void arb_get_ubound_fmpz(fmpz_t out, arb_t in, slong prec) {
    arf_t t;
    arf_init(t);
    arb_get_ubound_arf(t, in, prec);
    arf_get_fmpz(out, t, ARF_RND_UP);
    arf_clear(t);
}
