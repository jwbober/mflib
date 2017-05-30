#include "acb.h"

void acb_reasonable_sqrt(acb_t out, const acb_t in, slong prec) {
    if(arb_is_negative(acb_realref(in)) && arb_contains_zero(acb_imagref(in))) {
        acb_neg(out, in);
        acb_sqrt(out, out, prec);
        acb_mul_onei(out, out);
        if(arf_sgn(arb_midref(acb_imagref(in))) < 0) {
            acb_neg(out, out);
        }
    }
    else {
        acb_sqrt(out, in, prec);
    }
}
