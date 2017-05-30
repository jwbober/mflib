#include "cuspforms_acb.h"

void cuspforms_acb::evalpoly(fmpz_t out, fmpz_t t, fmpz_t n) {
    // translated from Ralph's code
    int k = weight - 2;
    if(k == 0) {fmpz_set_ui(out, 1); return;}
    if(k == 1) {fmpz_set(out, t); return;}

    fmpz_t val[2];
    fmpz_init(val[0]);
    fmpz_init(val[1]);

    fmpz_set_ui(val[0], 1);
    fmpz_set(val[1], t);

    fmpz_t x, y;
    fmpz_init(x);
    fmpz_init(y);
    for(int i = 2; i <= k; i++) {
        fmpz_mul(x, t, val[(i + 1) % 2]);
        fmpz_mul(y, n, val[i % 2]);
        fmpz_sub(val[i % 2], x, y);
    }
    fmpz_set(out, val[k % 2]);

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(val[0]);
    fmpz_clear(val[1]);
}
