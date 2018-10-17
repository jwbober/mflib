#include "cuspforms_weight1_modp.h"
#include "slint.h"

int cuspforms_weight1_modp::new_dimension() {
    int d = dimension();
    for(cuspforms_weight1_modp * S: subspaces) {
        int d2 = S->new_dimension();
        int M1 = S->level;
        d = d - ndivisors(level/M1)*d2;
    }
    return d;
}
