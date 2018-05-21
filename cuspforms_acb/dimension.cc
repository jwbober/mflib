#include <iostream>

#include "cuspforms_acb.h"

using std::cout;
using std::cerr;
using std::endl;

long cuspforms_acb::new_dimension() {
    if(_new_dimension != -1) return _new_dimension;
    if(traces_computed < 2) compute_traces(2);
    long a = arf_abs_bound_lt_2exp_si(arb_midref(acb_realref(traces[1])));
    if(a > 62) {
        cerr << "Something has gone wrong and we getting the wrong"
                "result for the dimension. Or (less likely) the dimension"
                "is too big to store as a long." << endl;
        exit(0);
    }
    fmpz_t z;
    fmpz_init(z);
    arb_get_unique_fmpz(z, acb_realref(traces[1]));
    long d = fmpz_get_si(z);
    fmpz_clear(z);
    _new_dimension = d;
    return d;
}

