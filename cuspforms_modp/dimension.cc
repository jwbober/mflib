#include "cuspforms_modp.h"

int cuspforms_modp::new_dimension() {
    return trace(1);
}

int cuspforms_modp::dimension() {
    compute_traces(2);
    return _dimension;
}
