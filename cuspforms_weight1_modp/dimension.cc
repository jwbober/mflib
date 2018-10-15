#include "cuspforms_weight1_modp.h"

int cuspforms_weight1_modp::dimension() {
    if(_dimension == -1) {
        compute_basis_data();
    }
    return _dimension;
}
