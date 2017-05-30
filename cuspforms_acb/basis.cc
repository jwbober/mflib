#include <vector>

#include "cuspforms_acb.h"

using std::vector;

void cuspforms_acb::newspace_basis(acb_mat_t B, int ncoeffs) {
    newspace_basis_data();
    int d = new_dimension();

    compute_traces(ncoeffs * basis_cols[d-1]);
    acb_mat_init(B, ncoeffs, d);
    for(int k = 0; k < d; k++) {
        int n = basis_cols[k];
        for(int m = 0; m < ncoeffs; m++) {
            trace_TnTm(acb_mat_entry(B, m, k), n, m);
        }
    }
}
const vector<int>& cuspforms_acb::newspace_basis_data() {
    if(basis_cols.size() > 0) return basis_cols;
    int d = new_dimension();
    if(d == 0) return basis_cols;
    basis_cols = modp_space->newspace_basis_data();
    return basis_cols;
}

