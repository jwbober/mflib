#include <iostream>
#include <cstdlib>

#include "flint/nmod_mat.h"

#include "characters.h"
#include "S2dimensions.h"
#include "trace-formula.h"

using std::vector;
using std::cerr;
using std::cout;
using std::endl;

void print_nmod_mat_t(nmod_mat_t M) {
    int nrows = nmod_mat_nrows(M);
    int ncols = nmod_mat_ncols(M);
    cout << "matrix(Integers(" << M->mod.n << "),[" << endl;
    for(int i = 0; i < nrows; i++) {
        cout << " [";
        for(int j = 0; j < ncols; j++) {
            cout << nmod_mat_entry(M, i, j);
            if(j + 1 < ncols) cout << ", ";
        }
        cout << "]";
        if(i + 1 < nrows) cout << ",";
        cout << endl;
    }
    cout << "])" << endl;
}

int main(int argc, char ** argv) {
    int level = atoi(argv[1]);
    int chi_number = atoi(argv[2]);
    int p0 = atoi(argv[3]);
    int start = atoi(argv[4]);
    int end = atoi(argv[5]);

    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chi_number);

    int * traces = new int[end - start];
    int * chi_values = new int[level];
    chi.values_mod_p(p0, chi_values);

    trace_Tn_modp_unsieved_weight2(traces, start, end, level, p0, chi_values, chi);
    cout << p0 << endl;
    for(int k = start; k < end; k++) {
        cout << k << " " << traces[k] << endl;
    }
    return 0;
}
