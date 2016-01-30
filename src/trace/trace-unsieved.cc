#include <iostream>
#include <cstdlib>
#include <complex>

#include "flint/nmod_mat.h"

#include "characters.h"
#include "S2dimensions.h"
#include "trace-formula.h"
#include "classnumbers.h"

using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::complex;

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
    if(argc < 6) {
        const char * usage =
            "./trace-unsieved level chi p start end\n"
            "\n"
            "Print some traces for chi mod level.\n";
        cout << usage;
        return 0;
    }
    int level = atoi(argv[1]);
    int chi_number = atoi(argv[2]);
    long p0 = atol(argv[3]);
    int start = atoi(argv[4]);
    int end = atoi(argv[5]);

    init_classnumbers();

    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chi_number);

    long * traces = new long[end - start];
    complex<double> * ztraces = new complex<double>[end - start];
    long * chi_values = new long[level];
    complex<double> * chi_zvalues = new complex<double>[level];
    chi.values_mod_p(p0, chi_values);
    for(int k = 0; k < level; k++) {
        chi_zvalues[k] = chi.value(k);
    }
    //for(int k = 0; k < level; k++) {
    //    cout << k << " " << chi_values[k] << " " << chi.value(k) << endl;
    //}

    trace_Tn_modp_unsieved_weight2(traces, start, end, level, p0, chi_values, chi);
    trace_Tn_unsieved_weight2(ztraces, start, end, level, chi_zvalues, chi);
    cout << p0 << endl;
    for(int m = start; m < end; m++) {
        int k = m - start;
        int ok = 0;
        if(imag(ztraces[k]) == 0.0 && (int)real(ztraces[k]) == real(ztraces[k])) {
            int z = (int)real(ztraces[k]);
            z = z % p0;
            if(z < 0) z += p0;
            if(z == traces[k]) ok = 1;
            else ok = -1;
        }
        cout << m << " " << traces[k] << " " << ztraces[k] << " ";
        if(ok == 1) cout << "ok" << endl;
        else if(ok == -1) cout << "ohno" << endl;
        else cout << "?" << endl;
    }
    return 0;
}
