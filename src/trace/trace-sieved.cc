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
    int level = atoi(argv[1]);
    int chi_number = atoi(argv[2]);
    long p = atol(argv[3]);
    int start = atoi(argv[4]);
    int end = atoi(argv[5]);

    init_classnumbers();

    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chi_number);

    if(!chi.is_even()) {
        return 0;
    }

    long primitive_index;
    int q = chi.conductor(&primitive_index);
    long ** chi_values = new long*[level + 1];
    complex<double> ** chi_zvalues = new complex<double>*[level + 1];

    vector<int> sublevels;                     //
    for(int M : divisors(level)) {             // We fill sublevels with all of the M such that
        if(M % q == 0) {                       // q | M | level.
            if(M > 1) sublevels.push_back(M);  //
            DirichletGroup G(M);
            DirichletCharacter psi = G.character(G.index_from_primitive_character(q, primitive_index));
            chi_values[M] = new long[M];
            chi_zvalues[M] = new complex<double>[M];
            psi.values_mod_p(p, chi_values[M]);
            for(int k = 0; k < M; k++) {
                chi_zvalues[M][k] = psi.value(k);
            }
        }
    }

    vector<long> * traces = new vector<long>[level + 1];
    vector<complex<double>> * ztraces = new vector<complex<double>>[level + 1];

    for(int M : sublevels) {                                                                            //
        traces[M] = vector<long>(end + 1);                                                         // We compute the (unsieved) traces
        ztraces[M] = vector<complex<double>>(end + 1);
        trace_Tn_modp_unsieved_weight2(traces[M].data(), 0, end + 1, M, p, chi_values[M], chi);   // for each sublevel.
        trace_Tn_unsieved_weight2(ztraces[M].data(), 0, end + 1, M, chi_zvalues[M], chi);
    }                                                                                                   //

    sieve_trace_Tn_modp_on_weight2_for_newspaces(traces, 0, end + 1, level, p, chi_values, chi);  // Then we sieve.
    sieve_trace_Tn_on_weight2_for_newspaces(ztraces, 0, end + 1, level, chi_zvalues, chi);  // Then we sieve.

    cout << p << endl;
    for(int M : sublevels) {
        cout << M << "\t";
        for(int k = start; k <= end; k++) {
            cout << traces[M][k];
            if(k < end) cout << "\t";
        }
        cout << endl;
    }

    //for(int M : sublevels) {
    //    cout << M << "\t";
    //    for(int k = start; k <= end; k++) {
    //        cout << ztraces[M][k];
    //        if(k < end) cout << "\t";
    //    }
    //    cout << endl;
    //}
    return 0;
}
