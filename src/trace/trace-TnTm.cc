#include <iostream>
#include <cstdlib>
#include <string>

#include "flint/nmod_mat.h"

#include "characters.h"
#include "S2dimensions.h"
#include "trace-formula.h"

using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

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
        string usage =  "./trace-TnTm level chi p N M\n"
                        "\n"
                        "For n <= N, m <= M, and every sublevel of S_2(level, chi), print out the\n"
                        "trace of TnTm.\n";
        cout << usage;
        return 0;
    }
    int level = atoi(argv[1]);
    int chi_number = atoi(argv[2]);
    int p = atoi(argv[3]);
    int NN = atoi(argv[4]);
    int MM = atoi(argv[5]);

    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chi_number);

    if(!chi.is_even()) {
        return 0;
    }

    long primitive_index;
    int q = chi.conductor(&primitive_index);
    int ** chi_values = new int*[level + 1];

    vector<int> sublevels;                     // 
    for(int M : divisors(level)) {             // We fill sublevels with all of the M such that 
        if(M % q == 0) {                       // q | M | level.
            if(M > 1) sublevels.push_back(M);  //
            DirichletGroup G(M);
            DirichletCharacter psi = G.character(G.index_from_primitive_character(q, primitive_index));
            chi_values[M] = new int[M];
            psi.values_mod_p(p, chi_values[M]);
        }
    }
    
    vector<int> * traces = new vector<int>[level + 1];

    for(int M : sublevels) {                                                                      //
        traces[M] = vector<int>(NN*MM + 1);                                                       // We compute the (unsieved) traces
        trace_Tn_modp_unsieved_weight2(traces[M].data(), 0, NN*MM + 1, M, p, chi_values[M], chi); // for each sublevel.
    }                                                                                             //


    sieve_trace_Tn_modp_on_weight2_for_newspaces(traces, 0, NN*MM + 1, level, p, chi_values, chi);// Then we sieve.
    
    cout << p << endl;
    for(int M : sublevels) {
        cout << M << endl;
        for(int n = 0; n <= NN; n++) {
            for(int m = 0; m <= MM; m++) {
            cout << trace_TmTn_mod_p(traces[M].data(), chi_values[M], 2, M, n, m, p);
            if(m < MM) cout << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    return 0;
}
