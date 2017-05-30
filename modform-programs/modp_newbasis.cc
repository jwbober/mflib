#include "cuspforms_modp.h"
#include "classnumbers.h"

#include <iostream>

using namespace std;

int main(int argc, char ** argv) {
    int level;
    int weight;
    int chi_number;
    int p;
    int ncoeffs;

    init_classnumbers();

    level = atoi(argv[1]);
    weight = atoi(argv[2]);
    chi_number = atoi(argv[3]);
    p = atoi(argv[4]);
    ncoeffs = atoi(argv[5]);

    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chi_number);

    cuspforms_modp * S = get_cuspforms_modp(chi, weight, p);
    p = S->p;
    vector<int> basis_rows = S->newspace_basis_data();
    //for(int k : basis_rows) cout << k << " ";
    //cout << endl;
    nmod_mat_t basis;
    S->newspace_basis(basis, ncoeffs);
    cout << p << endl;
    for(int k : basis_rows) cout << k << " ";
    cout << endl;
    for(int k = 0; k < S->new_dimension(); k++) {
        for(int j = 0; j < ncoeffs; j++) {
            long a = nmod_mat_entry(basis, k, j);
            if(a > p/2) a-= p;
            cout << a << " ";
        }
        cout << endl;
    }
    cout << nmod_mat_rank(basis) << endl;
    //for(int n = 0; n < 15; n++) {
    //    for(int m = 0; m < 15; m++) {
    //        cout << S->trace_TnTm(n, m) << " ";
    //    }
    //    cout << endl;
    //}

    return 0;
}
