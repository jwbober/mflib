#include <iostream>

#include "cuspforms_weight1_modp.h"

#include "classnumbers.h"

using namespace std;

int main(int argc, char ** argv) {
    init_classnumbers();
    load_factor_table();
    if(argc < 5) {
        cout << "Usage: " << argv[0] << " level chi p ncoeffs" << endl;
        return 0;
    }
    int level = atoi(argv[1]);
    int chinumber = atoi(argv[2]);
    int p = atoi(argv[3]);
    int ncoeffs = atoi(argv[4]);

    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chinumber);

    cuspforms_weight1_modp * S = get_cuspforms_weight1_modp(chi, p);
    p = S->p;


    nmod_mat_t basis;

    S->basis(basis, ncoeffs);

    int d = S->dimension();
    for(int k = 0; k < d; k++) {
        cout << nmod_mat_entry(basis, k, 0);
        for(int j = 1; j < ncoeffs; j++) {
            cout << " " << nmod_mat_entry(basis, k, j);
        }
        cout << endl;
    }

    return 0;
}
