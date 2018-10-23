#include <iostream>

#include "cuspforms_weight1_modp.h"
#include "classnumbers.h"

using namespace std;

int main(int argc, char ** argv) {
    init_classnumbers();
    load_factor_table();

    int level = atoi(argv[1]);
    int chinumber = atoi(argv[2]);
    int ncoeffs = atoi(argv[3]);

    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chinumber);

    cuspforms_weight1_modp * S = get_cuspforms_weight1_modp(chi, 50000);
    int p = S->p;

    nmod_mat_t newforms;

    int newform_count = S->newspace_basis(newforms, ncoeffs);
    int d = S->new_dimension();

    cout << p << endl;
    for(int k = 0; k < d; k++) {
        if(k < newform_count) cout << "* ";
        cout << nmod_mat_entry(newforms, d - k - 1, 0);
        for(int j = 1; j < ncoeffs; j++) {
            cout << " " << nmod_mat_entry(newforms, k, j);
        }
        cout << endl;
    }

    return 0;
}
