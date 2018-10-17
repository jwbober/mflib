#include <iostream>

#include "cuspforms_weight1_modp.h"

#include "classnumbers.h"

using namespace std;

int main(int argc, char ** argv) {
    int level = atoi(argv[1]);
    int chinumber = atoi(argv[2]);
    int p = atoi(argv[3]);
    int n1 = atoi(argv[4]);
    int n2 = atoi(argv[5]);

    init_classnumbers();
    load_factor_table();

    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chinumber);

    cuspforms_weight1_modp * S = get_cuspforms_weight1_modp(chi, p);
    p = S->p;

    nmod_mat_t Tp;
    for(int n = n1; n <= n2; n++) {
        S->hecke_matrix(Tp, n);
        nmod_mat_print_pretty(Tp);
        nmod_mat_clear(Tp);
    }
}
