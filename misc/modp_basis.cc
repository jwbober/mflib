#include <iostream>

#include "cuspforms_modp.h"
#include "classnumbers.h"

using namespace std;

int main(int argc, char ** argv) {
    init_classnumbers();
    int level = atoi(argv[1]);
    int weight = atoi(argv[2]);
    int chi_number = atoi(argv[3]);
    long p = atoi(argv[4]);
    int ncoeffs = atoi(argv[5]);

    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chi_number);

    cuspforms_modp * cuspforms = get_cuspforms_modp(chi, weight, p);
    nmod_mat_t basis;
    cuspforms->basis(basis, ncoeffs);
    nmod_mat_print_pretty(basis);
    return 0;
}
