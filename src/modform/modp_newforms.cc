#include "modform_modp.h"
#include "classnumbers.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

int main(int argc, char ** argv) {
    int level;
    int chi_number;
    int weight;
    int ncoeffs;
    long p;

    init_classnumbers();

    level = atoi(argv[1]);
    weight = atoi(argv[2]);
    chi_number = atoi(argv[3]);
    ncoeffs = atoi(argv[4]);
    p = atol(argv[5]);

    cout << level << endl
         << weight << endl
         << chi_number << endl
         << ncoeffs << endl
         << p << endl;

    int verbose = 0;
    if(argc > 6) verbose = atoi(argv[6]);

    DirichletGroup G(level);
    if(GCD(level, chi_number) != 1) return 0;
    DirichletCharacter chi = G.character(chi_number);
    if(chi.is_even() && weight % 2 == 1) return 0;
    if(!chi.is_even() && weight % 2 == 0) return 0;

    cuspforms_modp * S = get_cuspforms_modp(chi, weight, p, verbose);
    nmod_mat_t newforms;
    S->newforms(newforms, ncoeffs);

    return 0;
}
