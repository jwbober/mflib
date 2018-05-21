#include "cuspforms_acb.h"
#include "arb-extras.h"
#include "classnumbers.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>

#include "mag.h"

#include "mfformat.h"

using namespace std;

long arb_abs_prec_approx(arb_t x) {
    if(arb_is_exact(x)) return ARF_PREC_EXACT;
    double mag_bits = mag_get_d_log2_approx(arb_radref(x));
    return floor(mag_bits);
}

int main(int argc, char ** argv) {
    srand(time(NULL));
    int level;
    int chi_number;
    int weight;
    int ntraces;

    if(argc < 5) {
        cout << "usage: ./print-traces level weight chi ntraces prec [nthreads]" << endl;
        return 0;
    }
    init_classnumbers();
    load_factor_table();

    level = atoi(argv[1]);
    weight = atoi(argv[2]);
    chi_number = atoi(argv[3]);
    ntraces = atoi(argv[4]);

    int prec = atoi(argv[5]);
    int nthreads = 1;
    if(argc > 6) nthreads = atoi(argv[6]);
    int verbose = 0;

    DirichletGroup G(level, prec);
    if(GCD(level, chi_number) != 1) return 0;
    DirichletCharacter chi = G.character(chi_number);
    if(chi.is_even() && weight % 2 == 1) return 0;
    if(!chi.is_even() && weight % 2 == 0) return 0;

    cuspforms_acb * S = get_cuspforms_acb(chi, weight, nthreads, verbose);
    S->compute_traces(ntraces + 1);
    for(int k = 1; k <= ntraces; k++) {
        cout << k << " ";
        acb_printd(S->traces[k], 10);
        cout << endl;
    }

    return 0;
}
