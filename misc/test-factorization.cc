#include <iostream>
#include "slint.h"

using namespace std;

int main(int argc, char ** argv) {
    load_factor_table();
    int n = atoi(argv[1]);
    int_factorization_t fac;
    factor(n, fac);
    for(int k = 0; k < fac.nfactors; k++) {
        cout << fac.factors[k].p << "^" << fac.factors[k].e;
        if(k < fac.nfactors - 1) cout << " * ";
    }
    cout << endl;
    return 0;
}
