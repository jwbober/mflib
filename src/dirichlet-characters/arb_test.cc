#include <iostream>
#include <iomanip>

#define USE_ARB
#include "characters.h"

using namespace std;
int main(int argc, char ** argv) {
    int q = atoi(argv[1]);
    int n = atoi(argv[2]);

    cout << setprecision(10);
    DirichletGroup G(q);
    acb_t z;
    acb_t x;
    arb_t a;
    acb_init(z);
    acb_init(x);
    arb_init(a);
    for(int k = 0; k < q; k++) {
        G.chi(z, n, k);
        complex<double> zz = G.chi(n, k);
        acb_set_d_d(x, zz.real(), zz.imag());
        acb_sub(x, z, x, 100);
        acb_abs(a, x, 100);
        double b = arf_get_d(arb_midref(a), ARF_RND_NEAR);
        if(b > 1e-14) cout << k << " " << b << endl;
        //arb_printd(a, 10);
        //cout << endl;
    }
    return 0;
}
