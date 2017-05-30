#include <iostream>

#include "characters.h"

using namespace std;

int main(int argc, char ** argv) {
    if(argc < 2) {
        cout << "Usage:" << endl;
        cout << "./chimaxsum q [chi]" << endl;
        return 0;
    }
    long q = atol(argv[1]);
    long n = 0;
    if(argc > 2)
        n = atol(argv[2]);

    n %= q;
    DirichletGroup G(q);
    if(n < 0) {
        n += q;
    }
    if(GCD(q, n) == 1) {
        long index;
        complex<double> max_sum = G.character(n).max(&index);
        cout << q << " " << n << " " << abs(max_sum) << " " << abs(max_sum) << " " << index << endl;
    }
    else {
        for(int k = 0; k < q; k++) if(GCD(k, q) == 1){
            long index;
            complex<double> max_sum = G.character(k).max(&index);
            cout << q << " " << k << " " << abs(max_sum) << " " << abs(max_sum) << " " << index << endl;

        }
    }
    return 0;
}
