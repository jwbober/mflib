#include <iostream>

#include "characters.h"

using namespace std;

int main(int argc, char ** argv) {
    if(argc < 2) {
        cout << "Usage:" << endl;
        cout << "./chipartialsums q a b k" << endl;
        return 0;
    }
    long q = atol(argv[1]);
    long a = atol(argv[2]);
    long b = atol(argv[3]);
    long n = (q*a)/b;
    long k = atol(argv[4]);

    complex<double> * sums = new complex<double>[q];
    DirichletGroup G(q);
    G.all_sums(sums, n);
    double S = 0;
    for(int j = 2; j < q; j++) {
        if(G.chi(j, q-1).real() - 1 > -.5) {
            S = S + pow(abs(sums[j]), 2*k)/q;
        }
    }

    cout << S << endl;

    return 0;
}
