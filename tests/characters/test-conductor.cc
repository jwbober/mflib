#include <iostream>
#include <cstdlib>

#include "characters.h"

using namespace std;

int main(int argc, char ** argv) {
    long q = atol(argv[1]);
    DirichletGroup G(q);
    for(int k = 1; k < q; k++) {
        if(GCD(k, q) != 1) continue;
        long q0;
        long n0;

        DirichletCharacter chi = G.character(k);
        q0 = chi.conductor(&n0);
        cout << k << " " << q0 << " " << n0 << endl;
    }
    return 0;
}
