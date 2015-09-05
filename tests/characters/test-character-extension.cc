#include <iostream>
#include <cstdlib>

#include "characters.h"

int main(int argc, char ** argv) {
    long q = atol(argv[1]);
    DirichletGroup G(q);
    for(long q0 = 1; q0 <= q; q0++) {
        if(q % q0 != 0) continue;
        DirichletGroup G0(q0);
        for(int n = 1; n <= q0; n++) {
            if(GCD(n, q0) != 1) continue;
            if(!G0.character(n).is_primitive()) continue;
            long primitive_index = n;
            long induced_index = G.index_from_primitive_character(q0, n);
            long primitive_index2;
            long conductor2 = G.character(induced_index).conductor(&primitive_index2);
            if(conductor2 != q0) {
                cout << q0 << " " << n << " "
                     << conductor2 << " " << primitive_index2
                     << " conductor was wrong" << endl;
            }
            if(primitive_index2 != primitive_index && q0 != 1) {
                cout << "index was wrong" << endl;
            }
        }
    }
    return 0;
}
