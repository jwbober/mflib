#include <iostream>

#include "characters.h"

using namespace std;

void print_character_info(DirichletGroup& G, long n) {
    long q = G.q;
    DirichletCharacter chi = G.character(n);

    long primitive_index;
    long conductor = chi.conductor(&primitive_index);
    cout << "Dirichlet character (" << q << ", " << chi.m << ") of conductor " << conductor << endl;
    cout << "induced by primitive character (" << conductor << ", " << primitive_index << ")." << endl;
    cout << "Order: " << order_mod(chi.m, q) << endl;
    cout << "Parity: " << (chi.is_even() ? "even" : "odd") << endl;
    cout << endl;
    cout << "Retrictions to all lower levels:" << endl;
    for(int M : divisors((int)q)) {
        if(M % conductor != 0) continue;
        DirichletGroup G0(M);
        DirichletCharacter psi = G0.character(G0.index_from_primitive_character(conductor, primitive_index));
        cout << "    (" << M << ", " << psi.m << ")" << endl;
    }
    cout << endl;
    auto orbit = chi.galois_orbit();
    cout << orbit.size() << " characters in Galois orbit:";
    for(long m : orbit) {
        cout << " " << m;
    }
    cout << endl;
}


int main(int argc, char ** argv) {
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
