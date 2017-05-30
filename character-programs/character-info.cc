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
    if(argc < 2) {
        cout << "Usage:" << endl
             << "./character-info level [chi]" << endl;
        return 0;
    }
    long q = atol(argv[1]);
    long n = 0;
    if(argc > 2)
        n = atol(argv[2]);

    DirichletGroup G(q);
    if(n < 0) {
        n %= q;
        n += q;
    }
    if(GCD(q, n) == 1) {
        print_character_info(G, n);
    }
    else {
        auto orbits = G.galois_orbits();
        cout << "For modulus " << q << ", there are " << orbits.size() << " Galois orbits." << endl;
        cout << "Printing out some information about a representative from each orbit:" << endl;
        for(auto orbit : orbits) {
            long n = *(orbit.begin());
            print_character_info(G, n);
            cout << endl << "------------------------------------------------------------" << endl << endl;
        }
    }
    return 0;
}
