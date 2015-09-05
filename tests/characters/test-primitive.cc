#include <iostream>
#include "characters.h"

using namespace std;

int main() {
    int q = 11;
    cout << "q: ";
    cin >> q;

    DirichletGroup G(q);
    for(int n = 0; n < q; n++) {
        if(!G.is_coprime_to_q(n))
            continue;

        DirichletCharacter chi = G.character(n);
        if(chi.is_primitive())
            cout << n << endl;
    }
    return 0;
}
