#include <iostream>
#include "characters.h"

using namespace std;

int main() {
    int q = 11;
    int m = 3;
    cout << "q: ";
    cin >> q;
    cout << "m: ";
    cin >> m;
    
    cout << q << endl;
    DirichletGroup G(q);
    for(int n = 0; n < q; n++) {
        cout << n << " " << G.chi(m,n) << endl;
    }
    return 0;
}
