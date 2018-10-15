#include <iostream>

#include "cuspforms_weight1_modp.h"
#include "classnumbers.h"

using namespace std;

int main(int argc, char ** argv) {
    init_classnumbers();
    load_factor_table();

    if(argc < 2) {
        cout << "usage: " << argv[0] << " level [p0] [chi]" << endl
             << endl
             << "For each character mod level, or for just the character chi, this program" << endl
             << "prints out a line" << endl
             << endl
             << "level chi p dimension" << endl
             << endl
             << "where dimension is (an upper bound for) the dimension of the space of" << endl
             << "weight 1 cuspforms mod p and p is an appropriately chosen prime >= p0" << endl
             << endl
             << "WARNING: Don't expect this to work for small p0." << endl;
         return 0;
    }

    int level = atoi(argv[1]);
    int chinumber = 0;
    int p0 = 10000;
    if(argc > 2) {
        p0 = atoi(argv[2]);
    }
    if(argc > 3) {
        chinumber = atoi(argv[3]);
    }
    DirichletGroup G(level);
    if(chinumber == 0) {
        for(int k = 1; k < level; k++) {
            if(GCD(k, level) != 1) continue;
            DirichletCharacter chi = G.character(k);
            cuspforms_weight1_modp * S = get_cuspforms_weight1_modp(chi, p0, 2);
            int p = S->p;
            cout << level << " " << k << " " << p << " " << S->dimension() << endl;
        }
    }
    else {
        if(GCD(chinumber, level) != 1) {
            cout << level << " " << chinumber << " -1 0" << endl;
        }
        else {
            DirichletCharacter chi = G.character(chinumber);
            cuspforms_weight1_modp * S = get_cuspforms_weight1_modp(chi, p0, 2);
            int p = S->p;
            cout << level << " " << chinumber << " " << p << " " << S->dimension() << endl;
        }
    }
}
