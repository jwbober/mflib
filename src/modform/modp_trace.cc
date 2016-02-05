#include "modform_modp.h"
#include "classnumbers.h"

#include <iostream>

using namespace std;

int main(int argc, char ** argv) {
    int level;
    int chi_number;
    int weight;
    int p;
    int number_of_traces;

    init_classnumbers();

    level = atoi(argv[1]);
    chi_number = atoi(argv[2]);
    weight = atoi(argv[3]);
    p = atoi(argv[4]);
    number_of_traces = atoi(argv[5]);

    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chi_number);

    cuspforms_modp * S = get_cuspforms_modp(chi, weight, p);
    p = S->p;
    cout << p << endl;
    S->compute_traces(number_of_traces + 1);
    for(cuspforms_modp * S2: S->subspaces) {
        if(S2->level == 1) continue;
        cout << S2->level << "\t";
        for(int k = 1; k <= number_of_traces; k++) {
            cout << S2->traces[k];
            if(k < number_of_traces) cout << "\t";
        }
        cout << endl;
    }
    cout << S->level << endl << endl;
    for(int k = 1; k <= number_of_traces; k++) {
        long a = S->traces[k];
        if(a > p/2) a -= p;
        cout << a << endl;
        //cout << S->traces[k];
        //if(k < number_of_traces) cout << "\t";
    }
    cout << endl;
    return 0;
}
