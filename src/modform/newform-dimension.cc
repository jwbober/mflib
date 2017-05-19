#include "modform_acb.h"
#include "arb_extras.h"
//#include "modform_cc.h"
#include "classnumbers.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>

#include "mag.h"

#include "mfformat.h"

using namespace std;

int main2(int argc, char ** argv) {
    int level1 = atoi(argv[1]);
    int level2 = atoi(argv[2]);
    int weight1 = atoi(argv[3]);
    int weight2 = atoi(argv[4]);
    int prec = 500;

    init_classnumbers();
    load_factor_table();

    for(int level = level1; level <= level2; level++) {
        DirichletGroup G(level, prec);
        for(int weight = weight1; weight <= weight2; weight++) {
            for(int n = 1; n < (level == 1 ? 2 : level); n++) {
                if(GCD(level, n) != 1) continue;
                DirichletCharacter chi = G.character(n);
                if(chi.is_even() && weight % 2 == 1) continue;
                if(!chi.is_even() && weight % 2 == 0) continue;
                int dimension = get_cuspforms_acb(chi, weight)->new_dimension();
                if(dimension != 0)
                    cout << level << " " << weight << " " << n << " " << dimension << endl;
            }
        }
        if(level % 10 == 0) {
            clear_cuspform_cache();
            clear_cuspforms_modp();
        }
    }

    return 0;
}

int main(int argc, char ** argv) {
    if(argc > 4) {
        return main2(argc, argv);
    }
    int level;
    int chi_number;
    int weight;

    init_classnumbers();
    load_factor_table();

    level = atoi(argv[1]);
    weight = atoi(argv[2]);
    if(argc > 3)
        chi_number = atoi(argv[3]);
    else
        chi_number = 0;
    int prec = 500;
    DirichletGroup G(level, prec);
    if(chi_number == 0) {
        for(int n = 1; n < (level == 1 ? 2 : level); n++) {
            if(GCD(level, n) != 1) continue;
            DirichletCharacter chi = G.character(n);
            if(chi.is_even() && weight % 2 == 1) continue;
            if(!chi.is_even() && weight % 2 == 0) continue;
            int dimension = get_cuspforms_acb(chi, weight)->new_dimension();
            if(dimension != 0)
                cout << level << " " << weight << " " << n << " " << dimension << endl;
        }
    }
    else {
        if(GCD(level, chi_number) != 1) return 0;
        DirichletCharacter chi = G.character(chi_number);
        if(chi.is_even() && weight % 2 == 1) {cout << 0 << endl; return 0;}
        if(!chi.is_even() && weight % 2 == 0) {cout << 0 << endl; return 0;}
        cout << get_cuspforms_acb(chi, weight)->new_dimension() << endl;
    }

    return 0;
}
