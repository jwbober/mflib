#include <iostream>
#include <cstdlib>
#include "characters.h"

using namespace std;

int main(int argc, char ** argv) {
    if(argc < 2) {
        cout << "missing arguments." << endl;
        return -1;
    }
    int q = atoi(argv[1]);
    DirichletGroup G(q);
    G.print_dft_translation();
    return 0;
}
