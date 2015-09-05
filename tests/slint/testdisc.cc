#include <iostream>
#include <cstdlib>
#include <vector>

#include "slint.h"


using namespace std;

int main(int argc, char ** argv) {
    long end = atol(argv[1]);
    vector<long> discs;
    discriminant_range(&discs, end);
    long n = 1;
    for(auto i = discs.begin(); i != discs.end(); i++) {
        if(*i != n) cout << *i << " " << n << endl;
        do {
            n++;
        } while(!is_fundamental_discriminant(n));
    }
}
