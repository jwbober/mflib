#include <iostream>
#include <cstdlib>
#include <vector>

#include "slint.h"


using namespace std;

int main(int argc, char ** argv) {
    long end = atol(argv[1]);
    vector<long> squarefrees;
    squarefree_range(&squarefrees, end);
    long n = 1;
    for(auto i = squarefrees.begin(); i != squarefrees.end(); i++) {
        if(*i != n) cout << *i << " " << n << endl;
        do {
            n++;
        } while(!is_squarefree(n));
    }
}
