#include <iostream>
#include <fstream>
#include "slint.h"

using namespace std;

//struct int_factor_t {
//    int p;
//    int e;
//    int f; // f == p^e
//};

int_factor_t * _factor_table;
int _factor_table_size = 0;

void _build_factor_table(int size) {
    _factor_table = new int_factor_t[size];
    for(int k = 0; k < size; k++) {
        _factor_table[k] = {1,1,1};
    }
    std::vector<long> primes;
    prime_range(&primes, size);
    for(auto p : primes) {
        int e = 1;
        long f = p;
        while(f < size) {
            for(int n = f; n < size; n += f) {
                _factor_table[n] = {(int)p, e, (int)f};
            }
            f = f * p;
            e++;
        }
    }
    _factor_table_size = size;
}

int main(int argc, char ** argv) {
    int size = atoi(argv[1]);
    _build_factor_table(size);
    ofstream outfile(argv[2]);
    outfile.write( (char*) _factor_table, (long)_factor_table_size * (long)sizeof(int_factor_t));
    outfile.close();
    return 0;
}
