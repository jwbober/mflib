#include <iostream>
#include <cstdlib>
#include <vector>

#include "slint.h"


using namespace std;

int main(int argc, char ** argv) {
    long end = atol(argv[1]);
    vector<long> primes;
    prime_range(&primes, end);
    long p = 2;
    for(auto i = primes.begin(); i != primes.end(); i++) {
        if(*i != p) cout << *i << " " << p << endl;
        p = next_prime(p);
    }
}
