#include <iostream>
#include <cstdlib>
#include <ctime>

#include "slint.h"

using namespace std;

int main() {
    srand(time(NULL));
    for(int k = 0; k < 1000; k++) {
        vector<long> a;
        vector<long> n;
        long w = rand();
        long q = 1;
        long p = rand() % 500;
        for(int j = 0; j < 6; j++) {
            p = next_prime(p);
            q = q * p;
            n.push_back(p);
            a.push_back(w % p);
        }
        long w1 = CRT(a, n);
        if(w % q != w1) {
            cout << k << " failure " << w << " " << w1 << " " << q << endl;
        }
    }
    return 0;
}
