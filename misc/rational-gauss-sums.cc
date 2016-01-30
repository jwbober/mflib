#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iomanip>
#include "characters.h"

using namespace std;

const complex<double> I = complex<double>(0, 1);


void rational_in_interval(double start, double end, int &x, int &y) {
    int a = 0;
    int b = 1;
    int c = 1;
    int d = 1;

    if(start == 0) {
        if(1.0/end < ceil(1/end)) {
            x = 1;
            y = (int)(ceil(1/end));
            return;
        }
        else {
            x = 1;
            y = (int)(ceil(1/end) + 1);
        }
    }

    double z = (a + c)/(double)(b + d);
    while(1) {
        if(start < z && z < end) {
            x = a+c;
            y = b+d;
            return;
        }
        else if(start < z) {
            c = a+c;
            d = b+d;
        }
        else {
            a = a+c;
            b = b+d;
        }
        z = (a+c)/(double)(b+d);
    }
}

void do_q(int q) {
    DirichletGroup G(q);
    complex<double> * gauss_sums = new complex<double>[q]();
    complex<double> * a = new complex<double>[q]();

    // first compute the gauss sums

    for(int j = 0; j < q; j++) {
        a[j] = exp(2.0 * M_PI * I * (double)j/(double)q);
    }
    G.DFTsum(gauss_sums, a);
    for(int k = 2; k < q; k++) {
        if(GCD(k, q) > 1) continue;
        if(!G.character(k).is_primitive()) continue;
        double angle = arg(gauss_sums[k])/(2*M_PI);
        if(angle < 0) angle += 1;
        int a,b;
        if(angle < 1e-14) {a = 0; b = 1;}
        else rational_in_interval(angle - .00001, angle + .00001, a, b);
        if(abs(angle - (double)a/b) < 1e-14)
            //if(!G.character(k).is_primitive())
            //    cout << q << " " << k << " " << order_mod(k, q) << " " << angle << " " << a << " " << b << endl;
            //else
            //cout << q << " " << k << " " << order_mod(k, q) << " " << angle << " " << a << " " << b << endl;
        {
            int N0 = G.character(k*k % q).conductor();
            if(!is_squarefree(q/N0))
                cout << q << " " << k << " " << q/N0 << endl;
                //cout << q << " " << k << " " << N0 << " " << q/N0 << " " << is_squarefree(q/N0) << " " << GCD(N0, q/N0) << " " << endl;
        }
    }
    delete [] gauss_sums;
    delete [] a;
}

int main(int argc, char ** argv) {
    cout << setprecision(17);
    if(argc < 2) {
        cout << "Usage: rational_gauss_sums qstart [qend]" << endl;
        return 1;
    }
    int qstart = atoi(argv[1]);
    int qend = qstart;
    if(argc > 2) qend = atoi(argv[2]);
    for(int q = qstart; q <= qend; q++)
        do_q(q);
    return 0;
}
