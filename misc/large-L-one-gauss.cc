#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iomanip>
#include "characters.h"

using namespace std;

const complex<double> I = complex<double>(0, 1);

const double expgamma = 1.78107241799020;

void do_q(int q) {
    q = 3*q;
    DirichletGroup G(q);
    complex<double> * a = new complex<double>[q]();
    complex<double> * gauss_sums = new complex<double>[q]();
    complex<double> * even_sum1 = new complex<double>[q]();
    //complex<double> * odd_sum1 = new complex<double>[q]();
    complex<double> * Lone_values = new complex<double>[q]();

    // first compute the gauss sums

    for(int j = 0; j < q; j++) {
        a[j] = exp(2.0 * M_PI * I * (double)j/(double)q);
    }
    G.DFTsum(gauss_sums, a);
    double maxerror = 0;

    for(int j = 0; j < q; j++) {
        a[j] = log(sin(M_PI*j/(double)q));
    }
    G.DFTsum(even_sum1, a);

    //for(int j = 0; j < q; j++) {
    //    a[j] = (double)j/(double)q;
    //}
    //G.DFTsum(odd_sum1, a);

    for(int j = 2; j < q; j++) {
        if(GCD(j, q) == 1) {
            int inv = InvMod(j, q);
            if(G.character(j).is_even()) {
                Lone_values[j] = -1.0/gauss_sums[inv] * even_sum1[inv];
            }
            //else {
            //    Lone_values[j] = -I*M_PI/gauss_sums[inv] * odd_sum1[inv];
            //}
        }
    }

    int max_index = 0;
    double max = 0;


    for(unsigned int j = 2; j < q; j++) {
        if(GCD(j, q) > 1) continue;
        if(InvMod(j, q) < j) continue;

        DirichletCharacter chi = G.character(j);
        if(!chi.is_primitive()) continue;
        //if(chi.is_even()) continue;
        if(!chi.is_even()) continue;
        
        double x = abs(Lone_values[j]);

        if(x > max) {
            max_index = j;
            max = x;
        }
    }

    cout << q << " " << max_index << " " << arg(Lone_values[max_index]) << " " << arg(gauss_sums[max_index]) << endl;

    delete [] a;
    delete [] gauss_sums;
    //delete [] odd_sum1;
    delete [] even_sum1;
    delete [] Lone_values;
}


void do_q_odd(int q) {
    q = 3*q;
    DirichletGroup G(q);
    complex<double> * a = new complex<double>[q]();
    complex<double> * gauss_sums = new complex<double>[q]();
    //complex<double> * even_sum1 = new complex<double>[q]();
    complex<double> * odd_sum1 = new complex<double>[q]();
    complex<double> * Lone_values = new complex<double>[q]();

    // first compute the gauss sums

    for(int j = 0; j < q; j++) {
        a[j] = exp(2.0 * M_PI * I * (double)j/(double)q);
    }
    G.DFTsum(gauss_sums, a);
    double maxerror = 0;

    //for(int j = 0; j < q; j++) {
    //    a[j] = log(sin(M_PI*j/(double)q));
    //}
    //G.DFTsum(even_sum1, a);

    for(int j = 0; j < q; j++) {
        a[j] = (double)j/(double)q;
    }
    G.DFTsum(odd_sum1, a);

    for(int j = 2; j < q; j++) {
        if(GCD(j, q) == 1) {
            int inv = InvMod(j, q);
            if(G.character(j).is_even()) {
            //    Lone_values[j] = -1.0/gauss_sums[inv] * even_sum1[inv];
            }
            else {
                Lone_values[j] = -I*M_PI/gauss_sums[inv] * odd_sum1[inv];
            }
        }
    }

    int max_index = 0;
    double max = 0;


    for(unsigned int j = 2; j < q; j++) {
        if(GCD(j, q) > 1) continue;
        if(InvMod(j, q) < j) continue;

        DirichletCharacter chi = G.character(j);
        if(!chi.is_primitive()) continue;
        if(chi.is_even()) continue;
        //if(!chi.is_even()) continue;
        
        double x = abs(Lone_values[j]);

        if(x > max) {
            max_index = j;
            max = x;
        }
    }

    cout << q << " " << max_index << " " << arg(Lone_values[max_index]) << " " << arg(gauss_sums[max_index]) << endl;

    delete [] a;
    delete [] gauss_sums;
    delete [] odd_sum1;
    //delete [] even_sum1;
    delete [] Lone_values;
}



int main(int argc, char ** argv) {
    cout << setprecision(17);
    if(argc < 3) {
        cout << "Usage: l-one q_start q_end" << endl;
        return 1;
    }
    int q_start = atoi(argv[1]);
    int q_end = atoi(argv[2]);
    int q = next_prime(q_start - 1);

    while(q < q_end) {
        do_q_odd(q);
        q = next_prime(q);
    }
    return 0;
}
