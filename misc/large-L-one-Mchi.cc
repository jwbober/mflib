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

int main(int argc, char ** argv) {
    cout << setprecision(17);
    if(argc < 2) {
        cout << "Usage: l-one q" << endl;
        return 1;
    }
    int q = atoi(argv[1]);
    
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

    // a test...
    //for(int j = 2; j < q; j++) {
    //    double error = abs(gauss_sums[j] - G.character(j).gauss_sum());
    //    maxerror = max(error, maxerror);
    //    cout << j << " " << error << " " << maxerror << endl;
    //}

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

    //ofstream even_out("l-one-" + to_string(q) + "-even");
    //ofstream odd_out("l-one-" + to_string(q) + "-odd");

    //ofstream even_out("large-l-one-mchi-" + to_string(q) + "-even");
    ofstream odd_out("large-l-one-mchi-" + to_string(q) + "-odd");

    for(unsigned int j = 2; j < q; j++) {
        if(GCD(j, q) > 1) continue;
        if(InvMod(j, q) < j) continue;

        DirichletCharacter chi = G.character(j);
        if(!chi.is_primitive()) continue;
        if(chi.is_even()) continue;
        
        double x = abs(Lone_values[j]);
        double y = arg(Lone_values[j]);

        if(x > 2.5 * expgamma) {
            //cout << j << " " << x << endl;
            //if(chi.is_even()) {
            //    even_out << j << " " << x << endl;
            //}
            odd_out << j << " " << x << " " << abs(chi.max(0)) << endl;
            //odd_out << j << " " << x << endl;
        }
        //if(chi.is_even()) {
        //    even_out.write((const char *)&j, sizeof(j));
        //    even_out.write((const char *)&x, sizeof(x));
        //    even_out.write((const char *)&y, sizeof(y));
        //}
        //else {
        //    odd_out.write((const char *)&j, sizeof(j));
        //    odd_out.write((const char *)&x, sizeof(x));
        //    odd_out.write((const char *)&y, sizeof(y));
        //}
    }
    //even_out.close();
    odd_out.close();
    return 0;
}
