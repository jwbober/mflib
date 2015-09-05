#include <iostream>
#include <cstdlib>
#include <ctime>
#include "characters.h"

using namespace std;



double seconds_since_last_call() {
    static clock_t lasttime = 0;
    clock_t thistime = clock();
    double seconds = (thistime - lasttime)/(double)CLOCKS_PER_SEC;
    lasttime = thistime;
    return seconds;
}

void timesum(long qstart, long qend) {
    complex<double> S = 0.0;
    for(int q = qstart; q < qend; q++) {
        cout << q << endl;
        complex<double> * sums = new complex<double>[q];
        DirichletGroup G(q);
        
        G.all_sums(sums, q/3);
        S += G.chi(2,3);
        delete [] sums;
    }
    cout << S << endl;
}

void mpfi_c_set_complex(mpfi_c_t out, complex<double> in) {
    mpfi_c_set_d(out, in.real(), in.imag());
}

void mpfi_c_sub_complex(mpfi_c_ptr out, mpfi_c_ptr z1, complex<double> z2) {
    mpfi_sub_d(out->re, z1->re, z2.real());
    mpfi_sub_d(out->im, z1->im, z2.imag());
}

inline void init(mpfr_t X) {mpfr_init(X);}
inline void init(mpfi_t X) {mpfi_init(X);}
inline void init(mpfi_c_t X) {mpfi_c_init(X);}

double time_mp_fft(int q) {
    seconds_since_last_call();
    DirichletGroup G(q);

    mpfi_c_t * sums = new mpfi_c_t[q];
    for(int k = 0; k < q; k++) {
        mpfi_c_init(sums[k]);
    }

    for(int k = 0; k < 100; k++) {
        G.all_sums(sums, q-1);
    }

    for(int k = 0; k < q; k++) {
        mpfi_c_clear(sums[k]);
    }
    delete [] sums;

    return seconds_since_last_call()/100;
}

double time_mp_dftsum(int q) {
    seconds_since_last_call();
    DirichletGroup G(q);

    mpfi_c_t * sums = new mpfi_c_t[q];
    for(int k = 0; k < q; k++) {
        mpfi_c_init(sums[k]);
    }

    unsigned long * X = new unsigned long[q]();

    for(int k = 0; k < q; k++) {
        X[k] = 3;
    }

    for(int k = 0; k < 20; k++) {
        G.DFTsum_direct(sums, X);
    }

    delete [] X;

    for(int k = 0; k < q; k++) {
        mpfi_c_clear(sums[k]);
    }
    delete [] sums;
    return seconds_since_last_call()/20;
}


int main() {
    //timesum( 1 << 14, 1 << 15);
    int qstart = 0;
    double maxerror = 0;
    dft_init(200);

    for(int q = 1000; q < 1050; q++) {
        double t1 = time_mp_fft(q);
        double t2 = time_mp_dftsum(q);

        cout << q << " " << t1 << " " << t2 << " " << t2/t1 << endl;
    }

    mpfr_t mpfr1, mpfr2;
    mpfi_c_t mpc1;
    mpfi_t mp1;
    init(mpc1);
    init(mp1);
    init(mpfr1);
    init(mpfr2);
    for(int j = 1; j < 1000; j++) {
        int q = qstart + j;
        DirichletGroup G(q);
        complex<double> * a = new complex<double>[q]();
        complex<double> * X1 = new complex<double>[q]();
        complex<double> * X2 = new complex<double>[q]();

        for(int k = 0; k < q; k++) {
            a[k] = complex<double>(rand()/(double)RAND_MAX, rand()/(double)RAND_MAX);
        }

        //G.DFTsum(X1, a);
        //G.DFTsum_direct(X2, a);
        
        double error = 0;
        //for(int k = 0; k < q; k++) {
            //cout << X1[k] << " " << X2[k] << endl;
        //    error = max(error, abs(X1[k] - X2[k]));
        //}
        //maxerror = max(error, maxerror);
        //cout << q << " " << error << " " << maxerror << endl;

        complex<double> * sums1 = new complex<double>[q];
        complex<double> * sums2 = new complex<double>[q]();
        mpfi_c_t * sums3 = new mpfi_c_t[q];
        for(int k = 0; k < q; k++) {
            mpfi_c_init(sums3[k]);
        }

        mpfi_c_t * sums4 = new mpfi_c_t[q];
        for(int k = 0; k < q; k++) {
            mpfi_c_init(sums4[k]);
        }

        G.all_sums(sums1, q/3);
        G.all_sums(sums3, q/3);

        unsigned long * X;
        X = new unsigned long[q]();
        for(int k = 0; k <= q/3; k++) {
            X[k] = 1;
        }

        G.DFTsum_direct(sums4, X);

        delete [] X;

        //for(int k = 0; k < q; k++) {
        //    if(G.is_coprime_to_q(k)) {
        //        DirichletCharacter chi = G.character(k);
        //        sums2[k] = chi.sum(q/3);
        //    }
        //}
        error = 0;
        for(int k = 0; k < q; k++) {
            //error = max(error, abs(sums1[k] - sums2[k]));
            mpfi_c_sub(mpc1, sums3[k], sums4[k]);
            mpfi_c_abs(mp1, mpc1);
            mpfi_get_right(mpfr1, mp1);
            //mpfi_c_print(sums3[k]);
            error = max(error, mpfr_get_d(mpfr1, GMP_RNDN));
        }
        maxerror = max(error, maxerror);
        cout << q << " " << error << " " << maxerror << endl;

        delete [] a;
        delete [] X1;
        delete [] X2;
        delete [] sums1;
        delete [] sums2;

        for(int k = 0; k < q; k++) mpfi_c_clear(sums3[k]);
        for(int k = 0; k < q; k++) mpfi_c_clear(sums4[k]);
        delete [] sums3;
        delete [] sums4;
    }
    return 0;
}
