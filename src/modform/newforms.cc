#include "modform_cc.h"
#include "classnumbers.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

complex<double> Lsign(const Eigen::Matrix<complex<double>, 1, Eigen::Dynamic> &coeffs, int level, int verbose) {
    complex<double> S1 = 0.0, S2 = 0.0, S3 = 0.0, S4 = 0.0;
    complex<double> A = 1.0;
    complex<double> B = 1.1;
    for(int k = 1; k < coeffs.size(); k++) {
        complex<double> an = coeffs[k];
        double n = k;
        S1 += an/n * exp(-2 * M_PI * n * A/sqrt(level));
        S2 += conj(an)/n * exp(-2 * M_PI * n/(A*sqrt(level)));
        S3 += an/n * exp(-2 * M_PI * n * B/sqrt(level));
        S4 += conj(an)/n * exp(-2 * M_PI * n /(B*sqrt(level)));
    }
    complex<double> sign = (S1 - S3)/(S4 - S2);
    if(verbose)
        cerr << "computing sign = " << sign << ", |sign| = " << abs(sign) << endl;
    return sign;
}

int main(int argc, char ** argv) {
    int level;
    int chi_number;
    int weight;
    int ncoeffs;

    init_classnumbers();

    level = atoi(argv[1]);
    chi_number = atoi(argv[2]);
    weight = atoi(argv[3]);
    ncoeffs = atoi(argv[4]);

    int verbose = 0;
    if(argc > 5) verbose = atoi(argv[5]);

    DirichletGroup G(level);
    if(GCD(level, chi_number) != 1) return 0;
    DirichletCharacter chi = G.character(chi_number);
    if(chi.is_even() && weight % 2 == 1) return 0;
    if(!chi.is_even() && weight % 2 == 0) return 0;

    cuspforms_cc * S = get_cuspforms_cc(chi, weight, verbose);
    cmatrix_t newforms = S->newforms(ncoeffs);

    if(verbose > 1) {
        cout << S->newspace_basis(ncoeffs) << endl;
        cout << endl;
        cout << newforms << endl << endl;
    }
    int dim = newforms.rows();
    for(int k = 0; k < dim; k++) {
        string filename = "mf/" + to_string(level) + "." + to_string(chi_number) + "." + to_string(k) + ".lcalc";
        ofstream lcalc_file(filename);
        lcalc_file << setprecision(17);
        lcalc_file  << 3 << endl
                    << 2 << endl
                    << ncoeffs << endl
                    << 0 << endl
                    << 1 << endl
                    << 1 << endl
                    << .5 << " " << 0 << endl
                    << sqrt(level)/(2*M_PI) << endl;
        complex<double> sign = Lsign(newforms.row(k), level, verbose);
        //if(chi_number == 1) {
        //    sign = -newforms(k, level - 1);
        //}
        //else {
        //    sign = -chi.gauss_sum()/newforms(k, level - 1);
        //}
        lcalc_file  << real(sign) << " " << imag(sign) << endl
                    << 0 << endl;
        for(int j = 1; j < ncoeffs; j++) {
            complex<double> z = newforms(k,j)/sqrt(j);
            lcalc_file << real(z) << " " << imag(z) << endl;
        }
        lcalc_file.close();
    }


    return 0;
}
