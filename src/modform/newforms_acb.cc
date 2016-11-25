#include "modform_acb.h"
#include "arb_extras.h"
//#include "modform_cc.h"
#include "classnumbers.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>

using namespace std;

/*
complex<double> Lsign(const Eigen::Matrix<complex<double>, 1, Eigen::Dynamic> &coeffs, int level, int weight, int verbose) {
    complex<double> S1 = 0.0, S2 = 0.0, S3 = 0.0, S4 = 0.0;
    complex<double> A = 1.0;
    complex<double> B = 1.1;
    for(int k = 1; k < coeffs.size(); k++) {
        complex<double> an = coeffs[k];
        double n = k;
        double npow = pow(n, weight - 1);
        S1 += an/npow * exp(-2 * M_PI * n * A/sqrt(level));
        S2 += conj(an)/npow * exp(-2 * M_PI * n/(A*sqrt(level)));
        S3 += an/npow * exp(-2 * M_PI * n * B/sqrt(level));
        S4 += conj(an)/npow * exp(-2 * M_PI * n /(B*sqrt(level)));
    }
    complex<double> sign = (S1 - S3)/(S4 - S2);
    if(verbose)
        cerr << "computing sign = " << sign << ", |sign| = " << abs(sign) << endl;
    return sign;
}


static complex<double> acb_get_z(const acb_t in) {
    double x = arf_get_d(arb_midref(acb_realref(in)), ARF_RND_NEAR);
    double y = arf_get_d(arb_midref(acb_imagref(in)), ARF_RND_NEAR);
    return complex<double>(x,y);
}

void acb_sub_z(acb_t out, acb_t s, complex<double> z, int prec) {
    arb_t x;
    arb_init(x);

    arb_set_d(x, z.real());
    arb_sub(acb_realref(out), acb_realref(s), x, prec);

    arb_set_d(x, z.imag());
    arb_sub(acb_imagref(out), acb_imagref(s), x, prec);
    arb_clear(x);
}

static void acb_mat_printz(const acb_mat_t in) {
    int r = acb_mat_nrows(in);
    int c = acb_mat_ncols(in);
    cmatrix_t z(r, c);
    for(int j = 0; j < r; j++) {
        for(int k = 0; k < c; k++) {
            z(j, k) = acb_get_z(acb_mat_entry(in, j, k));
        }
    }
    cout << z << endl;
}

static cmatrix_t acb_mat_get_zmat(const acb_mat_t in) {
    int r = acb_mat_nrows(in);
    int c = acb_mat_ncols(in);
    cmatrix_t z(r, c);
    for(int j = 0; j < r; j++) {
        for(int k = 0; k < c; k++) {
            z(j, k) = acb_get_z(acb_mat_entry(in, j, k));
        }
    }
    return z;
}
*/

int main(int argc, char ** argv) {
    int level;
    int chi_number;
    int weight;
    int ncoeffs;

    if(argc < 7) {
        cout << "usage: ./newforms_acb level weight chi ncoeffs outpath prec [verbose]" << endl;
        return 0;
    }
    init_classnumbers();

    level = atoi(argv[1]);
    weight = atoi(argv[2]);
    chi_number = atoi(argv[3]);
    ncoeffs = atoi(argv[4]);
    string outpath(argv[5]);

    int prec = atoi(argv[6]);
    int verbose = 0;
    if(argc > 7) verbose = atoi(argv[7]);

    DirichletGroup G(level, prec);
    if(GCD(level, chi_number) != 1) return 0;
    DirichletCharacter chi = G.character(chi_number);
    if(chi.is_even() && weight % 2 == 1) return 0;
    if(!chi.is_even() && weight % 2 == 0) return 0;

    cuspforms_acb * S = get_cuspforms_acb(chi, weight, verbose);

    acb_mat_t newforms;
    int dim = S->new_dimension();
    cout << "computing a " << dim << " dimensional space." << endl;
    if(dim == 0) return 0;
    S->newforms(newforms, ncoeffs);
    acb_t z;
    acb_init(z);

    cout << "verifying multiplicativity." << endl;
#define newforms(i,j) acb_mat_entry(newforms, i, j)
    for(int k = 0; k < dim; k++) {
        for(int n = 1; n < ncoeffs; n++) {
            for(int m = 1; m*n < ncoeffs; m++) {
                if(GCD(m,n) != 0) continue;
                acb_mul(z, newforms(n, k), newforms(m, k), prec);
                if(!acb_overlaps(z, newforms(n*m, k))) {
                    cout << m << " " << n << " " << k <<  " multiplicativity check failed." << endl;
                    acb_printd(newforms(m,k), 10);
                    cout << endl;
                    acb_printd(newforms(n,k), 10);
                    cout << endl;
                    acb_printd(newforms(n*m,k), 10);
                    cout << endl;
                    acb_printd(z, 10);
                    cout << endl;
                    return 1;
                }
           }
        }
    }
    //cmatrix_t newforms2 = Scc->newforms(ncoeffs);

    //cout << (newforms1 - newforms2).transpose() << endl;
    //cout << endl;
    //cout << newforms1.transpose() << endl;
    //cout << endl;

    string fulloutpath = outpath + "/" + to_string(level)
                               + "/" + to_string(weight);
    string command = "mkdir -p " + fulloutpath;
    system(command.c_str());
    for(int k = 0; k < dim; k++) {
        string outfilename = fulloutpath + "/" 
                                   + to_string(level)
                                   + "." + to_string(weight)
                                   + "." + to_string(chi_number)
                                   + "." + to_string(k);
        
        cout << "writing data to " << outfilename << endl;
        FILE * outfile = fopen(outfilename.c_str(), "w");
        for(int j = 1; j < ncoeffs; j++) {
            fprintacb(outfile, newforms(j,k));
            fprintf(outfile, "\n");
            //char * realstring = arb_get_str(acb_realref(newforms(j, k)), prec/3, ARB_STR_NO_RADIUS);
            //char * imagstring = arb_get_str(acb_imagref(newforms(j, k)), prec/3, ARB_STR_NO_RADIUS);
            //cout << "(" << realstring << ", " << imagstring << ")" << endl;
            //free(realstring);
            //free(imagstring);
        }
        fclose(outfile);
    }
    //acb_mat_printd(newforms, 10);
    //cout << endl;
    //cmatrix_t newforms1 = acb_mat_get_zmat(newforms);
    //newforms2.transposeInPlace();

    /*
    for(int k = 0; k < ncoeffs; k++) {
        complex<double> S = 0;
        for(int j = 0; j < dim; j++) {
            S = S + newforms1(k, j) - newforms2(k, j);
        }
        cout << k << " " << abs(S) << endl;
    }
    */

    //cout << endl;
    //acb_mat_printz(newforms);
    //cout << endl;
    //cout << newforms2 << endl;

#undef newforms
    return 0;
}
