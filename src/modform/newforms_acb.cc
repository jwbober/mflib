#include "modform_acb.h"
#include "arb_extras.h"
//#include "modform_cc.h"
#include "classnumbers.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>

#include "mag.h"

#include "mfformat.h"

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

long arb_abs_prec_approx(arb_t x) {
    if(arb_is_exact(x)) return ARF_PREC_EXACT;
    double mag_bits = mag_get_d_log2_approx(arb_radref(x));
    return floor(mag_bits);
}

int main(int argc, char ** argv) {
    int level;
    int chi_number;
    int weight;
    int ncoeffs;

    if(argc < 7) {
        cout << "usage: ./newforms_acb level weight chi ncoeffs outpath prec targetprec[verbose]" << endl;
        return 0;
    }
    init_classnumbers();
    load_factor_table();

    level = atoi(argv[1]);
    weight = atoi(argv[2]);
    chi_number = atoi(argv[3]);
    ncoeffs = atoi(argv[4]);
    string outpath(argv[5]);

    int prec = atoi(argv[6]);
    long targetprec = atol(argv[7]);
    int verbose = 0;
    if(argc > 8) verbose = atoi(argv[8]);

    DirichletGroup G(level, prec);
    if(GCD(level, chi_number) != 1) return 0;
    DirichletCharacter chi = G.character(chi_number);
    if(chi.is_even() && weight % 2 == 1) return 0;
    if(!chi.is_even() && weight % 2 == 0) return 0;

    cuspforms_acb * S = get_cuspforms_acb(chi, weight, verbose);

    acb_mat_t newforms;
    int dim = S->new_dimension();
    if(verbose)
        cout << "computing a " << dim << " dimensional space." << endl;
    if(dim == 0) return 0;
    S->newforms(newforms, ncoeffs);
    acb_t z;
    acb_init(z);

    if(verbose)
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
    if(system(command.c_str()) != 0) {
        cerr << "error creating output directory. continuing, but don't expect much." << endl;
    }

    string outfilename = fulloutpath + "/"
                           + to_string(level)
                           + "." + to_string(weight)
                           + "." + to_string(chi_number)
                           + ".mfdb";

    sqlite3 * db;
    sqlite3_open(outfilename.c_str(), &db);
    sqlite3_exec(db, "CREATE TABLE modforms (level INTEGER, weight INTEGER, chi INTEGER, orbit INTEGER, j INTEGER, prec INTEGER, exponent INTEGER, ncoeffs INTEGER, coefficients BLOB)", NULL, 0, NULL);
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, 0, NULL);

    for(int k = 0; k < dim; k++) {
        long computed_precision = -ARF_PREC_EXACT;
        for(int j = 1; j < ncoeffs; j++) {
            long nx = arb_abs_prec_approx(acb_realref(newforms(j, k)));
            long ny = arb_abs_prec_approx(acb_imagref(newforms(j, k)));
            if(abs(nx) == ARF_PREC_EXACT) nx = -nx;
            if(abs(ny) == ARF_PREC_EXACT) ny = -ny;
            nx = max(nx, ny);
            computed_precision = max(computed_precision, nx);
        }
        if(computed_precision > targetprec) {
            cout << "ohno" << endl;
            cout << "target precision was " << targetprec << endl;
            cout << "computed precision was " << computed_precision << endl;
        }
        mfheader header;
        header.version = MFV2;
        header.level = level;
        header.weight = weight;
        header.chi = chi_number;
        header.orbit = 0;
        header.j = k;
        long fileprec = max(targetprec, computed_precision);
        if(fileprec == -ARF_PREC_EXACT)
            fileprec = MF_PREC_EXACT;
        header.prec = fileprec;
        if(fileprec == MF_PREC_EXACT)
            header.exponent = 0;
        else
            header.exponent = header.prec;
        header.ncoeffs = ncoeffs - 1;

        char * coeffblob;
        size_t coeffsize;
        FILE * coeff_psuedofile = open_memstream(&coeffblob, &coeffsize);

        int j = 0;
        int pp = prime_powers_table[j];
        while(pp <= header.ncoeffs) {
            long bytes_written = acb_write_mfcoeff(coeff_psuedofile, &header, newforms(pp,k));
            if(bytes_written == 0) cout << "error: no bytes written." << endl;
            j++;
            pp = prime_powers_table[j];
        }
        fclose(coeff_psuedofile);

        sqlite3_stmt * stmt;
        string sql = "INSERT INTO modforms (level, weight, chi, orbit, j, prec, exponent, ncoeffs, coefficients) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);";

        sqlite3_prepare_v2(db, sql.c_str(), sql.size(), &stmt, NULL);
        sqlite3_bind_int(stmt, 1, header.level);
        sqlite3_bind_int(stmt, 2, header.weight);
        sqlite3_bind_int(stmt, 3, header.chi);
        sqlite3_bind_int(stmt, 4, header.orbit);
        sqlite3_bind_int(stmt, 5, header.j);
        sqlite3_bind_int(stmt, 6, header.prec);
        sqlite3_bind_int(stmt, 7, header.exponent);
        sqlite3_bind_int(stmt, 8, header.ncoeffs);
        sqlite3_bind_blob(stmt, 9, (void *)coeffblob, coeffsize, free);

        sqlite3_step(stmt);
        sqlite3_finalize(stmt);
    }

    sqlite3_exec(db, "END TRANSACTION", NULL, 0, NULL);

#undef newforms
    return 0;
}
