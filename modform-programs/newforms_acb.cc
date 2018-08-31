#include "cuspforms_acb.h"
#include "arb-extras.h"
#include "classnumbers.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>

#include "mag.h"

#include "mfformat.h"

using namespace std;

long arb_abs_prec_approx(arb_t x) {
    if(arb_is_exact(x)) return ARF_PREC_EXACT;
    double mag_bits = mag_get_d_log2_approx(arb_radref(x));
    return floor(mag_bits);
}

int main(int argc, char ** argv) {
    srand(time(NULL));
    int level;
    int chi_number;
    int weight;
    int ncoeffs;

    if(argc < 7) {
        cout << "usage: ./newforms_acb level weight chi ncoeffs outpath prec targetprec [nthreads] [verbose]" << endl;
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
    int nthreads = 1;
    if(argc > 8) nthreads = atoi(argv[8]);
    int verbose = 0;
    if(argc > 9) verbose = atoi(argv[9]);

    DirichletGroup G(level, prec);
    if(GCD(level, chi_number) != 1) return 0;
    DirichletCharacter chi = G.character(chi_number);
    if(chi.is_even() && weight % 2 == 1) return 0;
    if(!chi.is_even() && weight % 2 == 0) return 0;

    cuspforms_acb * S = get_cuspforms_acb(chi, weight, nthreads, verbose);

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
    sqlite3_exec(db, "DELETE FROM modforms;", NULL, 0, NULL);

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

        if(computed_precision > 0) {
            cerr << "error1: We don't have enough precision, so aborting." << endl;
            sqlite3_exec(db, "ROLLBACK TRANSACTION", NULL, 0, NULL);
            return 1;
        }

        if(verbose && (computed_precision > targetprec)) {
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

        int precision_ok = 1;
        do {
            if(header.prec != MF_PREC_EXACT && header.prec > 0) {
                cerr << "error2: We don't have enough precision, so aborting." << endl;
                cerr << "precision was " << header.prec << endl;
                sqlite3_exec(db, "ROLLBACK TRANSACTION", NULL, 0, NULL);
                return 1;
            }
            precision_ok = 1;
            int j = 0;
            int pp = prime_powers_table[j];

            while(pp <= header.ncoeffs) {
                int ok_to_round = acb_attempt_rounding_mfcoeff(&header, newforms(pp,k));
                if(!ok_to_round) {
                    //acb_printd(newforms(pp, k), 10); cout << endl;
                    precision_ok = 0;
                    if(header.prec == MF_PREC_EXACT) {
                        cerr << "error: precision should be exact but we are having trouble rounding." << endl;
                        cerr << "This should not happen. Quitting." << endl;
                        sqlite3_exec(db, "ROLLBACK TRANSACTION", NULL, 0, NULL);
                        return 1;
                    }
                    header.prec++;
                    header.exponent++;
                    if(verbose) {
                        cout << "Rounding failed. Changing precision to " << header.prec << endl;
                    }
                    break;
                }
                j++;
                pp = prime_powers_table[j];
            }
        } while(!precision_ok);


        char * coeffblob;
        size_t coeffsize;
        FILE * coeff_psuedofile = open_memstream(&coeffblob, &coeffsize);

        int j = 0;
        int pp = prime_powers_table[j];
        while(pp <= header.ncoeffs) {
            long bytes_written = acb_write_mfcoeff(coeff_psuedofile, &header, newforms(pp,k));
            if(bytes_written == 0) {
                cerr << "error: no bytes written." << endl;
                sqlite3_exec(db, "ROLLBACK TRANSACTION", NULL, 0, NULL);
                return 1;
            }
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
