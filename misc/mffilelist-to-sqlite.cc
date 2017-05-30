#include <iostream>
#include <fstream>
#include <string>

#include <cstdio>

#include "mfformat.h"
#include "slint.h"

using namespace std;


int insert_mffile_into_sqlite(sqlite3 * db, string filename) {
    // expect a database with a table modforms
    // CREATE TABLE modforms (level INTEGER, weight INTEGER, chi INTEGER, orbit INTEGER, j INTEGER,
    //      prec INTEGER, exponent INTEGER, ncoeffs INTEGER, coefficients BLOB)

    FILE * infile = fopen(filename.c_str(), "r");
    fseek(infile, 0, SEEK_END);
    size_t filesize = ftell(infile);
    fseek(infile, 0, SEEK_SET);

    struct mfheader header;
    read_mfheader(infile, &header);

    char * coeffdata;
    size_t datasize;

    if(header.version == MFV2) {
        size_t filepos = ftell(infile);
        datasize = filesize - filepos;
        coeffdata = (char *)malloc(datasize);
        size_t bytes_read = fread((void *)coeffdata, 1, datasize, infile);
        fclose(infile);
        if(bytes_read != datasize) {
            cout << "ohno";
            return SQLITE_ERROR;
        }
    }
    else if(header.version == MFV1) {
        FILE * coeff_pseudofile = open_memstream(&coeffdata, &datasize);

        fmpz_t x;
        fmpz_t y;
        fmpz_init(x);
        fmpz_init(y);
        int k = 0;
        int p = prime_powers_table[k];
        for(int j = 0; j < header.ncoeffs; j++) {
            fmpz_inp_raw(x, infile);
            if(header.chi != 1) fmpz_inp_raw(y, infile);
            if(j + 1 == p) {
                fmpz_out_raw(coeff_pseudofile, x);
                if(header.chi != 1) fmpz_out_raw(coeff_pseudofile, y);
                k = k + 1;
                p = prime_powers_table[k];
            }
        }
        fclose(coeff_pseudofile);
        fclose(infile);
    }
    else {
        cerr << filename << " not recognized file type." << endl;
        return SQLITE_ERROR;
    }

    sqlite3_stmt * stmt;
    string sql = "INSERT INTO modforms (level, weight, chi, orbit, j, prec, exponent, ncoeffs, coefficients) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);";

    int result = sqlite3_prepare_v2(db, sql.c_str(), sql.size(), &stmt, NULL);
    if(result != SQLITE_OK) return result;
    sqlite3_bind_int(stmt, 1, header.level);
    sqlite3_bind_int(stmt, 2, header.weight);
    sqlite3_bind_int(stmt, 3, header.chi);
    sqlite3_bind_int(stmt, 4, header.orbit);
    sqlite3_bind_int(stmt, 5, header.j);
    sqlite3_bind_int(stmt, 6, header.prec);
    sqlite3_bind_int(stmt, 7, header.exponent);
    sqlite3_bind_int(stmt, 8, header.ncoeffs);
    sqlite3_bind_blob(stmt, 9, (void *)coeffdata, datasize, free);

    result = sqlite3_step(stmt);
    sqlite3_finalize(stmt);
    return result;
}

int main(int argc, char ** argv) {
    ifstream filelist_infile(argv[1]);
    sqlite3 * db;
    sqlite3_open(argv[2], &db);
    sqlite3_exec(db, "CREATE TABLE modforms (level INTEGER, weight INTEGER, chi INTEGER, orbit INTEGER, j INTEGER, prec INTEGER, exponent INTEGER, ncoeffs INTEGER, coefficients BLOB)", NULL, 0, NULL);

    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, 0, NULL);

    string filename;
    getline(filelist_infile, filename);
    long count = 0;
    while(filelist_infile) {
        cout << count << ": inserting " << filename << endl;
        //if(count % 100 == 0) {
        //    sqlite3_exec(db, "END TRANSACTION", NULL, 0, NULL);
        //    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, 0, NULL);
        //}
        if(insert_mffile_into_sqlite(db, filename) != SQLITE_DONE) {
            cout << "error inserting data." << endl;
            return 1;
        }
        count++;
        getline(filelist_infile, filename);
    }
    sqlite3_exec(db, "END TRANSACTION", NULL, 0, NULL);
    return 0;
}
