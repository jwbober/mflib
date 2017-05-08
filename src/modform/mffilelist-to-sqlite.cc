#include <iostream>
#include <fstream>
#include <string>

#include <cstdio>

#include "mfformat.h"

using namespace std;


int insert_v2file_into_sqlite(sqlite3 * db, string filename) {
    // expect a database with a table modforms
    // CREATE TABLE modforms (level INTEGER, weight INTEGER, chi INTEGER, orbit INTEGER, j INTEGER,
    //      prec INTEGER, exponent INTEGER, ncoeffs INTEGER, coefficients BLOB)

    FILE * infile = fopen(filename.c_str(), "r");
    fseek(infile, 0, SEEK_END);
    size_t filesize = ftell(infile);
    fseek(infile, 0, SEEK_SET);

    struct mfheader header;
    read_mfheader(infile, &header);

    if(header.version != MFV2) {
        return SQLITE_ERROR;
    }

    size_t filepos = ftell(infile);
    size_t datasize = filesize - filepos;

    char * coeffdata = (char *)malloc(datasize);
    size_t bytes_read = fread((void *)coeffdata, 1, datasize, infile);
    fclose(infile);
    if(bytes_read != datasize) {
        cout << "ohno";
        return 0;
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
        if(insert_v2file_into_sqlite(db, filename) != SQLITE_DONE) {
            cout << "error inserting data." << endl;
            return 1;
        }
        count++;
        getline(filelist_infile, filename);
    }
    sqlite3_exec(db, "END TRANSACTION", NULL, 0, NULL);
    return 0;
}
