#include <iostream>
#include <string>

#include "flint/fmpz_poly.h"

#include "mfformat.h"
#include "flint-extras.h"

using namespace std;

int main(int argc, char ** argv) {
    if(argc < 2) {
        cout << "usage: " << argv[0] << " dbname [fmt]" << endl;
        return 0;
    }
    int mfdatafmt = 0;
    if(argc == 3) {
        mfdatafmt = 1;
    }
    sqlite3 * db;
    sqlite3_open(argv[1], &db);


    const char * query_sql = "SELECT level, weight, chi, chiorbit, labelnumber, polynomial, traces from heckepolys"
                             " ORDER BY level, weight, chiorbit, labelnumber";

    sqlite3_stmt * stmt;

    sqlite3_prepare_v2(db, query_sql, -1, &stmt, NULL);
    int result = sqlite3_step(stmt);
    if(result != SQLITE_ROW) {
        cout << result << endl;
        return 1;
    }
    fmpz_poly_t f;
    fmpz_poly_init(f);
    while(result == SQLITE_ROW) {
        int level = sqlite3_column_int(stmt, 0);
        int weight = sqlite3_column_int(stmt, 1);
        int chi = sqlite3_column_int(stmt, 2);
        int chiorbit = sqlite3_column_int(stmt, 3);
        int labelnumber = sqlite3_column_int(stmt, 4);

        size_t polysize = sqlite3_column_bytes(stmt, 5);
        const void * polydata = sqlite3_column_blob(stmt, 5);

        x_fmpz_poly_read_raw(f, polydata, polysize);

        size_t tracesize = sqlite3_column_bytes(stmt, 6);
        const void * tracedata = sqlite3_column_blob(stmt, 6);
        long ntraces;
        fmpz * traces;
        x_fmpz_vec_read_raw(&ntraces, &traces, tracedata, tracesize);

        if(mfdatafmt) {
            cout << level << ":" << weight << ":" << chiorbit + 1 << ":[" << traces;
            for(int k = 1; k < ntraces; k++) {
                cout << ", " << traces + k;
            }
            cout << "]" << endl;
        }

        else {
            cout << level << " " << weight << " " << chi << " " << chiorbit << " " << labelnumber << ";";
            for(long k = 0; k < ntraces; k++) {
                cout << " " << traces + k;
            }
            cout << "; " << f << endl;
        }

        _fmpz_vec_clear(traces, ntraces);
        result = sqlite3_step(stmt);
    }
    sqlite3_finalize(stmt);

    return 0;
}
