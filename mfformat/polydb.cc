#include "mfformat.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>


extern "C" {
static void fmpz_poly_read_raw(fmpz_poly_t f, const void * data, size_t datasize) {
    FILE * pseudofile = fmemopen( (void *)data, datasize, "r");
    slong degree;
    fread((void*)&degree, sizeof(degree), 1, pseudofile);
    fmpz_t x;
    fmpz_init(x);
    fmpz_poly_zero(f);
    for(int k = degree; k >= 0; k--) {
        fmpz_inp_raw(x, pseudofile);
        fmpz_poly_set_coeff_fmpz(f, k, x);
    }
    fmpz_clear(x);
    fclose(pseudofile);
}

int polydb_get_entry(sqlite3 * db, fmpz_poly_t f, int * an, int ** orbit, int * orbitsize, int level, int weight, int chi, int j) {
    const char * query_sql = "SELECT an, orbit, polynomial from heckepolys WHERE "
                             "level=? and weight=? and chi=? and j=?";

    sqlite3_stmt * stmt;
    int result = sqlite3_prepare_v2(db, query_sql, -1, &stmt, NULL);
    if(result != SQLITE_OK) {
        fprintf(stderr, "ohno %d\n", result);
        exit(1);
    }
    sqlite3_bind_int(stmt, 1, level);
    sqlite3_bind_int(stmt, 2, weight);
    sqlite3_bind_int(stmt, 3, chi);
    sqlite3_bind_int(stmt, 4, j);

    result = sqlite3_step(stmt);
    if(result != SQLITE_ROW) {
        //fprintf(stderr, "ohno %d\n", result);
        return 0;
    }

    size_t polysize = sqlite3_column_bytes(stmt, 2);
    const void * polydata = sqlite3_column_blob(stmt, 2);

    fmpz_poly_read_raw(f, polydata, polysize);

    *an = sqlite3_column_int(stmt, 0);

    *orbitsize = sqlite3_column_bytes(stmt, 1)/sizeof(int);
    *orbit = (int*)malloc(*orbitsize * sizeof(int));
    memcpy( *orbit, sqlite3_column_blob(stmt, 1), *orbitsize * sizeof(int) );

    sqlite3_finalize(stmt);

    return 1;

}

int polydb_get_entries(sqlite3 * db,
                        fmpz_poly_t ** polys,
                        int * npolys,
                        int * an,
                        int ** chiorbit,
                        int * chiorbitsize,
                        int level,
                        int weight,
                        int chi) {
    const char * count_sql = "SELECT count(*) from heckepolys WHERE "
                             "level=? and weight=? and chi=?";

    const char * query_sql = "SELECT an, orbit, polynomial, j from heckepolys WHERE "
                             "level=? and weight=? and chi=?";

    sqlite3_stmt * stmt;
    int result = sqlite3_prepare_v2(db, count_sql, -1, &stmt, NULL);
    if(result != SQLITE_OK) {
        fprintf(stderr, "ohno %d\n", result);
        exit(1);
    }
    sqlite3_bind_int(stmt, 1, level);
    sqlite3_bind_int(stmt, 2, weight);
    sqlite3_bind_int(stmt, 3, chi);

    result = sqlite3_step(stmt);
    if(result != SQLITE_ROW) {
        return 0;
    }

    int count = sqlite3_column_int(stmt, 0);

    *npolys = count;
    sqlite3_finalize(stmt);

    *polys = (fmpz_poly_t*)malloc(count * sizeof(fmpz_poly_t));
    for(int k = 0; k < count; k++) {
        fmpz_poly_init((*polys)[k]);
    }

    if(count != 0) {
        sqlite3_prepare_v2(db, query_sql, -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, level);
        sqlite3_bind_int(stmt, 2, weight);
        sqlite3_bind_int(stmt, 3, chi);

        result = sqlite3_step(stmt);
        bool first = true;
        while(result == SQLITE_ROW) {
            if(first) {
                *an = sqlite3_column_int(stmt, 0);
                *chiorbitsize = sqlite3_column_bytes(stmt, 1)/sizeof(int);
                *chiorbit = (int*)malloc(*chiorbitsize * sizeof(int));
                memcpy( *chiorbit, sqlite3_column_blob(stmt, 1), *chiorbitsize * sizeof(int) );
                first = false;
            }
            int j = sqlite3_column_int(stmt, 3);

            size_t polysize = sqlite3_column_bytes(stmt, 2);
            const void * polydata = sqlite3_column_blob(stmt, 2);

            fmpz_poly_read_raw((*polys)[j-1], polydata, polysize);
            result = sqlite3_step(stmt);
        }
        sqlite3_finalize(stmt);
    }



    return count;

}

void polydb_init(sqlite3 * db) {
    const char sql1[] = "CREATE TABLE heckepolys (level INTEGER, weight INTEGER, chi INTEGER, j INTEGER,"
                        "orbitnumber INTEGER, operator BLOB, degree INTEGER, chiorbit BLOB, polynomial BLOB);";
    const char sql2[] = "CREATE INDEX degree_index on heckepolys (degree);";
    const char sql3[] = "CREATE INDEX level_weight_chi_j on heckepolys (level, weight, chi, j);";

    sqlite3_exec(db, sql1, NULL, NULL, NULL);
    sqlite3_exec(db, sql2, NULL, NULL, NULL);
    sqlite3_exec(db, sql3, NULL, NULL, NULL);
}

}
