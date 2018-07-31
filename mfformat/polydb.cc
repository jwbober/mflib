#include "mfformat.h"
#include "flint-extras.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
using namespace std;


extern "C" {


//    const char * table_sql =  "CREATE TABLE heckepolys ("
//                            "level          INTEGER,"
//                            "weight         INTEGER,"
//                            "chi            INTEGER,"
//                            "whatevernumber INTEGER,"
//                            "labelnumber    INTEGER,"
//                            "operator       BLOB,"
//                            "degree         INTEGER,"
//                            "mforbit        BLOB,"
//                            "polynomial     BLOB);";

int polydb_get_entry(sqlite3 * db, fmpz_poly_t f, int ** mforbit, size_t * orbitsize, int level, int weight, int chi, int j) {
    const char * query_sql = "SELECT mforbit, polynomial from heckepolys WHERE "
                             "level=? and weight=? and chi=? and whatevernumber=?";

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

    size_t polysize = sqlite3_column_bytes(stmt, 1);
    const void * polydata = sqlite3_column_blob(stmt, 1);
    x_fmpz_poly_read_raw(f, polydata, polysize);

    *orbitsize = sqlite3_column_bytes(stmt, 0)/sizeof(int);
    *mforbit = (int *)malloc(*orbitsize*sizeof(int));
    memcpy(*mforbit, sqlite3_column_blob(stmt, 0), *orbitsize*sizeof(int));

    sqlite3_finalize(stmt);

    return 1;

}

int polydb_get_entries(sqlite3 * db,
                        fmpz_poly_t ** polys,
                        int * npolys,
                        int *** mforbits,
                        size_t ** mforbitsizes,
                        int ** whatevernumbers,
                        int level,
                        int weight) {
    const char * count_sql = "SELECT count(*) from heckepolys WHERE "
                             "level=? and weight=?";

    const char * query_sql = "SELECT mforbit, polynomial, whatevernumber from heckepolys WHERE "
                             "level=? and weight=?";

    sqlite3_stmt * stmt;
    int result = sqlite3_prepare_v2(db, count_sql, -1, &stmt, NULL);
    if(result != SQLITE_OK) {
        fprintf(stderr, "ohno %d\n", result);
        exit(1);
    }
    sqlite3_bind_int(stmt, 1, level);
    sqlite3_bind_int(stmt, 2, weight);

    result = sqlite3_step(stmt);
    if(result != SQLITE_ROW) {
        return 0;
    }

    int count = sqlite3_column_int(stmt, 0);

    *npolys = count;
    sqlite3_finalize(stmt);

    *polys = (fmpz_poly_t*)malloc(count * sizeof(fmpz_poly_t));
    *whatevernumbers = (int *)malloc(count * sizeof(int));
    *mforbits = (int**)malloc(count*sizeof(int*));
    *mforbitsizes = (size_t*)malloc(count*sizeof(int));

    for(int k = 0; k < count; k++) {
        fmpz_poly_init((*polys)[k]);
    }

    cout << count << endl;
    if(count != 0) {
        sqlite3_prepare_v2(db, query_sql, -1, &stmt, NULL);
        sqlite3_bind_int(stmt, 1, level);
        sqlite3_bind_int(stmt, 2, weight);

        result = sqlite3_step(stmt);
        int k = 0;
        while(result == SQLITE_ROW) {
            cout << k << endl;
            (*mforbitsizes)[k] = (size_t)sqlite3_column_bytes(stmt, 0)/sizeof(int);
            (*mforbits)[k] = (int*)malloc((*mforbitsizes)[k]*sizeof(int));
            memcpy((*mforbits)[k], sqlite3_column_blob(stmt, 0), (*mforbitsizes)[k]*sizeof(int));

            size_t polysize = sqlite3_column_bytes(stmt, 1);
            const void * polydata = sqlite3_column_blob(stmt, 1);
            x_fmpz_poly_read_raw((*polys)[k], polydata, polysize);

            (*whatevernumbers)[k] = sqlite3_column_int(stmt, 2);

            result = sqlite3_step(stmt);
            k++;
        }
        sqlite3_finalize(stmt);
    }



    return count;

}

int polydb_insert(  sqlite3 * db,
                    fmpz_poly_t f,
                    int * hecke_operator,
                    size_t hecke_operator_size,
                    int * mforbit,
                    size_t orbitsize,
                    int level,
                    int weight,
                    int chi,
                    int whatevernumber,
                    int labelnumber) {
    const char * insert_sql = "INSERT INTO heckepolys VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);";
    int degree = fmpz_poly_degree(f);

    char * polydata;
    size_t polydatasize;
    x_fmpz_poly_write_raw(&polydata, &polydatasize, f);

    sqlite3_stmt * insert_stmt;
    sqlite3_prepare_v2(db, insert_sql, -1, &insert_stmt, NULL);
    sqlite3_bind_int(insert_stmt, 1, level);
    sqlite3_bind_int(insert_stmt, 2, weight);
    sqlite3_bind_int(insert_stmt, 3, chi);
    sqlite3_bind_int(insert_stmt, 4, whatevernumber);
    sqlite3_bind_int(insert_stmt, 5, labelnumber);
    sqlite3_bind_blob(insert_stmt, 6, (void *)hecke_operator, sizeof(int)*hecke_operator_size, SQLITE_STATIC);
    sqlite3_bind_int(insert_stmt, 7, degree);
    sqlite3_bind_blob(insert_stmt, 8, (void *)mforbit, sizeof(int)*orbitsize, SQLITE_TRANSIENT);
    sqlite3_bind_blob(insert_stmt, 9, (void *)polydata, polydatasize, free);
    sqlite3_step(insert_stmt);
    return sqlite3_finalize(insert_stmt);
}

void polydb_init(sqlite3 * db) {
    const char * table_sql =  "CREATE TABLE heckepolys ("
                            "level          INTEGER,"
                            "weight         INTEGER,"
                            "chi            INTEGER,"
                            "whatevernumber INTEGER,"
                            "labelnumber    INTEGER,"
                            "operator       BLOB,"
                            "degree         INTEGER,"
                            "mforbit        BLOB,"
                            "polynomial     BLOB);";

    const char sql2[] = "CREATE INDEX degree_index on heckepolys (degree);";
    const char sql3[] = "CREATE INDEX level_weight_chi_j on heckepolys (level, weight, chi, j);";

    sqlite3_exec(db, table_sql, NULL, NULL, NULL);
    sqlite3_exec(db, sql2, NULL, NULL, NULL);
    sqlite3_exec(db, sql3, NULL, NULL, NULL);
}

}
