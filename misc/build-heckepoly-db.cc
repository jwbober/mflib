#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstring>


#include "flint/fmpz_poly.h"

#include "arb_poly.h"
#include "acb.h"

#include "sqlite3.h"

using namespace std;

size_t prec = 5000;

std::vector<std::string> split(const std::string &text, char sep) {
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }
  tokens.push_back(text.substr(start));
  return tokens;
}

void bober_fmpz_poly_set_str_pretty(fmpz_poly_t poly, const char * str, size_t length) {
    FILE * pseudofile = fmemopen( (void * )str, length + 1, "r");
    char * x;
    fmpz_poly_fread_pretty(pseudofile, poly, &x);
    fclose(pseudofile);
    flint_free(x);
}

void fmpz_poly_write_raw( char ** data, size_t * datasize, fmpz_poly_t f) {
    FILE * pseudofile = open_memstream(data, datasize);

    slong degree = fmpz_poly_degree(f);
    fwrite((void*)&degree, sizeof(degree), 1, pseudofile);
    fmpz_t x;
    fmpz_init(x);
    for(int k = degree; k >= 0; k--) {
        fmpz_poly_get_coeff_fmpz(x, f, k);
        fmpz_out_raw(pseudofile, x);
    }
    fmpz_clear(x);
    fclose(pseudofile);
}

void arb_poly_gcd(arb_poly_t out, arb_poly_t f1, arb_poly_t f2, int expected_degree, size_t prec) {
    arb_poly_t _g1;
    arb_poly_t _g2;

    arb_poly_struct * g1 = _g1;
    arb_poly_struct * g2 = _g2;

    arb_poly_t g3;

    arb_poly_init(g1);
    arb_poly_init(g2);
    arb_poly_init(g3);

    arb_poly_set(g1, f1);
    arb_poly_set(g2, f2);

    arb_t x;
    arb_t y;
    arb_init(x);
    arb_init(y);

    int d1 = arb_poly_degree(g1);
    int d2 = arb_poly_degree(g2);
    while(d1 > 0 && arb_contains_zero(arb_poly_get_coeff_ptr(g1, d1))) {
        d1--;
    }
    while(d2 > 0 && arb_contains_zero(arb_poly_get_coeff_ptr(g2, d2))) {
        d2--;
    }
    if(d1 < d2) {
        arb_poly_struct * g;
        g = g1;
        g1 = g2;
        g2 = g;

        int d = d1;
        d1 = d2;
        d2 = d;
    }
    while(d2 > expected_degree) {
        arb_set(x, arb_poly_get_coeff_ptr(g1, d1));
        arb_set(y, arb_poly_get_coeff_ptr(g2, d2));

        arb_poly_scalar_div(g1, g1, x, prec);
        arb_poly_scalar_div(g2, g2, y, prec);
        arb_poly_set_coeff_si(g1, d1, 1);
        arb_poly_set_coeff_si(g2, d2, 1);
        //arb_poly_scalar_mul(g1, g1, y, prec);
        //arb_poly_scalar_mul(g3, g2, x, prec);
        arb_poly_shift_left(g3, g2, d1 - d2);
        arb_poly_sub(g1, g1, g3, prec);

        d1 = arb_poly_degree(g1);

        while(d1 > 0 && arb_contains_zero(arb_poly_get_coeff_ptr(g1, d1))) {
            d1--;
        }
        while(d2 > 0 && arb_contains_zero(arb_poly_get_coeff_ptr(g2, d2))) {
            d2--;
        }

        if(d1 <= d2) {
            arb_poly_struct * g;
            g = g1;
            g1 = g2;
            g2 = g;

            int d = d1;
            d1 = d2;
            d2 = d;
        }

        //cout << d1 << " " << d2 << endl;
    }

    arb_poly_set(out, g2);
    arb_poly_clear(g1);
    arb_poly_clear(g2);
    arb_poly_clear(g3);
}

int main(int argc, char ** argv) {
    if(argc < 3) {
        cout << "usage: ./build-heckepoly-db [polyfilelist] [dbfilename]" << endl;
        return 0;
    }
    ifstream filelist(argv[1]);
    string filename;
    filelist >> filename;
    sqlite3 * db;
    sqlite3_open(argv[2], &db);

    const char * table_sql =  "CREATE TABLE heckepolys ("
                            "level  INTEGER,"
                            "weight INTEGER,"
                            "chi    INTEGER,"
                            "j      INTEGER,"
                            "an     INTEGER,"
                            "orbit  BLOB,"
                            "polynomial BLOB);";

    const char * insert_sql = "INSERT INTO heckepolys VALUES (?, ?, ?, ?, ?, ?, ?);";

    sqlite3_exec(db, table_sql, NULL, 0, NULL);
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, 0, NULL);

    while(!filelist.eof()) {
        ifstream heckefile(filename);
        heckefile.peek();
        if(!heckefile.eof()) {
            cout << filename << endl;
            string line;
            getline(heckefile, line);
            vector<string> firstline = split(line, ' ');
            vector<int> characters;

            int level = stoi(firstline[0]);

            int weight = stoi(firstline[1]);
            for(auto x = firstline.begin() + 2; x != firstline.end(); x++) {
                characters.push_back(stoi(*x));
            }

            //for(auto x : characters) {
            //    cout << x << endl;
            //}

            getline(heckefile, line);
            int an = atoi(line.c_str());
            getline(heckefile, line);
            vector<string> degree_strings = split(line, ' ');

            int npolys = degree_strings.size();

            vector<string> polystrings;

            for(int k = 0; k < degree_strings.size(); k++) {
                getline(heckefile, line);
                fmpz_poly_t poly;
                fmpz_poly_init(poly);
                bober_fmpz_poly_set_str_pretty(poly, line.c_str(), line.size());

                char * polydata;
                size_t polydatasize;
                fmpz_poly_write_raw(&polydata, &polydatasize, poly);

                sqlite3_stmt * insert_stmt;
                sqlite3_prepare_v2(db, insert_sql, -1, &insert_stmt, NULL);
                sqlite3_bind_int(insert_stmt, 1, level);
                sqlite3_bind_int(insert_stmt, 2, weight);
                sqlite3_bind_int(insert_stmt, 3, characters[0]);
                sqlite3_bind_int(insert_stmt, 4, k + 1);
                sqlite3_bind_int(insert_stmt, 5, an);
                sqlite3_bind_blob(insert_stmt, 6, (void *)characters.data(), sizeof(int)*characters.size(), SQLITE_TRANSIENT);
                sqlite3_bind_blob(insert_stmt, 7, (void *)polydata, polydatasize, free);
                sqlite3_step(insert_stmt);
                sqlite3_finalize(insert_stmt);
            }
        }
        filelist >> filename;
        heckefile.close();
    }

    sqlite3_exec(db, "END TRANSACTION", NULL, 0, NULL);

    return 0;
}
