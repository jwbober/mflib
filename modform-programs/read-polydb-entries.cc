#include <iostream>
#include <string>

#include "mfformat.h"

#include "flint/fmpz_poly.h"

using namespace std;

int main(int argc, char ** argv) {
    if(argc != 5) {
        cout << "usage: " << argv[0] << " dbname level weight chi" << endl;
        return 0;
    }
    sqlite3 * db;
    sqlite3_open(argv[1], &db);
    int level = atoi(argv[2]);
    int weight = atoi(argv[3]);
    int chi = atoi(argv[4]);

    fmpz_poly_t * polys;
    int npolys;
    int an;
    int * orbit;
    int orbitsize;
    int result = polydb_get_entries(db, &polys, &npolys, &an, &orbit, &orbitsize, level, weight, chi);
    if(!result)
        return 0;

    cout << an;
    for(int k = 0; k < orbitsize; k++) {
        cout << " " << orbit[k];
    }
    cout << endl;
    for(int k = 0; k < npolys; k++) {
        fmpz_poly_print_pretty(polys[k], "x");
        cout << endl;
    }

    for(int k = 0; k < npolys; k++) {
        fmpz_poly_clear(polys[k]);
    }
    free(polys);
    free(orbit);
    sqlite3_close(db);
    return 0;
}
