#include <iostream>
#include <string>

#include "mfformat.h"

#include "flint/fmpz_poly.h"

using namespace std;

int main(int argc, char ** argv) {
    if(argc < 6) {
        cout << "usage: " << argv[0] << " dbname level weight chi j" << endl;
        return 0;
    }
    sqlite3 * db;
    sqlite3_open(argv[1], &db);
    int level = atoi(argv[2]);
    int weight = atoi(argv[3]);
    int chi = atoi(argv[4]);
    int j = atoi(argv[5]);

    fmpz_poly_t f;
    fmpz_poly_init(f);
    int an;
    int * orbit;
    int orbitsize;
    int result = polydb_get_entry(db, f, &an, &orbit, &orbitsize, level, weight, chi, j);
    if(!result)
        return 0;

    cout << an;
    for(int k = 0; k < orbitsize; k++) {
        cout << " " << orbit[k];
    }
    cout << " ";
    fmpz_poly_print_pretty(f, "x");
    cout << endl;
    free(orbit);
    sqlite3_close(db);
    return 0;
}
