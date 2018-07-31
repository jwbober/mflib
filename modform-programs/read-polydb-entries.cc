#include <iostream>
#include <string>

#include "mfformat.h"

#include "flint/fmpz_poly.h"

using namespace std;

int main(int argc, char ** argv) {
    if(argc != 4) {
        cout << "usage: " << argv[0] << " dbname level weight" << endl;
        return 0;
    }
    sqlite3 * db;
    sqlite3_open(argv[1], &db);
    int level = atoi(argv[2]);
    int weight = atoi(argv[3]);

    fmpz_poly_t * polys;
    int npolys;
    int ** mforbits;
    size_t * mforbitsizes;
    int * whatevernumbers;

    int result = polydb_get_entries(db, &polys, &npolys, &mforbits, &mforbitsizes, &whatevernumbers, level, weight);
    if(!result)
        return 0;

    cout  << "here" << endl;

    for(int k = 0; k < npolys; k++) {
        cout << whatevernumbers[k] << " ";
        fmpz_poly_print_pretty(polys[k], "x");
        for(int j = 0; j < mforbitsizes[k]/2; j++) {
            cout << " " << mforbits[k][2*j] << "." << mforbits[k][2*j+1];
        }
        cout << endl;
    }
    // some things to free. can't be bothered right now.
    sqlite3_close(db);
    return 0;
}
