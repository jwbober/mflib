#include <iostream>
#include <string>

#include "mfformat.h"

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
    struct mfheader header;
    acb_ptr coeffs;

    int result = mfdb_get_entry(db, &header, &coeffs, level, weight, chi, j);
    if(!result)
        return 0;

    for(int k = 0; k < header.ncoeffs; k++) {
        acb_printd(coeffs + k, 10);
        cout << endl;
    }
    sqlite3_close(db);
    return 0;
}
