#include <iostream>
#include "mfformat.h"

#include "acb.h"

using namespace std;

int count = 0;

int callback(struct mfheader * header, int coeff_datasize, const void * coeff_data, acb_ptr coeffs) {
    cout << count << " " << header->level << " " << header->weight << " " << header->chi << " " << header->j << " ";
    acb_printd(coeffs+18, 10);
    cout << endl;
    count++;
    return 0;
}

int main(int argc, char ** argv) {
    iterate_through_sqlitefile(argv[1], callback, 1);
    return 0;
}
