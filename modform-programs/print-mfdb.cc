#include <iostream>
#include <complex>

#include "mfformat.h"

#include "acb.h"

using namespace std;

int count = 0;

static complex<double> acb_get_z(const acb_t in) {
    double x = arf_get_d(arb_midref(acb_realref(in)), ARF_RND_NEAR);
    double y = arf_get_d(arb_midref(acb_imagref(in)), ARF_RND_NEAR);
    return complex<double>(x,y);
}

unsigned int print_start = 1;
unsigned int print_end = 10;

int callback(struct mfheader * header, int coeff_datasize, const void * coeff_data, acb_ptr coeffs) {
    cout << header->level << "." << header->weight << "." << header->chi << "." << header->j
        << " " << header->ncoeffs << " " << header->prec;

    unsigned int start = max(0u, print_start - 1);
    unsigned int end = min(print_end, header->ncoeffs);
    for(int k = start; k < end; k++) {
        cout << " " << acb_get_z(coeffs + k);
    }
    cout << endl;
    count++;
    return 0;
}

int main(int argc, char ** argv) {

    if(argc < 2) {
        cout << "Usage:" << endl
             << "./print-mfdb filename [start] [end]" << endl;
        return 0;
    }

    if(argc > 2) {
        print_start = atoi(argv[2]);
    }
    if(argc > 3) {
        print_end = atoi(argv[3]);
    }

    iterate_through_sqlitefile(argv[1], callback, 1);
    return 0;
}
