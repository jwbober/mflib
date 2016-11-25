#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <complex>

#include "acb.h"

using namespace std;

void acb_set_zzzzzzzz(acb_t out, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d,
                                 fmpz_t e, fmpz_t f, fmpz_t g, fmpz_t h) {
    arf_set_fmpz_2exp(arb_midref(acb_realref(out)), a, b);
    mag_set_fmpz_2exp_fmpz(arb_radref(acb_realref(out)), c, d);
    arf_set_fmpz_2exp(arb_midref(acb_imagref(out)), e, f);
    mag_set_fmpz_2exp_fmpz(arb_radref(acb_imagref(out)), g, h);
}

int parse_textfile(acb_ptr out, string filename, int maxlines) {
    ifstream infile(filename);
    if(!infile) return 0;

    fmpz_t a, b, c, d, e, f, g, h;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(d);
    fmpz_init(e);
    fmpz_init(f);
    fmpz_init(g);
    fmpz_init(h);

    int count = 0;
    while(infile && count < maxlines) {
        string x;
        infile >> x;
        if(!infile) break;
        fmpz_set_str(a, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(b, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(c, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(d, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(e, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(f, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(g, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(h, x.c_str(), 10);

        acb_set_zzzzzzzz(out + count, a, b, c, d, e, f, g, h);
        count++;
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(d);
    fmpz_clear(e);
    fmpz_clear(f);
    fmpz_clear(g);
    fmpz_clear(h);
    return count;
}

const int lines_expected = 4000;

int main(int argc, char ** argv) {
    if(argc < 2) {
        cout << "./print-mftextfile filename [prec]" << endl;
        return 0;
    }
    string filename(argv[1]);
    int prec = 10;
    if(argc > 2) {
        prec = atoi(argv[2]);
    }

    acb_t z;
    acb_init(z);
    acb_ptr coeffs = _acb_vec_init(lines_expected);
    int lines = parse_textfile(coeffs, filename, lines_expected);
    for(int k = 0; k < lines; k++) {
        //double x = arf_get_d(arb_midref(acb_realref(coeffs + k)), ARF_RND_NEAR);
        //double y = arf_get_d(arb_midref(acb_imagref(coeffs + k)), ARF_RND_NEAR);
        //if(abs(x) < 1e-20) x = 0;
        //if(abs(y) < 1e-20) y = 0;
        //complex<double> z(x,y);
        cout << k + 1 << "\t\t";
        acb_printd(coeffs + k, prec);
        cout  << endl;
        //cout << endl;
        //acb_add(z, z, coeffs + k, 1000);
    }
    //acb_printd(z, 10);
    _acb_vec_clear(coeffs, lines_expected);
    return 0;
}
