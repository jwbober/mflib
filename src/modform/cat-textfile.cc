#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

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
    int level = atoi(argv[1]);
    int weight = atoi(argv[2]);
    int chi_number = atoi(argv[3]);

    string filename = "mf/" + to_string(level) + "/" + to_string(weight) + "/" +
                      to_string(level) + "." +
                      to_string(weight) + "." +
                      to_string(chi_number) + ".0";

    acb_t z;
    acb_init(z);
    acb_ptr coeffs = _acb_vec_init(lines_expected);
    int lines = parse_textfile(coeffs, filename, lines_expected);
    for(int k = 0; k < lines; k++) {
        //acb_printd(coeffs + k, 10);
        //cout << endl;
        acb_add(z, z, coeffs + k, 1000);
    }
    acb_printd(z, 10);
    _acb_vec_clear(coeffs, lines_expected);
    return 0;
}
