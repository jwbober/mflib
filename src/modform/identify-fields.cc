#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#define USE_ARB
#include "characters.h"

#include "acb.h"
#include "acb_poly.h"

#include "flint/fmpz_poly.h"

using namespace std;

const int max_dimension = 2000;
acb_t roots[max_dimension];

void build_polynomial(acb_poly_t out, int start, int end) {
    if(end == start) {
        acb_poly_set_coeff_si(out, 1, 1);
        acb_poly_set_coeff_acb(out, 0, roots[start]);
        return;
    }
    acb_poly_t a;
    acb_poly_init(a);
    int mid = (start + end)/2;
    build_polynomial(out, start, mid);
    build_polynomial(a, mid + 1, end);
    acb_poly_mul(out, out, a, 2000);
    //for(int k = start; k < end; k++) {
    //    acb_poly_set_coeff_acb(a, 0, roots[k]);
    //    acb_poly_mul(out, out, a, 4000);
    //}
    acb_poly_clear(a);
}

void acb_set_zzzzzzzz(acb_t out, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d,
                                 fmpz_t e, fmpz_t f, fmpz_t g, fmpz_t h) {
    arf_set_fmpz_2exp(arb_midref(acb_realref(out)), a, b);
    mag_set_fmpz_2exp_fmpz(arb_radref(acb_realref(out)), c, d);
    arf_set_fmpz_2exp(arb_midref(acb_imagref(out)), e, f);
    mag_set_fmpz_2exp_fmpz(arb_radref(acb_imagref(out)), g, h);
}

int main(int argc, char ** argv) {
    int level = atoi(argv[1]);
    int weight = atoi(argv[2]);
    int chi_number = atoi(argv[3]);

    for(int k = 0; k < max_dimension; k++) {
        acb_init(roots[k]);
    }
    string datapath = "mf/" + to_string(level) + "/" + to_string(weight) + "/";
    acb_t a2;
    acb_init(a2);
    fmpz_t a, b, c, d, e, f, g, h;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(d);
    fmpz_init(e);
    fmpz_init(f);
    fmpz_init(g);
    fmpz_init(h);



    DirichletGroup G(level);
    for(auto orbit : G.galois_orbits()) {
        int k = 0;
        for(auto m : orbit) {
            string datafile_base =  to_string(level)
                                    + "." + to_string(weight)
                                    + "." + to_string(m) + '.';
            int n = 0;
            string filename = datapath + datafile_base + to_string(n);
            ifstream infile(filename);
            if(!infile) {
                break;
            }
            while(infile) {
                string line;
                getline(infile, line);

                string x;
                infile >> x;
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

                // cout << filename << endl;
                acb_set_zzzzzzzz(roots[k], a, b, c, d, e, f, g, h);
                acb_neg(roots[k], roots[k]);
                infile.close();
                n++;
                k++;
                filename = datapath + datafile_base + to_string(n);
                infile.open(filename);
            }
        }

        int nroots = k;

        if(nroots == 0) continue;

        acb_poly_t charpoly;
        acb_poly_init(charpoly);
        build_polynomial(charpoly, 0, nroots - 1);

        fmpz_poly_t charpolyz;
        fmpz_poly_init(charpolyz);
        int result = acb_poly_get_unique_fmpz_poly(charpolyz, charpoly);
        cout << level << " " << weight << " " << *(orbit.begin()) << " ";
        if(result) {
            // cout << "factoring" << endl;
            fmpz_poly_factor_t factors;
            fmpz_poly_factor_init(factors);
            //fmpz_poly_print(charpolyz);
            fmpz_poly_factor_zassenhaus(factors, charpolyz);
            for(int k = 0; k < factors[0].num; k++) {
                fmpz_poly_print_pretty(&factors[0].p[k], "x");
                if(k < factors[0].num - 1) cout << " ";
            }
            //fmpz_poly_factor_print(factors);
            cout << endl;
        }
        else {
            cout << "error" << endl;
            acb_poly_printd(charpoly, 10);
            cout << endl;
        }
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(d);
    fmpz_clear(e);
    fmpz_clear(f);
    fmpz_clear(g);
    fmpz_clear(h);
    acb_clear(a2);
    return 0;
}
