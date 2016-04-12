#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#define USE_ARB
#include "characters.h"

#include "acb.h"
#include "acb_poly.h"

#include "flint/fmpz_poly.h"
#include "flint/NTL-interface.h"

#include "NTL/ZZX.h"
#include "NTL/ZZXFactoring.h"

using namespace std;

const int max_dimension = 2000;
acb_t roots[max_dimension];

void build_polynomial(acb_poly_t out, int start, int end) {
    if(end == start) {
        acb_poly_zero(out);
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
    acb_t z;
    acb_init(z);

    acb_poly_t charpoly;
    acb_poly_init(charpoly);

    fmpz_poly_t charpolyz;
    fmpz_poly_init(charpolyz);

    fmpz_poly_t hh;
    fmpz_poly_init(hh);

    for(auto orbit : G.galois_orbits()) {
        for(int k = 0; k < max_dimension; k++) {
            acb_zero(roots[k]);
        }
        int skip = 1;
        bool finished = false;
        while(!finished) {
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
                    for(int l = 0; l < skip; l++) getline(infile, line);

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

                    acb_set_zzzzzzzz(z, a, b, c, d, e, f, g, h);
                    acb_sub(roots[k], roots[k], z, 2000);
                    infile.close();
                    n++;
                    k++;
                    filename = datapath + datafile_base + to_string(n);
                    infile.open(filename);
                }
            }

            int nroots = k;

            if(nroots == 0) {
                finished = true;
                continue;
            }

            build_polynomial(charpoly, 0, nroots - 1);
            int result = acb_poly_get_unique_fmpz_poly(charpolyz, charpoly);
            if(!result) {
                for(auto m : orbit) {
                    cout << level << "." << weight << "." << m << " ?" << endl;
                }
                finished = true;
            }
            else {
                if(fmpz_poly_is_squarefree(charpolyz)) {
                    NTL::ZZ content;
                    NTL::ZZX ff;
                    NTL::vec_pair_ZZX_long factors;
                    fmpz_poly_get_ZZX(ff, charpolyz);
                    NTL::factor(content, factors, ff);
                    if(factors.length() == 1) {
                        for(auto m : orbit) {
                            cout << level << "." << weight << "." << m << " ";
                            fmpz_poly_print_pretty(charpolyz, "x");
                            cout << endl;
                        }
                    }
                    else {
                        // now we need to match up the eigenvalues with the polynomials.

                        // negate all of the roots so that they are actually roots...
                        for(int l = 0; l < nroots; l++) {
                            acb_neg(roots[l], roots[l]);
                        }
                        acb_poly_t g;
                        acb_poly_init(g);
                        int root_to_polynomial[nroots];
                        for(int l = 0; l < nroots; l++) {
                            root_to_polynomial[l] = -1;
                        }
                        for(int l = 0; l < factors.length(); l++) {
                            fmpz_poly_set_ZZX(hh, factors[l].a);
                            acb_poly_set_fmpz_poly(g, hh, 2000);
                            int nroots_found = 0;
                            for(int k = 0; k < nroots; k++) {
                                acb_poly_evaluate(z, g, roots[k], 2000);
                                if(acb_contains_zero(z)) {
                                    nroots_found++;
                                    if(root_to_polynomial[k] != -1) {
                                        // something is the root of two polynomials, so something
                                        // went wrong. print question marks, break out of the loop,
                                        // and continue with the next orbit.
                                        for(auto m : orbit) {
                                            cout << level << "." << weight << "." << m << " ??" << endl;
                                        }
                                        k = nroots;
                                        l = factors.length();
                                        finished = true;
                                    }
                                    else {
                                        root_to_polynomial[k] = l;
                                    }
                                }
                            }
                            if(nroots_found != acb_poly_degree(g)) {
                                for(auto m : orbit) {
                                    cout << level << "." << weight << "." << m << " ???" << endl;
                                }
                                l = factors.length();
                                finished = true;
                            }
                        }
                        int single_orbit_dimension = nroots/orbit.size();
                        int k = 0;
                        for(auto m : orbit) {
                            for(int l = 0; l < single_orbit_dimension; l++) {
                                cout << level << "." << weight << "." << m << "." << l << " ";
                                fmpz_poly_set_ZZX(hh, factors[root_to_polynomial[k]].a);
                                fmpz_poly_print_pretty(hh, "x");
                                cout << endl;
                                k++;
                            }
                        }

                        //cout << level << "." << weight << "." << *(orbit.begin()) << " ";
                        //for(int k = 0; k < factors[0].num; k++) {
                        //    fmpz_poly_print_pretty(&factors[0].p[k], "x");
                        //    if(k < factors[0].num - 1) cout << " ";
                        //}
                        //cout << endl;
                        acb_poly_clear(g);
                    }
                    finished = true;
                }
                else {
                    skip++;
                }
            }
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
