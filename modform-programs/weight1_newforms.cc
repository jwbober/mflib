#include <iostream>
#include <set>
#include <map>

#include "cuspforms_weight1_modp.h"
#include "classnumbers.h"

using namespace std;

int main(int argc, char ** argv) {
    init_classnumbers();
    load_factor_table();

    if(argc < 4) {
        cout << "Usage: " << argv[0] << " level chi ncoeffs" << endl;
        return 0;
    }

    int level = atoi(argv[1]);
    int chinumber = atoi(argv[2]);
    int ncoeffs = atoi(argv[3]);

    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chinumber);
    int chiorder = order_mod(chinumber, level);
    int h = LCM(2*chiorder, 120);
    long p = 10000;
    p = 10000 - (p % h) + 1;
    while(!is_prime(p)) {
        p = p + 2*h;
    }

    long * chivalues = new long[level];
    chi.values_mod_p(p, chivalues);

    set<long> _chivalueset;
    for(int k = 0; k < level; k++) {
        _chivalueset.insert(chivalues[k]);
    }
    vector<long> chivalueset;
    for(long z : _chivalueset) {
        chivalueset.push_back(z);
    }
    map<long, map<long,long>> satake_tables;

    long a = primitive_root(p);
    long zeta_h = PowerMod(a, (p - 1)/h, p);
    long zeta_h_inv = InvMod(zeta_h, p);
    long zeta = 1;
    long zeta_inv = 1;

    long total_count = 0;
    long repeat_count = 0;
    for(int k = 0; k < h; k++) {
        for(long chivalue : chivalueset) {
            total_count++;
            long z = (zeta + chivalue * zeta_inv) % p;
            if(satake_tables[chivalue].find(z) != satake_tables[chivalue].end()) {
                repeat_count++;
            }
            else
                satake_tables[chivalue][z] = k;
        }
        zeta = zeta * zeta_h % p;
        zeta_inv = zeta_inv * zeta_h_inv % p;
    }
    //cout << chivalueset.size() << endl;
    //cout << total_count << endl;
    //cout << repeat_count << endl;

    cuspforms_weight1_modp * S = get_cuspforms_weight1_modp(chi, p);
    if(p != S->p) cout << "ohno" << endl;

    nmod_mat_t newforms;

    int newform_count = S->newspace_basis(newforms, ncoeffs);
    int d = S->new_dimension();

    cout << p << endl;
    for(int k = 0; k < newform_count; k++) {
        for(int q = 2; q < ncoeffs; q = next_prime(q)) {
            long z = nmod_mat_entry(newforms, k, q);
            if(level % q != 0) {
                int alpha_numerator = satake_tables[chivalues[q % level]][z];
                int d = GCD(alpha_numerator, h);
                cout << " " << alpha_numerator/d << "/" << h/d;
            }
            else {
                if(z == 0) cout << " []";
                else {
                    int alpha_numerator = h/order_mod(z, p);
                    int d = GCD(alpha_numerator, h);
                    cout << " [" << alpha_numerator/d << "/" << h/d << "]";
                }
            }
        }

        //if(k < newform_count) cout << "* ";
        //cout << nmod_mat_entry(newforms, d - k - 1, 0);
        //for(int j = 1; j < ncoeffs; j++) {
        //    cout << " " << nmod_mat_entry(newforms, k, j);
        //}
        cout << endl;
    }

    return 0;
}
