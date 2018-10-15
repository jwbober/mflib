#include "cuspforms_weight1_modp.h"

#include "slint.h"

cuspforms_weight1_modp::cuspforms_weight1_modp(DirichletCharacter &_chi, long _p, int _verbose) {
    verbose = _verbose;
    chi = _chi.m;
    level = _chi.parent->q;
    long primitive_index;
    conductor = _chi.conductor(&primitive_index);
    p = _p;

    if(_chi.is_even()) {
        _dimension = 0;
        zero_for_trivial_reasons = 1;
        return;
    }

    q3 = level; if(q3 % 3 != 0) q3 *= 3;
    q4 = level; if(q4 % 2 != 0) q4 *= 2; if(q4 % 4 != 0) q4 *= 2;

    DirichletGroup G3(q3);
    DirichletGroup G4(q4);

    chi3 = G3.character(G3.index_from_primitive_character(3, 2)
           * G3.index_from_primitive_character(conductor, primitive_index) % q3);
    chi4 = G4.character(G4.index_from_primitive_character(4, 3)
           * G4.index_from_primitive_character(conductor, primitive_index) % q4);

    S3 = get_cuspforms_modp(chi3, 2, p);
    S4 = get_cuspforms_modp(chi4, 2, p);

    //int order = order_mod(_chi.m, level);

    // We're not going to touch p here. We have to assume that
    // it was computed correctly by the caller (which should probably be
    // the get_cuspforms_weight1_modp() function.)

    //if(p == 0) {
    //    p = order1 + 1;
    //    while(p < 10000) p += order1;
    //}
    //else {
    //    if( (p - 1) % order1 != 0) {
    //        p += order1 - (p - 1) % order1;
    //    }
    //}
    //while(!is_prime(p)) {p += order1;}

    nmod_init(&modp, p);

    phi_table = new int[level + 1];
    gcd_tables = new int*[level + 1];

    for(int k = 1; k < level + 1; k++) {
        if(level % k == 0) {
            gcd_tables[k] = new int[k];
            for(int j = 0; j < k; j++) {
                gcd_tables[k][j] = GCD(k,j);
            }
            phi_table[k] = euler_phi(k);
        }
    }

    chi_values = new long[level + 1];
    divisor_counts = new int[level + 1];
    _chi.values_mod_p(p, chi_values);

    DirichletGroup G(conductor);
    DirichletCharacter chip = G.character(primitive_index);
    chip_values = new long[conductor + 1];
    chip.values_mod_p(p, chip_values);

    divisors_of_level = divisors(level);
    for(int M : divisors_of_level) {
        divisor_counts[M] = ndivisors(M);
        if(M % conductor == 0 && M != level) {
            sublevels.push_back(M);
            DirichletGroup G(M);
            DirichletCharacter psi = G.character(G.index_from_primitive_character(conductor, primitive_index));
            cuspforms_weight1_modp * S2 = get_cuspforms_weight1_modp(psi, p, verbose);
            subspaces.push_back(S2);
        }
    }
}

cuspforms_weight1_modp::~cuspforms_weight1_modp() {
    if(zero_for_trivial_reasons) return;
    for(int k = 1; k < level + 1; k++) {
        if(level % k == 0) {
            delete [] gcd_tables[k];
        }
    }
    delete [] gcd_tables;
    delete [] phi_table;
    delete [] chi_values;
    delete [] divisor_counts;
    delete [] chip_values;
}
