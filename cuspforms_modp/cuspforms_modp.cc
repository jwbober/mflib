#include "cuspforms_modp.h"

static long psi(int N) {
    int_factorization_t f;
    factor(N, f);

    long PSI = 1;
    for(int k = 0; k < f.nfactors; k++) {
        PSI *= f.factors[k].f/f.factors[k].p * (long)(f.factors[k].p + 1);
    }
    return PSI;
}

cuspforms_modp::cuspforms_modp(DirichletCharacter &_chi, const int w, long _p, int _verbose) {
    verbose = _verbose;
    chi = _chi.m;
    level = _chi.parent->q;
    long primitive_index;
    conductor = _chi.conductor(&primitive_index);
    weight = w;
    p = _p;

    int order = order_mod(_chi.m, level);
    if(p == 0) {
        p = order + 1;
        while(p < 10000) p += order;
    }
    else {
        if( (p - 1) % order != 0) {
            p += order - (p - 1) % order;
        }
    }
    while(!is_prime(p)) {p += order;}

    nmod_init(&modp, p);

    psi_table = new long[level + 1];
    phi_table = new int[level + 1];
    gcd_tables = new int*[level + 1];

    for(int k = 1; k < level + 1; k++) {
        if(level % k == 0) {
            psi_table[k] = psi(k);
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
            cuspforms_modp * S2 = get_cuspforms_modp(psi, weight, p, verbose);
            subspaces.push_back(S2);
        }
    }
}

cuspforms_modp::~cuspforms_modp() {
    delete [] psi_table;
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
