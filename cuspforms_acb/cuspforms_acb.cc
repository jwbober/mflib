#include "cuspforms_acb.h"

static long psi(int N) {
    int_factorization_t f;
    factor(N, f);

    long PSI = 1;
    for(int k = 0; k < f.nfactors; k++) {
        PSI *= f.factors[k].f/f.factors[k].p * (long)(f.factors[k].p + 1);
    }
    return PSI;
}

cuspforms_acb::cuspforms_acb(DirichletCharacter &_chi, const int w, int _verbose) {
    verbose = _verbose;
    chi = _chi.m;
    level = _chi.parent->q;
    prec = _chi.parent->prec;
    long primitive_index;
    conductor = _chi.conductor(&primitive_index);
    weight = w;

    psi_table = new int[level + 1];
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

    modp_space = get_cuspforms_modp(_chi, weight, 10000000l, verbose);

    chi_values = new acb_t[level + 1];
    for(int k = 0; k <= level; k++) {
        acb_init(chi_values[k]);
        _chi.value(chi_values[k], k);
    }
    divisor_counts = new int[level + 1];

    DirichletGroup G(conductor, prec);
    DirichletCharacter chip = G.character(primitive_index);
    chip_values = new acb_t[conductor + 1];
    for(int k = 0; k <= conductor; k++) {
        acb_init(chip_values[k]);
        chip.value(chip_values[k], k);
    }

    divisors_of_level = divisors(level);
    for(int M : divisors_of_level) {
        divisor_counts[M] = ndivisors(M);
        if(M % conductor == 0 && M != level) {
            sublevels.push_back(M);
            DirichletGroup G(M, prec);
            DirichletCharacter psi = G.character(G.index_from_primitive_character(conductor, primitive_index));
            cuspforms_acb * S2 = get_cuspforms_acb(psi, weight, verbose);
            subspaces.push_back(S2);
        }
    }
}

cuspforms_acb::~cuspforms_acb() {
    delete [] psi_table;
    for(int k = 1; k < level + 1; k++) {
        if(level % k == 0) {
            delete [] gcd_tables[k];
        }
    }
    delete [] gcd_tables;
    delete [] phi_table;
    for(int k = 0; k <= level; k++) {
        acb_clear(chi_values[k]);
    }
    delete [] chi_values;
    delete [] divisor_counts;
    for(int k = 0; k <= conductor; k++) {
        acb_clear(chip_values[k]);
    }
    delete [] chip_values;

    if(traces_size > 0) {
        for(int k = 0; k < traces_size; k++) {
            acb_clear(traces[k]);
        }
        delete [] traces;
    }
}
