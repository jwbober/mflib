#include "cuspforms_modp.h"

long cuspforms_modp::trace(int n) {
    if(n > 1 && new_dimension() == 0) return 0;
    if(n < traces.size()) return traces[n];

    int new_end = std::max( (int)(n * 1.1), n + 10 );
    compute_traces(new_end);
    return traces[n];
}

long cuspforms_modp::trace_TnTm(int n, int m) {
    //
    // Return the trace mod p of TmTn acting on S_k(chi).
    //

    if(m*n == 0) return 0;
    int g = GCD(m,n);
    long TnTm = 0;

    for(int d = 1; d <= g; d++) {
        if(g % d != 0) continue; // this is silly, but this will usually
                                 // be a short sum
        if(GCD(d,level) != 1) continue;

        long dk = nmod_pow_ui(d, weight - 1, modp);
        long z = nmod_mul(trace(m/d * n/d), dk, modp);
        z = nmod_mul(z, chi_values[d % level], modp);
        TnTm = nmod_add(TnTm, z, modp);
    }
    return TnTm;
}

long cuspforms_modp::trace_TpTnTm(int p, int n, int m) {
    //
    // Return the trace mod p of TpTmTn acting on S_k(chi).
    //
    //TODO
    long trace = 0;
    if(n*m == 0) {return 0; }
     compute_traces(n*m*p + 1);
    int g = GCD(m,n);

    long z = 0;
    long z2 = 0;
    for(int d = 1; d <= g; d++) {
        if(g % d != 0) continue; // this is silly, but this will usually
                                 // be a short sum
        if(GCD(d,level) != 1) continue;

        z = nmod_pow_ui(d, weight - 1, modp);
        z = nmod_mul(z, chi_values[d % level], modp);
        z2 = trace_TnTm(m/d * n/d, p);
        z = nmod_mul(z, z2, modp);
        trace = nmod_add(trace, z, modp);
    }
    return trace;
}
