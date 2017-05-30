#include "cuspforms_acb.h"

void cuspforms_acb::trace(acb_t out, int n) {
    if(n > 1 && new_dimension() == 0) { acb_set_ui(out, 0); return; }
    if(n < traces_computed) {acb_set(out, traces[n]); return; }

    int new_end = std::max( (int)(n * 1.1), n + 10 );
    compute_traces(new_end);
    acb_set(out, traces[n]);
}

void cuspforms_acb::trace_TnTm(acb_t out, int n, int m) {
    //
    // Return the trace mod p of TmTn acting on S_k(chi).
    //
    //TODO
    acb_set_ui(out, 0);
    if(n*m == 0) {return; }

    if(traces_computed < n*m + 1) {
         compute_traces(n*m + 1);
    }
    int g = GCD(m,n);

    acb_t z; acb_init(z);
    acb_set_ui(z, 0);
    for(int d = 1; d <= g; d++) {
        if(g % d != 0) continue; // this is silly, but this will usually
                                 // be a short sum
        if(GCD(d,level) != 1) continue;

        //long dk = PowerMod(d, k-1, p);
        //double dk = pow(d, weight - 1);
        acb_set_si(z, d);
        arb_pow_ui(acb_realref(z), acb_realref(z), weight - 1, prec);
        //acb_pow_ui(z, z, weight - 1, prec);
        acb_mul(z, z, chi_values[d % level], prec);
        acb_mul(z, z, traces[m/d * n/d], prec);
        acb_add(out, out, z, prec);
        //acb_addmul(out, traces[m/d * n/d], z, prec);
        //complex<double> z = trace(m/d * n/d) * dk * chi_values[d % level];
        //TnTm += z;
    }
    acb_clear(z);
}

void cuspforms_acb::trace_TpTnTm(acb_t out, int p, int n, int m) {
    //
    // Return the trace mod p of TmTn acting on S_k(chi).
    //
    //TODO
    acb_set_ui(out, 0);
    if(n*m == 0) {return; }

    if(traces_computed < n*m + 1) {
         compute_traces(n*m + 1);
    }
    int g = GCD(m,n);

    acb_t z; acb_init(z);
    acb_t z2; acb_init(z2);
    acb_set_ui(z, 0);
    for(int d = 1; d <= g; d++) {
        if(g % d != 0) continue; // this is silly, but this will usually
                                 // be a short sum
        if(GCD(d,level) != 1) continue;

        //long dk = PowerMod(d, k-1, p);
        //double dk = pow(d, weight - 1);
        acb_set_si(z, d);
        arb_pow_ui(acb_realref(z), acb_realref(z), weight - 1, prec);
        //acb_pow_ui(z, z, weight - 1, prec);
        acb_mul(z, z, chi_values[d % level], prec);
        trace_TnTm(z2, m/d * n/d, p);
        //acb_mul(z, z, traces[m/d * n/d], prec);
        acb_mul(z, z, z2, prec);
        acb_add(out, out, z, prec);
        //acb_addmul(out, traces[m/d * n/d], z, prec);
        //complex<double> z = trace(m/d * n/d) * dk * chi_values[d % level];
        //TnTm += z;
    }
    acb_clear(z);
    acb_clear(z2);
}
