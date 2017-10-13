#include <iostream>

#include "slint.h"
#include "classnumbers.h"
#include "cuspforms_modp.h"

using std::cerr;
using std::endl;


static void square_divisors_mod03(int D, int * divisors, int& k) {
    // For f0, f1, ..., f_{j-1} such that
    //
    //  (f_i)^2 divides D AND D/(f_i)^2 == 0 or 3 mod 4
    //
    // we set divisors[k + i] = f_i, and then set k := k + j. If one of these
    // f_i is equal to 1, then we choose f_0 = 1, so divisors[k] = 1 after the
    // function exits. These f are otherwise unsorted.
    //
    // (We assume, of course, that divisors[k] through divisors[k + j - 1]
    // are valid memory locations...)
    //
    // What's going on here is that -D is (probably) a discriminant, and we
    // are finding all the f such that -D/f^2 is also a discriminant.
    //
    // For reference: from a separate check we know that an array of size
    // ceil(0.83 * X + 10) will be large enough to hold all the f for D == 0 or 3 mod 4,
    // D < x (That's a bit generous, even, but I don't want to bother being too precise.)

    int_factorization_t fac;
    factor(D, fac);
    int e[fac.nfactors];
    for(int j = 0; j < fac.nfactors; j++) {e[j] = 0;}
    int f = 1;                                                      // We'll write f = prod(fac.factors[k].p^(e[k]/2))
    do {
        if(D/(f*f) % 4 == 0 || D/(f*f) % 4 == 3) {
            divisors[k] = f;
            k++;
        }
        int j = 0;                                                  // This bit is a bit messy.
        while(j < fac.nfactors && e[j] + 2 > fac.factors[j].e) {    //
            while(e[j] > 0) {                                       // We're iterating through the f such that f^2
                e[j] -= 2;                                          // divides n by using the known factorization
                f = f/fac.factors[j].p;                             // of it. This means iterating through all the
            }                                                       // tuples (e_0, e_1, ..., e_k) where e_j is even
            j++;                                                    // and e_i <= fac.factors[j].e.
        }                                                           //
        if(j == fac.nfactors) f = 0;
        else {
            e[j] += 2;
            f *= fac.factors[j].p;
        }
    } while(f != 0);
}

void cuspforms_modp::compute_traces(int end) {
    if(end < traces.size()) return;
    for(cuspforms_modp *subspace : subspaces) {
        subspace->compute_traces(end);
    }
    int start = traces.size();
    traces.resize(end, 0);

    long one_over_twelve = nmod_inv(12, modp);
    long one_over_three = nmod_inv(3, modp);
    long one_over_two = nmod_inv(2, modp);

    long A1 = nmod_mul(psi_table[level], one_over_twelve, modp);
    A1 = nmod_mul(A1, weight - 1, modp);
    int z = (int)sqrt(start);                               //
    while(z*z < start) z++;                                 // A1 is only nonzero when n is a square,
    for( ; z*z < end; z++) {                                // and in the weight 2 case it hardly
        //traces[z*z] += A1 * chi_values[z % level] % p;    // depends on n at all.
        long a = z % p;
        a = nmod_pow_ui(a, weight - 2, modp);
        //a = PowerMod(a, weight - 2, p);
        a = nmod_mul(a, A1, modp);
        NMOD_ADDMUL(traces[z*z], a, chi_values[z % level], modp); // traces[z*z] += A1 * chi(z) mod p
    }
    // computation of A2

    long print_interval = 0;
    if(verbose > 1) {
        cerr << "cuspforms_modp: level " << level << endl;
        cerr << "cuspforms_modp: start = " << start << " end = " << end << endl;
    }
    if(start == end) return;
    if(verbose > 1) {
        cerr << "cuspforms_modp: A2:"; print_interval = round(2 * sqrt(4*end)/70); print_interval = std::max(print_interval, 1l);
    }

    int * square_divisors = new int[ (int)(4*end*.83 + 11) ];   // To avoid repeated factorizations in the following loop,
    int * square_divisors_indices = new int[4*end + 1];         // we create some arrays such that for each D == 0 or 3 mod 4,
    square_divisors_indices[3] = 0;                             // D < 4*end,
    square_divisors[0] = 1;                                     //
    for(int k = 1, j = 1; j < end; j++) {                       // square_divisors[square_divisors_indices[D]]
        square_divisors_indices[4*j] = k;                       //
        square_divisors_mod03(4*j, square_divisors, k);         // is the beginning of a 1-terminated list of f such that
        square_divisors_indices[4*j + 3] = k;                   // f^2 divides D and D/f^2 == 0 or 3 mod 4.
        square_divisors_mod03(4*j + 3, square_divisors, k);
        square_divisors[k] = 1;     // This is about to get
                                    // overwritten, except for
                                    // the last iteration.
    }

    int t = sqrt(4*end);
    if(t*t == 4*end) t--;
    for(t = -t; (long)t*t < 4l*end; t++) {
        if(verbose > 1 && t % print_interval == 0)  {
            cerr << '.';
            cerr.flush();
        }
        for(long x = 0; x < level; x++) {
            long chi = chi_values[x];
            if(chi == 0) continue;
            // These values of x and t will contribute to those n for
            // which the level divides x^2 - tx + n.
            int starting_n = std::max(start, t*t/4 + 1);
            if((x*x - t*x + starting_n) % level != 0) {
                starting_n += (level - (x*x - t*x + starting_n) % level);
            }
            for(long n = starting_n; n < end; n += level) {
                // Now for this (t,n,x) we will add a term to traces[n]
                // for each f such that f*f divides (t^2 - 4n) AND
                // gcd(f, N) * N divides x*x - t*x + n
                //
                // So we find the square divisors of t^2 - 4n that we've
                // already found
                long D = 4*n - t*t;
                int k = square_divisors_indices[D];
                int f = square_divisors[k];
                long polyterm = evalpoly(t, n);
                do {
                    int g;
                    if(f == 1) {
                        g = 1;
                    }
                    else if(f < level) {
                        g = gcd_tables[level][f];
                    }
                    else {
                        g = gcd_tables[level][f % level];
                    }
                    if(g == 1 || ((long)x*x - t*x + n) % (g * level) == 0) {
                        long z = psi_table[level]/psi_table[level/g];
                        if(z >= p) z %= p;
                        z = nmod_mul(z, chi, modp);
                        z = nmod_mul(z, classnumbers[D/(f*f)], modp);
                        if(D/(f*f) == 3) z = nmod_mul(z, one_over_three, modp); // (z * one_over_three) % p;
                        if(D/(f*f) == 4) z = nmod_mul(z, one_over_two, modp);   //(z * one_over_two) % p;
                        z = nmod_mul(z, one_over_two, modp);
                        z = nmod_mul(z, polyterm, modp);
                        traces[n] = nmod_sub(traces[n], z, modp);
                    }
                    k++;
                    f = square_divisors[k];
                } while(f != 1);
            }
        }
    }
    delete [] square_divisors;
    delete [] square_divisors_indices;
    if(verbose > 1) { cerr << endl; cerr.flush(); }

    if(verbose > 1) { cerr << "cuspforms_modp: A3:"; print_interval = (long)std::max(round(sqrt(end)/70), 1.0); }
    for(int d = 1; d*d < end; d++) {
        if(verbose > 1 && d % print_interval == 0) {
            cerr << '.';
            cerr.flush();
        }
        int starting_n = start;
        if(d*d > starting_n) {
            starting_n = d*d;
        }
        else if(starting_n % d != 0) {
            starting_n += (d - starting_n % d);
        }                                              // In A3, we'll get a contribution
        for(int n = starting_n; n < end; n += d) {     // to traces[n] for each d that divides n
            long a = 0;                                // as long as d^2 <= n
            for(int c : divisors_of_level) {
                int z = (n/d - d) % (level/conductor);
                if(z < 0) z += (level/conductor);
                int g = gcd_tables[c][level/c % c];
                if(gcd_tables[level/conductor][z] % g == 0) {
                    long y;
                    if(g == 1) {
                        y = CRT(d, n/d, c, level/c);
                        a = nmod_add(a, chi_values[y], modp);
                    }
                    else {
                        unsigned long g1, u, v;
                        if(c >= level/c)
                            g1 = n_xgcd(&u, &v, c, level/c);
                        else {
                            g1 = n_xgcd(&v, &u, level/c, c);
                            u %= level/c;                               // We are doing a CRT lift
                            u = level/c - u;                            // here, but the moduli might
                        }
                        y = (d + c * u * (n/d - d)/g1) % (level/g);     // not be coprime. Hence the
                        if(y < 0) y += (level/g);                       // complication.

                                                                    // XXX: I've got this computation wrong
                        //a = a + phi_table[g] * chi_values[y];     // a few times. I hope it is right now.
                        NMOD_ADDMUL(a, phi_table[g], chi_values[y], modp);
                        //a %= p;

                    }
                }
            }
            if(d*d == n) a = nmod_mul(a, one_over_two, modp);
            //long dpow = PowerMod(d, weight - 1, p);
            long dpow = nmod_pow_ui(d, weight - 1, modp);
            a = nmod_mul(a, dpow, modp);
            traces[n] = nmod_sub(traces[n], a, modp);       // traces -= d*a mod p
        }
    }
    if(verbose > 1) { cerr << endl; cerr.flush(); }
    // computation of A4
    if(chi == 1 && weight == 2) {
        for(int t = 1; t < end; t++) {
            int starting_n = start;
            if(starting_n % t != 0) {
                starting_n += (t - starting_n % t);
            }
            for(int n = starting_n; n < end; n += t) {
                if(gcd_tables[level][n/t % level] == 1) {
                    traces[n] = nmod_add(traces[n], t, modp);
                    //traces[n] += t;
                    //traces[n] %= p;
                }
            }
        }
    }

    // before sieving, record the dimension of the full space,
    // provided that the dimension hasn't already been set.
    if(_dimension == -1)
        _dimension = traces[1];
    // now we sieve.

    int q = conductor;
    for(int k = 0; k < sublevels.size(); k++) {
        cuspforms_modp * subspace = subspaces[k];
        int N = level;
        int M = sublevels[k];
        for(int n = start; n < end; n++) {
            if(chi_values[n % level] != 0) {     // This is just a check on the GCD.
                long z = nmod_mul(divisor_counts[N/M], subspace->traces[n], modp);
                traces[n] = nmod_sub(traces[n], z, modp);
            }
            else {
                // I'm sure that we can this more efficiently, but let's
                // try the slow way first. This is not the slow part of
                // the program anyway...
                int x = N/M; while(GCD(x, n) > 1) {x = x/GCD(x,n);}     // now x = (N/M) / GCD(N/M, n^oo)
                int y = N; while(GCD(y, M) > 1) {y = y/GCD(y,M);}       // now y = N / GCD(N, M^oo)
                for(int b : divisors(y)) {
                    if(b == 1 && N == M) continue;
                    if(n % (b*b) == 0) {
                        int mu = mobius(b);
                        //long bpow = PowerMod(b, weight - 1, p);
                        long bpow = nmod_pow_ui(b, weight - 1, modp);
                        if(mu == 1) {
                            long z = nmod_mul(divisor_counts[x], chip_values[b % q], modp);
                            z = nmod_mul(z, bpow, modp);
                            z = nmod_mul(z, subspace->traces[n/(b*b)], modp);
                            traces[n] = nmod_sub(traces[n], z, modp);
                        }
                        else if(mu == -1) {
                            long z = nmod_mul(divisor_counts[x], chip_values[b % q], modp);
                            z = nmod_mul(z, bpow, modp);
                            z = nmod_mul(z, subspace->traces[n/(b*b)], modp);
                            traces[n] = nmod_add(traces[n], z, modp);
                        }
                    }
                }
            }
        }
    }

}
