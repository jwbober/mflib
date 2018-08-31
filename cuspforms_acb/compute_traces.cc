#include <iostream>
#include <thread>
#include <cstring>

#include "classnumbers.h"
#include "cuspforms_acb.h"

using std::cerr;
using std::cout;
using std::endl;
using std::memset;

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

void cuspforms_acb::compute_traces(int end) {
    if(end < traces_computed) return;
    for(cuspforms_acb *subspace : subspaces) {
        subspace->compute_traces(end);
    }
    int start = traces_computed;
    if(verbose) cout << "cuspforms_acb: level " << level << endl;
    if(verbose) cout << "cuspforms_acb: start = " << start << " end = " << end << endl;
    if(start == end) return;
    if(end > traces_size) {
        if(verbose > 1) cout << "resizing" << endl;
        acb_t * new_traces = new acb_t[end];
        for(int k = 0; k < traces_size; k++) {
            new_traces[k][0] = traces[k][0];
        }
        for(int k = traces_size; k < end; k++) {
            acb_init(new_traces[k]);
        }
        if(traces_size > 0)
            delete [] traces;
        traces = new_traces;
        traces_size = end;
        if(verbose > 1) cout << "resized" << endl;
    }

    if(end > 2) {
        if(new_dimension() == 0) {
            for(int k = start; k < end; k++) {
                acb_zero(traces[k]);
            }
            return;
        }
    }

    arb_t t1; arb_init(t1);
    acb_t s1; acb_init(s1);
    acb_t s2; acb_init(s2);
    fmpz_t x1; fmpz_init(x1);

    arb_t A1; arb_init(A1);
    arb_set_si(A1, psi_table[level]);
    arb_div_ui(A1, A1, 12, prec);
    arb_mul_si(A1, A1, weight - 1, prec);
    //complex<double> A1 = psi_table[level]/12.0 * (weight - 1.0);
    int z = (int)sqrt(start);
    while(z*z < start) z++;

    arb_t npow;
    arb_init(npow);
    for( ; z*z < end; z++) {
        arb_set_si(npow, z);
        arb_pow_ui(npow, npow, weight - 2, prec);
        acb_mul_arb(s1, chi_values[z % level], npow, prec);  // s1 = chi(z) * n^{weight - 2}
        acb_mul_arb(s1, s1, A1, prec);                       // s1 = A1 * n^(weight - 2) * chi(z)
        acb_add(traces[z*z], traces[z*z], s1, prec);         // traces[z*z] += s1
        //traces[z*z] += npow * A1 * chi_values[z % level];
    }
    arb_clear(npow);

    long print_interval = 0;
    if(verbose > 0) { cout << "A2:"; print_interval = round(2 * sqrt(4*end)/70)/nthreads; print_interval = std::max(print_interval, 1l); }

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

    int tstart = sqrt(4*end);
    if(tstart*tstart == 4*end) tstart--;
    tstart = -tstart;

    acb_ptr thread_contribution[nthreads];
    if(nthreads == 1) {
        thread_contribution[0] = &traces[start][0];
    }
    else {
        if(verbose) cerr << "initializing storage for " << nthreads << " threads." << endl;
        for(int j = 0; j < nthreads; j++) {
            //cerr << j << " " << (end - start) * sizeof(acb_struct) << endl;
            //thread_contribution[j] = _acb_vec_init(end - start);
            thread_contribution[j] = new acb_struct[end - start];
            memset(thread_contribution[j], 0, (end - start) * sizeof(acb_struct));
        }
    }

    auto A2thread_function = [start, end, tstart, square_divisors, square_divisors_indices, this, print_interval, &thread_contribution](int threadid) {
        arb_t t1; arb_init(t1);
        acb_t s1; acb_init(s1);
        fmpz_t polyterm;
        fmpz_init(polyterm);
        for(int t = tstart + threadid; (long)t*t < 4l*end; t += nthreads) {
            if(threadid == 0 && verbose > 0 && t % print_interval == 0)  {
                cout << '.';
                cout.flush();
            }

            for(long x = 0; x < level; x++) {
                if(GCD(x, level) != 1) continue;
                // These values of x and t will contribute to those n for
                // which the level divides x^2 - tx + n.
                int starting_n = std::max(start, t*t/4 + 1);
                if((x*x - t*x + starting_n) % level != 0) {
                    starting_n += (level - (x*x - t*x + starting_n) % level);
                }
                for(int n = starting_n; n < end; n += level) {
                    // Now for this (t,n,x) we will add a term to traces[n]
                    // for each f such that f*f divides (t^2 - 4n) AND
                    // gcd(f, N) * N divides x*x - t*x + n
                    //
                    // So we find the square divisors of t^2 - 4n that we've
                    // already found
                    int D = 4*n - t*t;
                    int k = square_divisors_indices[D];
                    int f = square_divisors[k];
                    fmpz_t tt = {t};
                    fmpz_t nn = {n};
                    evalpoly(polyterm, tt, nn);
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
                            arb_set_si(t1, psi_table[level]);
                            arb_div_si(t1, t1, psi_table[level/g], prec);
                            arb_mul_si(t1, t1, classnumbers[D/(f*f)], prec);

                            if(D/(f*f) == 3) arb_div_si(t1, t1, 3, prec);
                            if(D/(f*f) == 4) arb_div_si(t1, t1, 2, prec);
                            arb_div_si(t1, t1, 2, prec);
                            arb_mul_fmpz(t1, t1, polyterm, prec);
                            acb_mul_arb(s1, chi_values[x], t1, prec);
                            //acb_sub(traces[n], traces[n], s1, prec);
                            acb_sub(thread_contribution[threadid] + n - start, thread_contribution[threadid] + n - start, s1, prec);
                        }
                        k++;
                        f = square_divisors[k];
                    } while(f != 1);
                }
            }
        }
        arb_clear(t1);
        acb_clear(s1);
        fmpz_clear(polyterm);
    };

    if(verbose) cerr << "launching " << nthreads << " threads." << endl;
    std::thread threads[nthreads];
    for(int j = 0; j < nthreads; j++) {
        threads[j] = std::thread(A2thread_function, j);
    }
    for(int j = 0; j < nthreads; j++) {
        threads[j].join();
    }
    if(verbose) cerr << "Threads finished." << endl;
    if(nthreads != 1) {
        for(int n = start; n < end; n++) {
            for(int j = 0; j < nthreads; j++) {
                acb_add(traces[n], traces[n], thread_contribution[j] + n - start, prec);
            }
        }
        for(int j = 0; j < nthreads; j++) {
            for(int k = 0; k < end - start; k++) {
                acb_clear(thread_contribution[j] + k);
            }
            //_acb_vec_clear(thread_contribution[j], end - start);
            delete [] thread_contribution[j];
        }
    }

    delete [] square_divisors;
    delete [] square_divisors_indices;
    if(verbose > 0) { cout << endl; cout.flush(); }

    //cout << traces[1] << endl;

    if(verbose > 0) { cout << "A3:"; print_interval = (long)std::max(round(sqrt(end)/70), 1.0); }
    for(int d = 1; d*d < end; d++) {
        if(verbose > 0 && d % print_interval == 0) {
            cout << '.';
            cout.flush();
        }
        int starting_n = start;
        if(d*d > starting_n) {
            starting_n = d*d;
        }
        else if(starting_n % d != 0) {
            starting_n += (d - starting_n % d);
        }                                              // In A3, we'll get a contribution
        for(int n = starting_n; n < end; n += d) {     // to traces[n] for each d that divides n
            //complex<double> a = 0;                     // as long as d^2 <= n
            acb_set_ui(s1, 0);
            for(int c : divisors_of_level) {
                int z = (n/d - d) % (level/conductor);
                if(z < 0) z += (level/conductor);
                int g = gcd_tables[c][level/c % c];
                if(gcd_tables[level/conductor][z] % g == 0) {
                    long y;
                    if(g == 1) {
                        y = CRT(d, n/d, c, level/c);
                        acb_add(s1, s1, chi_values[y], prec);
                        //a += chi_values[y];
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
                        //a = a + phi_table[g] * chi_values[y];       // a few times. I hope it is right now.

                        acb_mul_si(s2, chi_values[y], phi_table[g], prec);
                        acb_add(s1, s1, s2, prec);
                        //a += (double)phi_table[g] * chi_values[y];
                        //NMOD_ADDMUL(a, phi_table[g], chi_values[y], modp);
                        //a %= p;

                        //cout << "y = " << y << " " << d << " " << n/d << " " << c << " " << level/c << endl;
                    }
                }
            }
            if(d*d == n) acb_div_ui(s1, s1, 2, prec); //a /= 2.0; //a = nmod_mul(a, one_over_two, modp);
            //a *= (double)d; //a = nmod_mul(a, d, modp);

            fmpz_set_si(x1, d);
            fmpz_pow_ui(x1, x1, weight - 1);
            acb_mul_fmpz(s1, s1, x1, prec);
            acb_sub(traces[n], traces[n], s1, prec);

            //a *= pow(d, weight - 1);
            //traces[n] -= a;
            //traces[n] = nmod_sub(traces[n], a, modp);       // traces -= d*a mod p
        }
    }
    if(verbose > 0) { cout << endl; cout.flush(); }
    //cout << traces[1] << endl;

    //cout << "computing A4 for level " << level << endl;
    // computation of A4
    if(chi == 1 && weight == 2) {
        for(int t = 1; t < end; t++) {
            int starting_n = start;
            if(starting_n % t != 0) {
                starting_n += (t - starting_n % t);
            }
            for(int n = starting_n; n < end; n += t) {
                if(gcd_tables[level][n/t % level] == 1) {
                    acb_add_si(traces[n], traces[n], t, prec);
                    //traces[n] += (double)t;
                    //traces[n] = nmod_add(traces[n], t, modp);
                    //traces[n] += t;
                    //traces[n] %= p;
                }
            }
        }
    }

    // now we sieve.

    int q = conductor;
    for(unsigned int k = 0; k < sublevels.size(); k++) {
        cuspforms_acb * subspace = subspaces[k];
        int N = level;
        int M = sublevels[k];
        for(int n = start; n < end; n++) {
            if(GCD(n, level) == 1) {
            //if(chi_values[n % level] != 0.0) {     // This is just a check on the GCD.
                acb_submul_ui(traces[n], subspace->traces[n], divisor_counts[N/M], prec);
                //traces[n] -= (complex<double>)divisor_counts[N/M] * subspace->traces[n];
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
                        arb_set_si(t1, b);
                        arb_pow_ui(t1, t1, weight - 1, prec);
                        arb_mul_si(t1, t1, mu * divisor_counts[x], prec);
                        acb_mul_arb(s1, subspace->traces[n/(b*b)], t1, prec);
                        acb_submul(traces[n], chip_values[b % q], s1, prec);
                        //traces[n] -= mu * (double)divisor_counts[x] * chip_values[b % q] * pow((double)b, weight - 1) * subspace->traces[n/(b*b)];
                    }
                }
            }
        }
    }
    traces_computed = end;
}
