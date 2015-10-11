#include <iostream>
#include <vector>
#include <complex>

#include "flint/flint.h"
#include "flint/nmod_mat.h"

#include "slint.h"
#include "characters.h"
#include "classnumbers.h"

#include "trace-formula.h"


using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::vector;
using std::complex;


//
// Uncomment for range checking on nmod_mat_t accesses.
//
//#undef nmod_mat_entry
//
//static mp_limb_t& nmod_mat_entry(nmod_mat_t mat, int i, int j) {
//    if(i >= mat->r) { cout << "nmod_mat_entry() called with i out of range."; exit(0); }
//    if(j >= mat->c) { cout << "nmod_mat_entry() called with j out of range."; exit(0); }
//    return mat->rows[i][j];
//}

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

long trace_TmTn_mod_p(
        long * traces,
        long * chi_values,
        int k,
        int level,
        int m,
        int n,
        long p) {
    //
    // Return the trace mod p of TmTn acting on S_k(chi).
    //

    nmod_t modp;            // we really should pass an
    nmod_init(&modp, p);    // nmod_t to this function.
    int g = GCD(m,n);
    long TnTm = 0;

    for(int d = 1; d <= g; d++) {
        if(g % d != 0) continue; // this is silly, but this will usually
                                 // be a short sum
        if(GCD(d,level) != 1) continue;
        
        //long dk = PowerMod(d, k-1, p);
        long dk = nmod_pow_ui(d, k - 1, modp);
        long z = nmod_mul(traces[m/d * n/d], dk, modp);
        z = nmod_mul(z, chi_values[d % level], modp);
        TnTm = nmod_add(TnTm, z, modp);
        //TnTm += ((traces[m/d * n/d] * dk) % p) * chi_values[d % level] % p;
        //TnTm %= p;
    }
    return TnTm;
}




int psi(int N) {
    int_factorization_t f;
    factor(N, f);

    //int PSI = N;
    //for(int k = 0; k < f.nfactors; k++) {
    //    PSI += PSI/f.factors[k].p;
    //}
    int PSI = 1;
    for(int k = 0; k < f.nfactors; k++) {
        PSI *= f.factors[k].f/f.factors[k].p * (f.factors[k].p + 1);
    }
    return PSI;
}

void trace_Tn_unsieved_weight2(complex<double> * traces, int start, int end, int level, DirichletCharacter& chi, int verbose) {
    traces = traces - start; // We adjust the traces pointer so that we
                             // can work with traces[start] through
                             // trace[end-1] instead of subtracting start all
                             // the time.

    int ** gcd_tables = new int*[level + 1];
    int * psi_table = new int[level + 1];

    for(int k = 1; k < level + 1; k++) {
        if(level % k == 0) {
            psi_table[k] = psi(k);
            gcd_tables[k] = new int[k];
            for(int j = 0; j < k; j++) {
                gcd_tables[k][j] = GCD(k,j);
            }
        }
    }

    int conductor = chi.conductor();

    vector<int> divisors_of_level = divisors(level);

    // Start by zeroing out the range.
    for(int n = start; n < end; n++) {
        traces[n] = 0.0;
    }
    
    complex<double> A1 = (complex<double>)(psi_table[level]/12.0);
    int z = (int)sqrt(start);                               //
    while(z*z < start) z++;                                 // A1 is only nonzero when n is a square,
    for( ; z*z < end; z++) {                                // and in the weight 2 case it hardly
        traces[z*z] += A1 * chi.value(z);    // depends on n at all.
    }

    cout << traces[1] << endl;

    // computation of A2
    int t = sqrt(4*end);
    if(t*t == 4*end) t--;
    for(t = -t; (long)t*t < 4l*end; t++) {
        for(int x = 0; x < level; x++) {
            if(chi.value(x) == 0.0) continue;
            // These values of x and t will contribute to those n for
            // which the level divides x^2 - tx + n.
            int starting_n = max(start, t*t/4 + 1);
            if((x*x - t*x + starting_n) % level != 0) {
                starting_n += (level - (x*x - t*x + starting_n) % level);
            }
            for(int n = starting_n; n < end; n += level) {
                // Now for this (t,n,x) we will add a term to traces[n]
                // for each f such that f*f divides (t^2 - 4n) AND
                // gcd(f, N) * N divides x*x - t*x + n
                //
                // So we find the square divisors of t^2 - 4n...
                int_factorization_t fac;
                factor(4*n - t*t, fac);                                         //
                int D = 1;                                                      //
                for(int k = 0; k < fac.nfactors; k++) {                         // Set D = squarefree_part(4n - t^2)
                    if(fac.factors[k].e % 2 == 1) D *= fac.factors[k].p;        //
                }                                                               //
                if(D % 4 == 2 || D % 4 == 3) {                                  // and then set D to be the negative of the
                    D *= 4;                                                     // discriminant of Q(t^2 - 4n)
                }                                                               //
                D = 4*n - t*t;
                int e[fac.nfactors];
                for(int k = 0; k < fac.nfactors; k++) {e[k] = 0;}
                int f = 1;                                                      // We'll write f = prod(fac.factors[k].p^(e[k]/2))
                do {                                                            // (But we don't directly compute it that way.)
                    if( (4l*n - (long)t*t)/((long)f*f) % 4 == 0 || (4l*n - (long)t*t)/((long)f*f) % 4 == 3) {
                        int g = gcd_tables[level][f % level];
                        if( (x*x - t*x + n) % (g * level) == 0) {
                            complex<double> z = (complex<double>)(psi_table[level]/psi_table[level/g]) * chi.value(x) * (complex<double>)classnumbers[D/(f*f)];
                            if(D/(f*f) == 3) z /= 3.0;
                            if(D/(f*f) == 4) z /= 4.0;
                            traces[n] -= z/2.0;
                        }
                    }
                    int j = 0;                                                  // This bit is a bit messy.
                    while(j < fac.nfactors && e[j] + 2 > fac.factors[j].e) {    //
                        while(e[j] > 0) {                                       // We're iterating through the f such that f^2
                            e[j] -= 2;                                          // divides 4n - t^2 by using the known factorization
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
        }
    }
    cout << traces[1] << endl;

    // computation of A3
    for(int d = 1; d*d < end; d++) {            // In A3, we'll get a contribution
        for(int n = d*d; n < end; n += d) {     // to traces[n] for each d that divides n
            complex<double> a = 0.0;                          // as long as d^2 <= n
            for(int c : divisors_of_level) {
                int z = (n/d - d) % (level/conductor);
                if(z < 0) z += (level/conductor);
                int g = gcd_tables[c][level/c % c];
                if(gcd_tables[level/conductor][z] % g == 0) {
                    int y;
                    if(g == 1) {
                        y = CRT(d, n/d, c, level/c);
                        a = a + chi.value(y);
                    }
                    else {
                        // XXX THIS IS WRONG
                        if(GCD(c/g, level/c) == 1) {            // We are doing a CRT lift
                            y = CRT(d, n/d, c/g, level/c);      // here, but the moduli might
                        }                                       // not be coprime. Hence the
                        else { // GCD(c, (N/c)/g) == 1          // complication.
                            y = CRT(d, n/d, c/g, (level/c)/g);
                        }
                        a = a + (complex<double>)euler_phi(g) * chi.value(y);
                        if(n == 1) cout << "* y = " << y << " " << d << " " << n/d << " " << c << " " << level/c << endl;
                    }
                }
            }
            if(d*d == n) a /= 2.0;
            traces[n] -=  (complex<double>)d * a;
        }
    }
    cout << traces[1] << endl;
    
    // computation of A4
    if(chi.m == 1) {
        for(int t = 1; t < end; t++) {
            int starting_n = start;
            if(starting_n % t != 0) {
                starting_n += (t - starting_n % t);
            }
            for(int n = starting_n; n < end; n += t) {
                if(gcd_tables[level][n/t % level] == 1) {
                    traces[n] += t;
                }
            }
        }
    }
    cout << traces[1] << endl;

    delete [] psi_table;
    for(int k = 1; k < level + 1; k++) {
        if(level % k == 0) {
            delete [] gcd_tables[k];
        }
    }
    delete [] gcd_tables;
}

void trace_Tn_modp_unsieved_weight2(long * traces, int start, int end, int level, long p, long * chi_values, DirichletCharacter& chi, int verbose) {
    // We expect that p is not 2 or 3.

    traces = traces - start; // We adjust the traces pointer so that we
                             // can work with traces[start] through
                             // trace[end-1] instead of subtracting start all
                             // the time.

    nmod_t modp;
    nmod_init(&modp, p);
    long one_over_twelve = nmod_inv(12, modp);
    long one_over_three = nmod_inv(3, modp);
    long one_over_two = nmod_inv(2, modp);

    int ** gcd_tables = new int*[level + 1];
    int * psi_table = new int[level + 1];
    int * phi_table = new int[level + 1];

    for(int k = 1; k < level + 1; k++) {
        if(level % k == 0) {
            psi_table[k] = psi(k);
            gcd_tables[k] = new int[k];
            for(int j = 0; j < k; j++) {
                gcd_tables[k][j] = GCD(k,j);
            }
            phi_table[k] = euler_phi(k) % p;
        }
    }

    int conductor = chi.conductor();

    vector<int> divisors_of_level = divisors(level);

    // Start by zeroing out the range.
    for(int n = start; n < end; n++) {
        traces[n] = 0;
    }

    long A1 = nmod_mul(psi_table[level], one_over_twelve, modp);
    int z = (int)sqrt(start);                               //
    while(z*z < start) z++;                                 // A1 is only nonzero when n is a square,
    for( ; z*z < end; z++) {                                // and in the weight 2 case it hardly
        //traces[z*z] += A1 * chi_values[z % level] % p;    // depends on n at all.
        NMOD_ADDMUL(traces[z*z], A1, chi_values[z % level], modp); // traces[z*z] += A1 * chi(z) mod p
    }

    //cout << traces[1] << endl;
    // computation of A2
    //cout << "computing A2 for level " << level << endl;
    
    long print_interval = 0;
    if(verbose > 0) { cerr << "A2:"; print_interval = round(2 * sqrt(4*end)/70); print_interval = max(print_interval, 1l); }

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
        if(verbose > 0 && t % print_interval == 0)  {
            cerr << '.';
            cerr.flush();
        }
        for(int x = 0; x < level; x++) {
            long chi = chi_values[x];
            if(chi == 0) continue;
            // These values of x and t will contribute to those n for
            // which the level divides x^2 - tx + n.
            int starting_n = max(start, t*t/4 + 1);
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
    if(verbose > 0) { cerr << endl; cerr.flush(); }
    //cout << traces[1] << endl;

    if(verbose > 0) { cerr << "A3:"; print_interval = (long)max(round(sqrt(end)/70), 1.0); }
    for(int d = 1; d*d < end; d++) {
        if(verbose > 0 && d % print_interval == 0) {
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
                        unsigned long g1, u, v;                     // We are doing a CRT lift
                        g1 = n_xgcd(&u, &v, c, level/c);            // here, but the moduli might
                        y = (d + c * u * (n/d - d)/g1) % (level/g); // not be coprime. Hence the
                        if(y < 0) y += (level/g);                   // complication.

                                                                    // XXX: I've got this computation wrong
                        //a = a + phi_table[g] * chi_values[y];       // a few times. I hope it is right now.
                        NMOD_ADDMUL(a, phi_table[g], chi_values[y], modp);
                        //a %= p;

                        //cout << "y = " << y << " " << d << " " << n/d << " " << c << " " << level/c << endl;
                    }
                }
            }
            if(d*d == n) a = nmod_mul(a, one_over_two, modp);
            a = nmod_mul(a, d, modp);
            traces[n] = nmod_sub(traces[n], a, modp);       // traces -= d*a mod p
        }
    }
    if(verbose > 0) { cerr << endl; cerr.flush(); }
    //cout << traces[1] << endl;
    
    //cout << "computing A4 for level " << level << endl;
    // computation of A4
    if(chi.m == 1) {
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
    //cout << traces[1] << endl;

    //cout << "cleaning up for level " << level << endl;
    delete [] psi_table;
    for(int k = 1; k < level + 1; k++) {
        if(level % k == 0) {
            delete [] gcd_tables[k];
        }
    }
    delete [] gcd_tables;
    delete [] phi_table;
}

void sieve_trace_Tn_modp_on_weight2_for_newspaces(vector<long> * traces, int start, int end, int level, long p, long ** chi_values, DirichletCharacter& chi, int verbose) {
    //
    // On input, we should have:
    //  - level is the (maximum) level that we are after
    //  - traces[M][n] = Trace(Tn|S_2(M, chi)) (mod p) for every M > 1 such that q | M | level and start <= n < end;
    //    (note the special exception that we don't require that traces[1] is filled with zeros.)
    //
    //  - chi_values[M][n] = The value mod p at n of the (restriction or extension) of chi to M, for 0 <= n < M;
    //
    //  - Is chi mod q or mod level? I don't know yet. What am I going to use chi for, anyway, other than for
    //    the computation of the conductor? Maybe I should compute chi_values in this function...
    
    nmod_t modp;
    nmod_init(&modp, p);

    int q = chi.conductor();
    if(q == level) {// There is nothing to do!
        return;
    }
        
    int * divisor_counts = new int[level + 1];
    vector<int> level_divisors = divisors(level);
    vector<int> sublevels;                          // We fill sublevels with all of the M such that
    for(int M : level_divisors) {                   // q | M | level.
        if(M % q == 0) {
            sublevels.push_back(M);
        }
        divisor_counts[M] = ndivisors(M);           // At the same time, we remember
                                                    // how many divisors every divisor of N has,
                                                    // since we'll want this later.
    }                                               

    for(unsigned int k = 0; k < sublevels.size(); k++) {
        int M = sublevels[k];
        if(M == 1) continue;
        for(unsigned int j = k; j < sublevels.size(); j++) {
            int N = sublevels[j];
            vector<int> Ndivisors = divisors(N);
            if(N % M != 0) continue;
            // At this point, we are assuming that traces[M] has already been sieved,
            // and we are removing the contribution of traces[M] from from traces[N].

            for(int n = max(1, start); n < end; n++) {
                if(chi_values[N][n % N] != 0) {     // This is just a check on the GCD.
                    if(N == M) continue;
                    long z = nmod_mul(divisor_counts[N/M], traces[M][n], modp);
                    traces[N][n] = nmod_sub(traces[N][n], z, modp);
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
                            if(mu == 1) {
                                long z = nmod_mul(divisor_counts[x], chi_values[q][b % q], modp);
                                z = nmod_mul(z, b, modp);
                                z = nmod_mul(z, traces[M][n/(b*b)], modp);
                                traces[N][n] = nmod_sub(traces[N][n], z, modp);
                            }
                            else if(mu == -1) {
                                long z = nmod_mul(divisor_counts[x], chi_values[q][b % q], modp);
                                z = nmod_mul(z, b, modp);
                                z = nmod_mul(z, traces[M][n/(b*b)], modp);
                                traces[N][n] = nmod_add(traces[N][n], z, modp);
                            }
                        }
                    }
                }
            }
        }
    }

    delete [] divisor_counts;
}

int newspace_bases_weight2_modp(nmod_mat_t * bases, int& ncoeffs, int level, long& p, DirichletCharacter& chi, int extra_rows, int verbose) {
    //
    // Compute ncoeffs coefficients of bases for S^new_2(M, chi) for all M with q | M | level,
    // and put these in bases[M].
    //
    // Through a quirk of convenience, we also return the dimension of
    // S_2(level, chi).
    //

    if(!chi.is_even()) {                                            // 
        int q = chi.conductor();                                    // 
        for(int M : divisors(level)) {                              // We do something somewhat sensible if
            if(M % q == 0 && M > 1) nmod_mat_init(bases[M], 0, 0, 2);// we know that the dimension of the space
        }                                                           // is zero.
        return 0;                                                   //
    }                                                               //

    long primitive_index;
    int q = chi.conductor(&primitive_index);
    long ** chi_values = new long*[level + 1];

    if(verbose > 0) {
        cerr << "Computing bases of" << endl;
    }
    vector<int> sublevels;                     // 
    for(int M : divisors(level)) {             // We fill sublevels with all of the M such that 
        if(M % q == 0) {                       // q | M | level.
            if(M > 1) sublevels.push_back(M);  //
            DirichletGroup G(M);
            DirichletCharacter psi = G.character(G.index_from_primitive_character(q, primitive_index));
            chi_values[M] = new long[M];
            psi.values_mod_p(p, chi_values[M]);
            if(verbose > 0) {
                cerr << "    " << "S_2^new(" << M << ", " << psi.m << ") mod " << p << endl;
            }
        }
    }

    long cuspform_dimension;
    trace_Tn_modp_unsieved_weight2(&cuspform_dimension, 1, 2, level, p, chi_values[level], chi);  // XXX
                                                                            // I'm assuming that the chosen prime
                                                                            // is large enough that the dimension
                                                                            // is actually the same as the dimension
                                                                            // mod p.

    if(ncoeffs == 0) ncoeffs = cuspform_dimension + extra_rows + 20;        // It is unlikely that we actually need this many
    int max_trace = (cuspform_dimension + extra_rows + 20)* ncoeffs;        // traces, but it doesn't hurt much, and maybe it
                                                                            // better to not have to go back for more.
    
    vector<long> * traces = new vector<long>[level + 1];

    for(int M : sublevels) {
        if(verbose > 0)
            cerr << "Computing " << max_trace << " traces of Tn on " << " S_2(" << M << ", chi)." << endl;
        traces[M] = vector<long>(max_trace + 1);                                                         // We compute the (unsieved) traces
        trace_Tn_modp_unsieved_weight2(traces[M].data(), 0, max_trace + 1, M, p, chi_values[M], chi, verbose);   // for each sublevel.
    }                                                                                                   //

    if(verbose > 0) cerr << "Sieving..." << endl;
    sieve_trace_Tn_modp_on_weight2_for_newspaces(traces, 0, max_trace + 1, level, p, chi_values, chi);  // Then we sieve.

    if(verbose > 0) cerr << "And computing bases..." << endl;
    for(int M : sublevels) {                                                // And now that we have the traces on the
        int new_dimension = traces[M][1];                                   // new subspace we can actually start writing
        int nrows = new_dimension + extra_rows;                             // down bases.
        nmod_mat_init(bases[M], nrows, ncoeffs, p);                         //
                                                                            // The inside of this loop is something that
                                                                            // could naturally be a separate function, but
                                                                            // if we try to do that it will get messy if
                                                                            // we need to backup and compute more traces.
        for(int m = 1; m <= nrows; m++) {
            for(int n = 1; n < ncoeffs; n++) {
                long z = trace_TmTn_mod_p(traces[M].data(), chi_values[M], 2, M, m, n, p);
                nmod_mat_entry(bases[M], m-1, n) = z;
            }
        }

        if(nmod_mat_rank(bases[M]) == new_dimension) {
            //cerr << "At level " << M << " character (" << level << ", " << chi.m << ") using p = " << p;
            //cerr << ". Used first " << new_dimension << " rows." << endl;
        }

        int compute_more_count = 0;
        if(nmod_mat_rank(bases[M]) < new_dimension) {
            if(verbose > 0) {
                cerr << "For level " << M << " the first rows were not enough." << endl
                     << "Now making sure that the rank increases with each row that we add." << endl;
                cerr.flush();
            }
            nmod_mat_zero(bases[M]);
            int rank = 0;
            // It is possible that we don't actualy have a basis yet. I think we
            // usually have a basis, but sometimes we'll need to build up the basis
            // more conservatively, making sure that the rank increases with each
            // new row, which we do below.
            //
            // XXX: The following code has only been lightly tested.
            //
            int j = 0; // We'll be putting the trace of TmTn into the nth column
            int m = 1; // in the jth row.
            while(j < nrows) {
                if(verbose > 1) {
                    cerr << "Trying to put " << m << "th possible basis element into the " << j << "th row." << endl;
                    cerr << "Looking for rank " << new_dimension << endl;
                    cerr.flush();
                }
                // We need to compute more traces now, since
                // need more than we originally expected.
                int max_trace_needed_next = m * ncoeffs;
                if(max_trace_needed_next > max_trace) {
                    cerr << "At level " << M << " character (" << level << ", " << chi.m << ") using p = " << p;
                    cerr << ", computing more traces: looking for rank " << new_dimension;
                    cerr << ", currently at rank " << rank << endl;
                    cerr << " already have " << max_trace << " traces. ncoeffs = " << ncoeffs << endl;
                    for(int M : sublevels) {
                        traces[M].resize(max_trace_needed_next + 1);
                        trace_Tn_modp_unsieved_weight2(&traces[M][max_trace + 1], max_trace + 1, max_trace_needed_next + 1, M, p, chi_values[M], chi);
                    }
                    sieve_trace_Tn_modp_on_weight2_for_newspaces(traces, max_trace + 1, max_trace_needed_next + 1, level, p, chi_values, chi);
                    max_trace = max_trace_needed_next;
                    compute_more_count++;
                }

                for(int n = 1; n < ncoeffs; n++) {
                    long z = trace_TmTn_mod_p(traces[M].data(), chi_values[M], 2, M, m, n, p);
                    nmod_mat_entry(bases[M], j, n) = z;
                }
                rank = nmod_mat_rank(bases[M]);
                if(rank == j+1 || rank == new_dimension) {
                    //cerr << "At level " << M << " character (" << level << ", " << chi.m << ") using p = " << p;
                    //cerr << ". Using row " << m << endl;
                    //cerr << "rank is " << rank << endl;

                    j++;
                }
                m++;
                if(compute_more_count > 100) {
                    cerr << "Giving up on finding full rank." << endl;
                    j = nrows;
                }
            }
        }
    }
 
    for(int M : sublevels) {
        delete [] chi_values[M];
    }

    delete [] chi_values;
    delete [] traces;

    return cuspform_dimension;
}

void cuspform_basis_weight2_modp(nmod_mat_t basis, int ncoeffs, int level, long& p, DirichletCharacter& chi, int verbose) {
    nmod_mat_t * bases_for_new_spaces = new nmod_mat_t[level + 1];
    int cuspform_dimension = newspace_bases_weight2_modp(bases_for_new_spaces, ncoeffs, level, p, chi, 0, verbose);
    nmod_mat_init(basis, cuspform_dimension, ncoeffs, p);

    int q = chi.conductor();
    vector<int> sublevels;                     // 
    for(int M : divisors(level)) {             // We fill sublevels with all of the M such that 
        if(M % q == 0) {                       // q | M | level.
            if(M > 1) sublevels.push_back(M);  //
        }
    }

    int r = 0;
    for(int M1 : sublevels) {
        for(int M2 : divisors(level/M1)) {
            for(int l = 0; l < nmod_mat_nrows(bases_for_new_spaces[M1]); l++) {
                for(int j = 1; j*M2 < ncoeffs; j++) {
                    nmod_mat_entry(basis, r, j * M2) = nmod_mat_entry(bases_for_new_spaces[M1], l, j);
                }
                r++;
            }
        }
    }
    for(int M : sublevels) {
        nmod_mat_clear(bases_for_new_spaces[M]);
    }
    delete [] bases_for_new_spaces;
}
