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


//#undef nmod_mat_entry
//static mp_limb_t& nmod_mat_entry(nmod_mat_t mat, int i, int j) {
//    if(i >= mat->r) { cout << "nmod_mat_entry() called with i out of range."; exit(0); }
//    if(j >= mat->c) { cout << "nmod_mat_entry() called with j out of range."; exit(0); }
//    return mat->rows[i][j];
//}

int trace_TmTn_mod_p(int * traces,
                     int * chi_values,
                     int k,
                     int level,
                     int m,
                     int n,
                     int p) {
    //
    // Return the trace mod p of TmTn acting on S_k(chi).
    //

    int g = GCD(m,n);
    int TnTm = 0;

    for(int d = 1; d <= g; d++) {
        if(g % d != 0) continue; // this is silly, but this will usually
                                 // be a short sum
        if(GCD(d,level) != 1) continue;
        
        long dk = PowerMod(d, k-1, p);
        TnTm += ((traces[m/d * n/d] * dk) % p) * (long)chi_values[d % level] % p;
        TnTm %= p;
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

void trace_Tn_unsieved_weight2(complex<double> * traces, int start, int end, int level, DirichletCharacter& chi) {
    // We expect that p is not 2 or 3.

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
                        if(GCD(c/g, level/c) == 1) {            // We are doing a CRT lift
                            y = CRT(d, n/d, c/g, level/c);      // here, but the moduli might
                        }                                       // not be coprime. Hence the
                        else { // GCD(c, (N/c)/g) == 1          // complication.
                            y = CRT(d, n/d, c/g, (level/c)/g);
                        }
                        a = a + (complex<double>)euler_phi(g) * chi.value(y);
                    }
                }
            }
            if(d*d == n) a /= 2.0;
            traces[n] -=  (complex<double>)d * a;
        }
    }
    
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

    delete [] psi_table;
    for(int k = 1; k < level + 1; k++) {
        if(level % k == 0) {
            delete [] gcd_tables[k];
        }
    }
    delete [] gcd_tables;
}

void trace_Tn_modp_unsieved_weight2(int * traces, int start, int end, int level, int p, int * chi_values, DirichletCharacter& chi) {
    // We expect that p is not 2 or 3.

    traces = traces - start; // We adjust the traces pointer so that we
                             // can work with traces[start] through
                             // trace[end-1] instead of subtracting start all
                             // the time.

    int one_over_twelve = InvMod(12, p);
    int one_over_three = InvMod(3, p);
    int one_over_two = InvMod(2, p);

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
        traces[n] = 0;
    }

    int A1 = ((psi_table[level] % p) * one_over_twelve) % p;
    int z = (int)sqrt(start);                               //
    while(z*z < start) z++;                                 // A1 is only nonzero when n is a square,
    for( ; z*z < end; z++) {                                // and in the weight 2 case it hardly
        traces[z*z] += A1 * chi_values[z % level] % p;    // depends on n at all.
    }


    // computation of A2
    //cout << "computing A2 for level " << level << endl;
    
    int t = sqrt(4*end);
    if(t*t == 4*end) t--;
    for(t = -t; (long)t*t < 4l*end; t++) {
        for(int x = 0; x < level; x++) {
            int chi = chi_values[x];
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
                            int z = (long)psi_table[level]/psi_table[level/g] * (long)chi * (long)classnumbers[D/(f*f)] % p;
                            if(D/(f*f) == 3) z = (z * one_over_three) % p;
                            if(D/(f*f) == 4) z = (z * one_over_two) % p;
                            traces[n] -= (z * one_over_two) % p;
                            traces[n] %= p;
                            if(traces[n] < 0) traces[n] += p;
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

    //cout << "computing A3 for level " << level << endl;
    // computation of A3
    for(int d = 1; d*d < end; d++) {            // In A3, we'll get a contribution
        int starting_n = start;
        if(d*d > starting_n) {
            starting_n = d*d;
        }
        else if(starting_n % d != 0) {
            starting_n += (d - starting_n % d);
        }
        for(int n = starting_n; n < end; n += d) {     // to traces[n] for each d that divides n
            int a = 0;                          // as long as d^2 <= n
            for(int c : divisors_of_level) {
                int z = (n/d - d) % (level/conductor);
                if(z < 0) z += (level/conductor);
                int g = gcd_tables[c][level/c % c];
                if(gcd_tables[level/conductor][z] % g == 0) {
                    int y;
                    if(g == 1) {
                        y = CRT(d, n/d, c, level/c);
                        a = a + chi_values[y];
                        a %= p;
                    }
                    else {
                        if(GCD(c/g, level/c) == 1) {            // We are doing a CRT lift
                            y = CRT(d, n/d, c/g, level/c);      // here, but the moduli might
                        }                                       // not be coprime. Hence the
                        else { // GCD(c, (N/c)/g) == 1          // complication.
                            y = CRT(d, n/d, c, (level/c)/g);    //
                        }                                       // XXX: I've got this computation wrong
                        a = a + euler_phi(g) * chi_values[y];   // a few times. I hope it is right now.
                        a %= p;
                    }
                }
            }
            if(d*d == n) a = (a * one_over_two) % p;
            traces[n] -=  d * a;
            traces[n] %= p;
            if(traces[n] < 0) traces[n] += p;
        }
    }
    
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
                    traces[n] += t;
                    traces[n] %= p;
                }
            }
        }
    }

    //cout << "cleaning up for level " << level << endl;
    delete [] psi_table;
    for(int k = 1; k < level + 1; k++) {
        if(level % k == 0) {
            delete [] gcd_tables[k];
        }
    }
    delete [] gcd_tables;
}

void sieve_trace_Tn_modp_on_weight2_for_newspaces(vector<int> * traces, int start, int end, int level, int p, int ** chi_values, DirichletCharacter& chi) {
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

    for(int k = 0; k < sublevels.size(); k++) {
        int M = sublevels[k];
        if(M == 1) continue;
        for(int j = k; j < sublevels.size(); j++) {
            int N = sublevels[j];
            vector<int> Ndivisors = divisors(N);
            if(N % M != 0) continue;
            // At this point, we are assuming that traces[M] has already been sieved,
            // and we are removing the contribution of traces[M] from from traces[N].

            for(int n = max(1, start); n < end; n++) {
                if(chi_values[N][n % N] != 0) {     // This is just a check on the GCD.
                    if(N == M) continue;
                    traces[N][n] -= divisor_counts[N/M] * traces[M][n];
                    traces[N][n] %= p;
                    if(traces[N][n] < 0) traces[N][n] += p;
                }
                else {
                        // I'm sure that we can this more efficiently, but let's
                            // try the slow way first.
                    int x = N/M; while(GCD(x, n) > 1) {x = x/GCD(x,n);}  // now x = (N/M) / GCD(N/M, n^oo)
                    //int y = n;   while(GCD(y, M) > 1) {y = y/GCD(y,M);}  // and y = n / GCD(n, M^oo)
                    int y = N; while(GCD(y, M) > 1) {y = y/GCD(y,M);}  // now y = N / GCD(N, M^oo)
                    //for(int b : Ndivisors) {
                    for(int b : divisors(y)) {
                        if(b == 1 && N == M) continue;
                        if(n % (b*b) == 0) {
                            traces[N][n] -= divisor_counts[x] * chi_values[q][b % q] * mobius(b) * b * traces[M][n/(b*b)];
                            traces[N][n] %= p;
                            if(traces[N][n] < 0) traces[N][n] += p;
                            /*
                            int_factorization_t bfac;
                            factor(b, bfac);
                            int P = 1;
                            for(int k = 0; k < bfac.nfactors; k++) {
                                int p = bfac.factors[k].p;
                                int e = bfac.factors[k].e;
                                int a = 1;
                                if(e != 1) {
                                    a = 0; int x = n; while(x % p == 0) {x /= p; a++;}
                                }
                                if(!(e == 1 || e*2 == a)) {
                                    P = 0;
                                    break;
                                }
                                if(e == 1 && a == 2) {
                                    P *= 2;
                                }
                                P *= -1;
                            }
                            traces[N][n] -= chi_values[q][b % q] * P * traces[M][n/(b*b)];
                            traces[N][n] %= p;
                            if(traces[N][n] < 0) traces[N][n] += p;
                            */
                        }
                    }
                }
            }
            //for(long M2 : divisors(N/M)) {  // Might need to do something like this.
            //}
        }
    }

    delete [] divisor_counts;
}

int newspace_bases_weight2_modp(nmod_mat_t * bases, int& ncoeffs, int level, int& p, DirichletCharacter& chi, int extra_rows) {
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
            if(M % q == 0) nmod_mat_init(bases[level], 0, 0, 2);    // we know that the dimension of the space
        }                                                           // is zero.
        return 0;                                                   //
    }                                                               //

    long primitive_index;
    int q = chi.conductor(&primitive_index);
    int ** chi_values = new int*[level + 1];

    vector<int> sublevels;                     // 
    for(int M : divisors(level)) {             // We fill sublevels with all of the M such that 
        if(M % q == 0) {                       // q | M | level.
            if(M > 1) sublevels.push_back(M);  //
            DirichletGroup G(M);
            DirichletCharacter psi = G.character(G.index_from_primitive_character(q, primitive_index));
            chi_values[M] = new int[M];
            psi.values_mod_p(p, chi_values[M]);
        }                                           
    } 
    
    int cuspform_dimension;
    trace_Tn_modp_unsieved_weight2(&cuspform_dimension, 1, 2, level, p, chi_values[level], chi);  // XXX
                                                                            // I'm assuming that the chosen prime
                                                                            // is large enough that the dimension
                                                                            // is actually the same as the dimension
                                                                            // mod p.

    if(ncoeffs == 0) ncoeffs = cuspform_dimension + extra_rows + 20;        // It is unlikely that we actually need this many
    int max_trace = (cuspform_dimension + extra_rows + 20)* ncoeffs;        // traces, but it doesn't hurt much, and maybe it
                                                                            // better to not have to go back for more.
    
    vector<int> * traces = new vector<int>[level + 1];

    for(int M : sublevels) {                                                                            //
        traces[M] = vector<int>(max_trace + 1);                                                         // We compute the (unsieved) traces
        trace_Tn_modp_unsieved_weight2(traces[M].data(), 0, max_trace + 1, M, p, chi_values[M], chi);   // for each sublevel.
    }                                                                                                   //

    sieve_trace_Tn_modp_on_weight2_for_newspaces(traces, 0, max_trace + 1, level, p, chi_values, chi);  // Then we sieve.

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
                int z = trace_TmTn_mod_p(traces[M].data(), chi_values[M], 2, M, m, n, p);
                nmod_mat_entry(bases[M], m-1, n) = z;
            }
        }

        if(nmod_mat_rank(bases[M]) == new_dimension) {
            //cerr << "At level " << M << " character (" << level << ", " << chi.m << ") using p = " << p;
            //cerr << ". Used first " << new_dimension << " rows." << endl;
        }

        if(nmod_mat_rank(bases[M]) < new_dimension) {
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
                        trace_Tn_modp_unsieved_weight2(&traces[M][max_trace + 1], max_trace + 1, max_trace_needed_next + 1, level, p, chi_values[M], chi);
                    }
                    sieve_trace_Tn_modp_on_weight2_for_newspaces(traces, max_trace + 1, max_trace_needed_next + 1, level, p, chi_values, chi);
                    max_trace = max_trace_needed_next;
                }

                for(int n = 1; n < ncoeffs; n++) {
                    int z = trace_TmTn_mod_p(traces[M].data(), chi_values[M], 2, M, m, n, p);
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

void cuspform_basis_weight2_modp(nmod_mat_t basis, int ncoeffs, int level, int& p, DirichletCharacter& chi) {
    nmod_mat_t * bases_for_new_spaces = new nmod_mat_t[level + 1];
    int cuspform_dimension = newspace_bases_weight2_modp(bases_for_new_spaces, ncoeffs, level, p, chi);
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
                for(int j = 1; j < ncoeffs/M2; j++) {
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

/*
int main(int argc, char ** argv) {
    int level = atoi(argv[1]);
    int character = atoi(argv[2]);
    //int start = atoi(argv[3]);
    //int end = atoi(argv[4]);

    //build_factor_table(4*max(level, end) + 1);
    build_factor_table(2000000);

    int p = 0;
    int * chi_values = new int[level];
    DirichletGroup G(level);
    DirichletCharacter chi = G.character(character);

    nmod_mat_t basis;
    cuspform_basis_weight2_modp(basis, 0, level, p, chi);
    
    //nmod_mat_t * bases = new nmod_mat_t[level + 1];
    //newspace_bases_weight2_modp(bases, 0, level, p, chi);
    //int * traces = new int[11];
    //chi.values_mod_p(p, chi_values);
    //trace_Tn_modp_unsieved_weight2(traces, 0, 10, level, p, chi_values, chi);
    //for(int k = 0; k < 10; k++) cout << traces[k] << " ";
    //cout << endl;
    return 0;
}
*/
