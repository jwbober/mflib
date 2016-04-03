#include <iostream>
#include <map>
#include <tuple>
#include <vector>
#include <complex>

#include "modform_cc.h"
#include "classnumbers.h"
using std::map;
using std::pair;
using std::tuple;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::real;
using std::imag;
using std::complex;

using Eigen::ComplexEigenSolver;

typedef tuple<int,int,int> space_desc_t;
static map<space_desc_t, cuspforms_cc*> cache;

cuspforms_cc * get_cuspforms_cc(DirichletCharacter &chi, int weight, int verbose) {
    if(verbose > 2) cout << "getting space cuspforms_cc(" << chi.parent->q << ", " << weight << ", " << chi.m << ")" << endl;
    auto result = cache.find( space_desc_t(chi.parent->q, chi.m, weight) );
    if(result != cache.end()) return result->second;

    cuspforms_cc * S = new cuspforms_cc(chi, weight, verbose);
    cache[space_desc_t(chi.parent->q, chi.m, weight)] = S;
    return S;
}

bool has_unique_entries(const complex<double> * zz, int length, double eps = 1e-7) {
    complex<double> * z = new complex<double>[length];
    std::copy(zz, zz + length, z);
    std::sort(z, z + length,
        [](complex<double> x, complex<double> y){
            return real(x) < real(y);
        }
    );
    complex<double> prev = z[0];
    for(int k = 1; k < length; k++) {
        int j = k;
        while(j < length && abs(real(z[j]) - real(prev)) < eps) {
            if(abs(z[j] - prev) < eps) return false;
            j++;
        }
        prev = z[k];
    }
    delete [] z;
    return true;
}

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

double cuspforms_cc::evalpoly(long tt, long nn) {
    // translated from Ralph's code
    double t = tt;
    double n = nn;
    int k = weight - 2;
    switch(k) {
        case 0: return 1;
        case 1: return t;
        case 2: { // t^2 - n
            return t*t - n;
            }
        case 3: { // t^3 - 2tn
            return t*t*t - 2*t*n;
            }
        case 4: { // t^4 - 3nt^2 + n^2
            return t*t*(t*t - 3*n) + n*n;
        }
    case 5: {
            // =t^5 - 4*t^3*n + 3*t*n^2;
            return t*t*t*t*t - 4*t*t*t*n + 3*t*n*n;
        }
    //case 10:
    //    i64 tt = t*t;
    //    return - n*n*n*n*n + tt * (15*n*n*n*n + tt * (-35*n*n*n + tt * (28*n*n + tt * (-9*n + tt))));
    }

    double val[2];
    val[0] = t*t*t*t - 3*t*t*n + n*n;
    val[1] = t*t*t*t*t - 4*t*t*t*n + 3*t*n*n;

    for(int i = 6; i<= k; i++) {
        val[i % 2] = t*val[(i + 1) % 2] - n *val[i % 2];
    }
    return val[k % 2];
}

complex<double> cuspforms_cc::trace(int n) {
    if(n > 1 && new_dimension() == 0) return 0;
    if(n < traces.size()) return traces[n];

    int new_end = std::max( (int)(n * 1.1), n + 10 );
    compute_traces(new_end);
    return traces[n];
}

complex<double> cuspforms_cc::trace_TnTm(int n, int m) {
    //
    // Return the trace mod p of TmTn acting on S_k(chi).
    //

    if(n*m == 0) return 0.0;
    int g = GCD(m,n);
    complex<double> TnTm = 0;

    for(int d = 1; d <= g; d++) {
        if(g % d != 0) continue; // this is silly, but this will usually
                                 // be a short sum
        if(GCD(d,level) != 1) continue;

        //long dk = PowerMod(d, k-1, p);
        double dk = pow(d, weight - 1);
        complex<double> z = trace(m/d * n/d) * dk * chi_values[d % level];
        TnTm += z;
    }
    return TnTm;
}

int cuspforms_cc::new_dimension() {
    return (int)(round)(real(trace(1)));
}

cmatrix_t cuspforms_cc::newspace_basis(int ncoeffs) {
    newspace_basis_data();
    int d = new_dimension();

    compute_traces(ncoeffs * basis_rows[d-1]);
    cmatrix_t basis(d, ncoeffs);
    for(int k = 0; k < d; k++) {
        int n = basis_rows[k];
        for(int m = 0; m < ncoeffs; m++) {
            basis(k, m) = trace_TnTm(n, m);
        }
    }
    return basis;
}

cmatrix_t cuspforms_cc::newforms(int ncoeffs) {
    int dim = new_dimension();
    if(dim == 0) {return cmatrix_t(0,0);}
    int coefficients_needed_for_full_rank = newspace_basis_data()[dim - 1] + 1;
    ncoeffs = std::max(ncoeffs, coefficients_needed_for_full_rank);
    cmatrix_t basis = newspace_basis(ncoeffs);

    cmatrix_t transformed_basis = cmatrix_t::Zero(dim, coefficients_needed_for_full_rank);

    auto LU = basis.leftCols(coefficients_needed_for_full_rank).transpose().fullPivLu();
    bool found_unique_eigenvalues = false;
    int n = 2;
    cmatrix_t basis_transformation;

    while(!found_unique_eigenvalues) {
        int max_coeff_needed = coefficients_needed_for_full_rank * n;
        if(max_coeff_needed > ncoeffs) {
            cerr << "need to compute more coefficients." << endl;
            exit(0);
        }
        for(int j = 0; j < dim; j++) {
            for(int m = 1; m < coefficients_needed_for_full_rank; m++) {
                for(auto d : divisors(GCD(m, n))) {
                    transformed_basis(j, m) += chi_values[d % level] * pow((double)d, weight - 1) * basis(j, m*n/(d*d));
                }
            }
        }
        //cout << "transformed basis:" << endl;
        //cout << transformed_basis << endl;
        cmatrix_t Tn_sum = LU.solve(transformed_basis.transpose());
        //cout << "Tn_sum:" << endl;
        //cout << Tn_sum << endl;
        ComplexEigenSolver<cmatrix_t> ces;
        ces.compute(Tn_sum);
        //cout << "eigenvalues:" << endl;
        //cout << ces.eigenvalues() << endl;
        if(has_unique_entries(ces.eigenvalues().data(), ces.eigenvalues().size())) {
            found_unique_eigenvalues = true;
        }
        else {
           found_unique_eigenvalues = false;
        }
        if(found_unique_eigenvalues) {
           cmatrix_t z = ces.eigenvectors();
           //cout << "eigenvectors:" << endl;
           //cout << z << endl;
           z.transposeInPlace();
           cmatrix_t eigenforms = z * basis.leftCols(dim + 1);
           for(int k = 0; k < z.rows(); k++) {
               z.row(k) = z.row(k)/eigenforms(k,1);
           }
           basis_transformation = z;
        }
        n++;
    }

    return basis_transformation * basis;
}

const vector<int>& cuspforms_cc::newspace_basis_data() {
    if(basis_rows.size() > 0) return basis_rows;
    int d = new_dimension();
    if(d == 0) return basis_rows;
    basis_rows = modp_space->newspace_basis_data();
    return basis_rows;
}

void cuspforms_cc::compute_traces(int end) {
    if(end < traces.size()) return;
    for(cuspforms_cc *subspace : subspaces) {
        subspace->compute_traces(end);
    }
    int start = traces.size();
    traces.resize(end, 0.0);

    complex<double> A1 = psi_table[level]/12.0 * (weight - 1.0);
    int z = (int)sqrt(start);                               //
    while(z*z < start) z++;                                 // A1 is only nonzero when n is a square,
    for( ; z*z < end; z++) {                                // and in the weight 2 case it hardly
                                                            // depends on n at all.
        double npow = pow(z, weight - 2);
        traces[z*z] += npow * A1 * chi_values[z % level];
    }

    long print_interval = 0;
    if(verbose > 0) { cerr << "A2:"; print_interval = round(2 * sqrt(4*end)/70); print_interval = std::max(print_interval, 1l); }

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
        for(long x = 0; x < level; x++) {
            complex<double> chi_value = chi_values[x];
            if(chi_value == 0.0) continue;
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
                double polyterm = evalpoly(t, n);
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
                        complex<double> z = (complex<double>)(psi_table[level]/psi_table[level/g]);
                        z *= chi_value;
                        z *= classnumbers[D/(f*f)];
                        if(D/(f*f) == 3) z /= 3.0; // z = nmod_mul(z, one_over_three, modp); // (z * one_over_three) % p;
                        if(D/(f*f) == 4) z /= 2.0; // z = nmod_mul(z, one_over_two, modp);   //(z * one_over_two) % p;
                        z /= 2.0;
                        z *= polyterm;
                        traces[n] -= z;
                        //z = nmod_mul(z, one_over_two, modp);
                        //traces[n] = nmod_sub(traces[n], z, modp);
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

    if(verbose > 0) { cerr << "A3:"; print_interval = (long)std::max(round(sqrt(end)/70), 1.0); }
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
            complex<double> a = 0;                     // as long as d^2 <= n
            for(int c : divisors_of_level) {
                int z = (n/d - d) % (level/conductor);
                if(z < 0) z += (level/conductor);
                int g = gcd_tables[c][level/c % c];
                if(gcd_tables[level/conductor][z] % g == 0) {
                    long y;
                    if(g == 1) {
                        y = CRT(d, n/d, c, level/c);
                        a += chi_values[y];
                    }
                    else {
                        unsigned long g1, u, v;                     // We are doing a CRT lift
                        g1 = n_xgcd(&u, &v, c, level/c);            // here, but the moduli might
                        y = (d + c * u * (n/d - d)/g1) % (level/g); // not be coprime. Hence the
                        if(y < 0) y += (level/g);                   // complication.

                                                                    // XXX: I've got this computation wrong
                        //a = a + phi_table[g] * chi_values[y];       // a few times. I hope it is right now.

                        a += (double)phi_table[g] * chi_values[y];
                        //NMOD_ADDMUL(a, phi_table[g], chi_values[y], modp);
                        //a %= p;

                        //cout << "y = " << y << " " << d << " " << n/d << " " << c << " " << level/c << endl;
                    }
                }
            }
            if(d*d == n) a /= 2.0; //a = nmod_mul(a, one_over_two, modp);
            //a *= (double)d; //a = nmod_mul(a, d, modp);
            a *= pow(d, weight - 1);
            traces[n] -= a;
            //traces[n] = nmod_sub(traces[n], a, modp);       // traces -= d*a mod p
        }
    }
    if(verbose > 0) { cerr << endl; cerr.flush(); }
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
                    traces[n] += (double)t;
                    //traces[n] = nmod_add(traces[n], t, modp);
                    //traces[n] += t;
                    //traces[n] %= p;
                }
            }
        }
    }

    // now we sieve.

    int q = conductor;
    for(int k = 0; k < sublevels.size(); k++) {
        cuspforms_cc * subspace = subspaces[k];
        int N = level;
        int M = sublevels[k];
        for(int n = start; n < end; n++) {
            if(chi_values[n % level] != 0.0) {     // This is just a check on the GCD.
                traces[n] -= (complex<double>)divisor_counts[N/M] * subspace->traces[n];
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
                        traces[n] -= mu * (double)divisor_counts[x] * chip_values[b % q] * pow((double)b, weight - 1) * subspace->traces[n/(b*b)];
                    }
                }
            }
        }
    }

}
