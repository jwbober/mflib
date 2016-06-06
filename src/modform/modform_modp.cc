#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include "modform_modp.h"
#include "classnumbers.h"
using std::map;
using std::pair;
using std::tuple;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;

typedef tuple<int,int,int,long> space_desc_t;
static map<space_desc_t, cuspforms_modp*> cache;
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

int cuspforms_modp::new_dimension() {
    return trace(1);
}

long cuspforms_modp::evalpoly(long t, long n) {
    // translated from Ralph's code
    int k = weight - 2;
    t = t % p; if(t < 0) t += p;
    n = n % p; if(n < 0) n += p;
    switch(k) {
        case 0: return 1;
        case 1: return t;
        case 2: { // t^2 - n
            long z = nmod_mul(t, t, modp);
            return nmod_sub(z, n, modp);
            }
        case 3: { // t^3 - 2tn
            long t2 = nmod_mul(t, t, modp);
            long two_n = nmod_mul(n, 2, modp);
            long z = nmod_sub(t2, two_n, modp);
            return nmod_mul(t, z, modp);
            }
        case 4: { // t^4 - 3nt^2 + n^2
            long t2 = nmod_mul(t, t, modp);
            long t4 = nmod_mul(t2, t2, modp);
            long n2 = nmod_mul(n, n, modp);
            long z = nmod_mul(t2, 3, modp);
            z = nmod_mul(z, n, modp);
            z = nmod_sub(t4, z, modp);
            z = nmod_add(z, n2, modp);
            return z;
        }
    case 5: {
            long t2 = nmod_mul(t, t, modp);
            long t3 = nmod_mul(t2, t, modp);
            long t5 = nmod_mul(t2, t3, modp);
            long n2 = nmod_mul(n, n, modp);
            long z1 = nmod_mul(3, t, modp);
            z1 = nmod_mul(z1, n2, modp);
            long z2 = nmod_mul(4, t3, modp);
            z2 = nmod_mul(n, z2, modp);
            long z = nmod_sub(t5, z2, modp);
            return nmod_add(z, z1, modp);
        }
    //case 10:
    //    i64 tt = t*t;
    //    return - n*n*n*n*n + tt * (15*n*n*n*n + tt * (-35*n*n*n + tt * (28*n*n + tt * (-9*n + tt))));
    }

    long val[2];
    long t2 = nmod_mul(t, t, modp);
    long t4 = nmod_mul(t2, t2, modp);
    long n2 = nmod_mul(n, n, modp);
    long z = nmod_mul(t2, 3, modp);
    z = nmod_mul(z, n, modp);
    z = nmod_sub(t4, z, modp);
    z = nmod_add(z, n2, modp);
    //val[0]=t*t*t*t - 3*t*t*n + n*n;
    val[0] = z;

    long t3 = nmod_mul(t2, t, modp);
    long t5 = nmod_mul(t2, t3, modp);
    long z1 = nmod_mul(3, t, modp);
    z1 = nmod_mul(z1, n2, modp);
    long z2 = nmod_mul(4, t3, modp);
    z2 = nmod_mul(n, z2, modp);
    z = nmod_sub(t5, z2, modp);
    z = nmod_add(z, z1, modp);
    // val[1]=t*t*t*t*t - 4*t*t*t*n + 3*t*n*n;
    val[1] = z;
    for(int i = 6; i<= k; i++) {
        long z1 = nmod_mul(t, val[(i + 1) % 2], modp);
        long z2 = nmod_mul(n, val[i % 2], modp);
        val[i % 2] = nmod_sub(z1, z2, modp);
        //val[i % 2] = t*val[(i + 1) % 2] - n *val[i %2];
    }
    return val[k % 2];
}

cuspforms_modp * get_cuspforms_modp(DirichletCharacter &chi, int weight, long p, int verbose) {
    long order = order_mod(chi.m, chi.parent->q);
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

    if(verbose > 2) cout << "getting space cuspforms_modp(" << chi.parent->q << ", " << weight << ", " << chi.m << ")" << " with p == " << p << endl;
    auto result = cache.find( space_desc_t(chi.parent->q, chi.m, weight, p) );
    if(result != cache.end()) return result->second;

    cuspforms_modp * S = new cuspforms_modp(chi, weight, p, verbose);
    cache[space_desc_t(chi.parent->q, chi.m, weight, p)] = S;
    return S;
}

void clear_cuspforms_modp() {
    for(auto x : cache) {
        delete x.second;
    }
    cache.clear();
}

void cuspforms_modp::newforms(nmod_mat_t forms, int ncoeffs) {
    int dim = new_dimension();
    vector<int> basis_data = newspace_basis_data();
    int coefficients_needed_for_full_rank = basis_data[dim - 1] + 1;
    ncoeffs = std::max(ncoeffs, coefficients_needed_for_full_rank);
    nmod_mat_t basis;
    newspace_basis(basis, ncoeffs);
    nmod_mat_t smallbasis; // flint wants all the matrices to be square
                           // for linear equation solving, so we're going
                           // to pick out a piece of the basis that we know
                           // gives full rank. At the same time, we'll transpose
                           // it.
    nmod_mat_t Tn;
    nmod_mat_t transformed_basis;
    nmod_mat_init(smallbasis, dim, dim, p);
    nmod_mat_init(transformed_basis, dim, dim, p);
    nmod_mat_init(Tn, dim, dim, p);


    for(int k = 0; k < dim; k++) {
        for(int j = 0; j < dim; j++) {
            nmod_mat_entry(smallbasis, j, k) = nmod_mat_entry(basis, k, basis_data[j]);
        }
    }

    int n = 2;
    nmod_mat_t X, Y;
    nmod_mat_init(X, dim, dim, p);
    nmod_mat_init(Y, dim, dim, p);
    while(n < 20) {
        nmod_mat_zero(transformed_basis);
        int max_coeff_needed = coefficients_needed_for_full_rank * n;
        if(max_coeff_needed > ncoeffs) {
            cerr << "need to compute more coefficients." << endl;
            exit(0);
        }
        for(int j = 0; j < dim; j++) {
            for(int k = 0; k < dim; k++) {
                int m = basis_data[k];
                for(auto d : divisors(GCD(m, n))) {
                    d = d % p;
                    long z = chi_values[d % level];
                    z = nmod_mul(z, nmod_pow_ui(d, weight - 1, modp), modp);
                    z = nmod_mul(z, nmod_mat_entry(basis, j, m*n/(d*d)), modp);
                    nmod_mat_entry(transformed_basis, k, j) =
                        nmod_add(nmod_mat_entry(transformed_basis, k, j), z, modp);
                    //transformed_basis(j, m) += chi_values[d % level] * pow((double)d, weight - 1) * basis(j, m*n/(d*d));
                }
            }
        }
        int x = nmod_mat_solve(Tn, smallbasis, transformed_basis);
        //cout << x << endl;
        //nmod_mat_print_pretty(smallbasis);
        //nmod_mat_print_pretty(transformed_basis);
        cout << n << " ";
        //nmod_mat_print_pretty(Tn);
        for(int a = 0; a <= 2*n; a++) {
            for(int k = 0; k < dim; k++) {
                nmod_mat_entry(X, k, k) = a;
            }
            nmod_mat_sub(Y, Tn, X);
            if(nmod_mat_det(Y) == 0) cout << a << " ";
            if(a == 0) continue;
            a = p - a;
            for(int k = 0; k < dim; k++) {
                nmod_mat_entry(X, k, k) = a;
            }
            nmod_mat_sub(Y, Tn, X);
            if(nmod_mat_det(Y) == 0) cout << a - p << " ";
            a = p - a;
        }
        n++;
        cout << endl;
    }

    nmod_mat_clear(Tn);
    nmod_mat_clear(basis);
    nmod_mat_clear(smallbasis);
    nmod_mat_clear(transformed_basis);
}

void cuspforms_modp::hecke_matrix(nmod_mat_t Tp, int l) {
    int dim = new_dimension();
    vector<int> basis_data = newspace_basis_data(false);
    if(verbose) {
        cout << "using basis rows";
        for(int k : basis_data) {
            cout << " " << k;
        }
        cout << endl;
    }
    compute_traces(l * basis_data[dim - 1] * basis_data[dim - 1] + 5);

    nmod_mat_t basis;
    nmod_mat_init(basis, dim, dim, p);
    for(int k = 0; k < dim; k++) {
        for(int j = k; j < dim; j++) {
            int n = basis_data[k];
            int m = basis_data[j];
            long z = trace_TnTm(n, m);
            nmod_mat_entry(basis, j, k) = z;
            nmod_mat_entry(basis, k, j) = z;
        }
    }

    nmod_mat_t Tp_basis;
    nmod_mat_init(Tp_basis, dim, dim, p);

    for(int k = 0; k < dim; k++) {
        for(int j = k; j < dim; j++) {
            int n = basis_data[k];
            int m = basis_data[j];
            long z1 = trace_TpTnTm(l, n, m);
            nmod_mat_entry(Tp_basis, j, k) = z1;
            nmod_mat_entry(Tp_basis, k, j) = z1;
        }
    }

    nmod_mat_init(Tp, dim, dim, p);
    nmod_mat_solve(Tp, basis, Tp_basis);

    nmod_mat_clear(basis);
    nmod_mat_clear(Tp_basis);
}




void cuspforms_modp::newspace_basis(nmod_mat_t basis, int ncoeffs) {
    newspace_basis_data();
    int d = new_dimension();
    if(d == 0) {
        nmod_mat_init(basis, 0, 0, p);
        return;
    }

    compute_traces(ncoeffs * basis_rows[d-1]);
    nmod_mat_init(basis, d, ncoeffs, p);
    for(int k = 0; k < d; k++) {
        int n = basis_rows[k];
        for(int m = 0; m < ncoeffs; m++) {
            nmod_mat_entry(basis, k, m) = trace_TnTm(n, m);
        }
    }
}

const vector<int>& cuspforms_modp::newspace_basis_data(bool coprime_only) {
    if(basis_rows.size() > 0 && coprime_only) return basis_rows;
    if(basis_rows2.size() > 0 && !coprime_only) return basis_rows2;
    int d = new_dimension();
    if(d == 0) return basis_rows;
    nmod_mat_t basis;
    nmod_mat_init(basis, d, d, p);
    compute_traces(d*d);
    int n = 1;
    int row = 0;
    vector<int> *basis_rows_ptr;
    if(coprime_only) {
        basis_rows_ptr = &basis_rows;
    }
    else {
        basis_rows_ptr = &basis_rows2;
    }
    while(row < d) {
        if(coprime_only) {while(GCD(n, level) > 1) n++;}
        basis_rows_ptr->push_back(n);
        int m = 1;
        int col = 0;
        while(col < d) {
            if(coprime_only) {while(GCD(m, level) > 1) n++;}
            //for(int m = 1; m < d + 1; m++) {
            nmod_mat_entry(basis, row, col) = trace_TnTm(n, m);
            col++;
        }
        row++;
        n++;
    }

    if(nmod_mat_det(basis) != 0) {
        nmod_mat_clear(basis);
        return basis_rows;
    }

    int k = 0;
    n = 1;
    basis_rows_ptr->clear();

    while(k < d) {
        if(coprime_only) {while(GCD(n, level) != 1) n++;}
        for(int j = 0; j < k; j++) {
            //while(GCD(m, level) != 1) m++;
            int m = basis_rows_ptr->at(j);
            long t = trace_TnTm(n,m);
            nmod_mat_entry(basis, k, j) = t;
            nmod_mat_entry(basis, j, k) = t;
            m++;
        }
        nmod_mat_entry(basis, k, k) = trace_TnTm(n, n);
        nmod_mat_t topleft;
        nmod_mat_window_init(topleft, basis, 0, 0, k + 1, k + 1);
        //cout << endl;
        //for(int i = 0; i < nmod_mat_nrows(topleft); i++) {
        //    for(int j = 0; j < nmod_mat_ncols(topleft); j++) {
        //        cout << nmod_mat_entry(topleft, i, j) << " ";
        //    }
        //    cout << endl;
        //}
        //cout << endl;
        if(nmod_mat_det(topleft) != 0) {
            k++;
            basis_rows_ptr->push_back(n);
        }
        //cout << nmod_mat_rank(topleft) << " " << k << " " << n << endl;
        n++;
        nmod_mat_clear(topleft);
    }
    return *basis_rows_ptr;
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
    //cout << A1 << endl;
    A1 = nmod_mul(A1, weight - 1, modp);
    //cout << A1 << endl;
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
    //cout << "After A1: trace(T1) = " << traces[1] << endl;
    //cout << "After A1: trace(T2) = " << traces[2] << endl;
    // computation of A2
    //cout << "computing A2 for level " << level << endl;

    long print_interval = 0;
    if(verbose > 1) { 
        cerr << "level " << level << endl;
        cerr << "start = " << start << " end = " << end << endl;
        cerr << "A2:"; print_interval = round(2 * sqrt(4*end)/70); print_interval = std::max(print_interval, 1l); }

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
    //cout << traces[1] << endl;
    //cout << "After A2: trace(T1) = " << traces[1] << endl;
    //cout << "After A2: trace(T2) = " << traces[2] << endl;

    if(verbose > 1) { cerr << "A3:"; print_interval = (long)std::max(round(sqrt(end)/70), 1.0); }
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

                        //cout << "y = " << y << " " << d << " " << n/d << " " << c << " " << level/c << endl;
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
    //cout << traces[1] << endl;
    //cout << "After A3: trace(T1) = " << traces[1] << endl;
    //cout << "After A3: trace(T2) = " << traces[2] << endl;

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
                    traces[n] = nmod_add(traces[n], t, modp);
                    //traces[n] += t;
                    //traces[n] %= p;
                }
            }
        }
    }
    //cout << "After A4: trace(T1) = " << traces[1] << endl;
    //cout << "After A4: trace(T2) = " << traces[2] << endl;

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
