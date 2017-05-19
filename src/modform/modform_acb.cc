#include <iostream>
#include <map>
#include <tuple>
#include <vector>
#include <complex>

#include "modform_acb.h"
#include "arb_extras.h"

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

void acb_reasonable_sqrt(acb_t out, const acb_t in, slong prec) {
    if(arb_is_negative(acb_realref(in)) && arb_contains_zero(acb_imagref(in))) {
        acb_neg(out, in);
        acb_sqrt(out, out, prec);
        acb_mul_onei(out, out);
        if(arf_sgn(arb_midref(acb_imagref(in))) < 0) {
            acb_neg(out, out);
        }
    }
    else {
        acb_sqrt(out, in, prec);
    }
}

static complex<double> acb_get_z(const acb_t in) {
    double x = arf_get_d(arb_midref(acb_realref(in)), ARF_RND_NEAR);
    double y = arf_get_d(arb_midref(acb_imagref(in)), ARF_RND_NEAR);
    return complex<double>(x,y);
}
static void acb_set_z(acb_t out, complex<double> in) {
    arb_set_d(acb_realref(out), in.real());
    arb_set_d(acb_imagref(out), in.imag());
}
static double arb_get_d(arb_t in) {
    return arf_get_d(arb_midref(in), ARF_RND_NEAR);
}

//using Eigen::ComplexEigenSolver;

typedef tuple<int,int,int> space_desc_t;
static map<space_desc_t, cuspforms_acb*> cache;

cuspforms_acb * get_cuspforms_acb(DirichletCharacter &chi, int weight, int verbose) {
    if(verbose > 2) cout << "getting space cuspforms_acb(" << chi.parent->q << ", " << weight << ", " << chi.m << ")" << endl;
    auto result = cache.find( space_desc_t(chi.parent->q, chi.m, weight) );
    if(result != cache.end()) return result->second;

    cuspforms_acb * S = new cuspforms_acb(chi, weight, verbose);
    cache[space_desc_t(chi.parent->q, chi.m, weight)] = S;
    return S;
}

void clear_cuspform_cache() {
    for(auto S : cache) {
        delete S.second;
    }
    cache.clear();
}


static bool has_unique_entries(const complex<double> * zz, int length, double eps = 1e-7) {
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
        while(j < length && std::abs(real(z[j]) - real(prev)) < eps) {
            if(std::abs(z[j] - prev) < eps) return false;
            j++;
        }
        prev = z[k];
    }
    delete [] z;
    return true;
}

static bool has_unique_entries(const double * zz, int length, double eps = 1e-7) {
    double * z = new double[length];
    std::copy(zz, zz + length, z);
    std::sort(z,z + length);
    for(int k = 1; k < length; k++) {
        if(std::abs(z[k] - z[k-1]) < eps) return false;
    }
    delete [] z;
    return true;
}

static bool has_unique_entries(const arb_mat_t z) {
    int N = arb_mat_nrows(z) * arb_mat_ncols(z);
    if(N == 0) return true;
    arb_ptr x = new arb_struct[N];
    for(int k = 0; k < N; k++) {
        x[k] = z[0].entries[k];
    }
    std::sort(x,x + N,
        [](arb_struct x, arb_struct y) {
            return arf_cmp(arb_midref(&x), arb_midref(&y)) < 0;
        });
    bool unique = true;
    for(int k = 1; k < N; k++) {
        if(arb_overlaps(&x[k], &x[k-1])) {
            cout << "overlaps" << endl;
            unique = false;
            break;
        }

        // placeholder until things are implemented more reasonably...
        double a = arf_get_d(arb_midref(&x[k]), ARF_RND_NEAR);
        double b = arf_get_d(arb_midref(&x[k-1]), ARF_RND_NEAR);
        if(std::abs(a - b) < 1.0e-6) {
            unique = false;
            break;
        }
    }
    delete [] x;
    return unique;
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

void cuspforms_acb::evalpoly(fmpz_t out, fmpz_t t, fmpz_t n) {
    // translated from Ralph's code
    int k = weight - 2;
    if(k == 0) {fmpz_set_ui(out, 1); return;}
    if(k == 1) {fmpz_set(out, t); return;}

    fmpz_t val[2];
    fmpz_init(val[0]);
    fmpz_init(val[1]);

    fmpz_set_ui(val[0], 1);
    fmpz_set(val[1], t);

    fmpz_t x, y;
    fmpz_init(x);
    fmpz_init(y);
    for(int i = 2; i <= k; i++) {
        fmpz_mul(x, t, val[(i + 1) % 2]);
        fmpz_mul(y, n, val[i % 2]);
        fmpz_sub(val[i % 2], x, y);
    }
    fmpz_set(out, val[k % 2]);

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(val[0]);
    fmpz_clear(val[1]);
}

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


long cuspforms_acb::new_dimension() {
    if(traces_computed < 2) compute_traces(2);
    long a = arf_abs_bound_lt_2exp_si(arb_midref(acb_realref(traces[1])));
    if(a > 62) {
        cerr << "Something has gone wrong and we getting the wrong"
                "result for the dimension. Or (less likely) the dimension"
                "is too big to store as a long." << endl;
        exit(0);
    }
    fmpz_t z;
    fmpz_init(z);
    arb_get_unique_fmpz(z, acb_realref(traces[1]));
    long d = fmpz_get_si(z);
    fmpz_clear(z);
    return d;
    //return (int)(round)(real(trace(1)));
}

void cuspforms_acb::newspace_basis(acb_mat_t B, int ncoeffs) {
    newspace_basis_data();
    int d = new_dimension();

    compute_traces(ncoeffs * basis_cols[d-1]);
    acb_mat_init(B, ncoeffs, d);
    //cmatrix_t basis(d, ncoeffs);
    for(int k = 0; k < d; k++) {
        int n = basis_cols[k];
        for(int m = 0; m < ncoeffs; m++) {
            trace_TnTm(acb_mat_entry(B, m, k), n, m);
            //basis(k, m) = trace_TnTm(n, m);
        }
    }
}

void cuspforms_acb::newforms(acb_mat_t out, int ncoeffs) {
    int dim = new_dimension();
    if(verbose) {
        cout << "computing newform basis for space of dimension " << dim << endl;
    }
    if(dim == 0) {acb_mat_init(out, 0, 0);}
    acb_t z1;
    acb_init(z1);
    vector<int> basis_data = newspace_basis_data();
    if(verbose > 1) {
        cout << "basis n for newspace:";
        for(int n : basis_data) cout << " " << n;
        cout << endl;
    }
    int coefficients_needed_for_full_rank = basis_data[dim - 1];
    ncoeffs = std::max(ncoeffs, 2 * coefficients_needed_for_full_rank + 5);
    if(verbose) cout << "computing traces" << endl;
    compute_traces(ncoeffs * basis_data[dim-1]);
    if(verbose) cout << "finished initial computation of traces" << endl;

    int sqrt_lim = std::max(1000, basis_data[dim - 1] + 1); // XXX. This is a potential problem.
                                                            // There is no reason that p in the
                                                            // !found_unique_eigenvalues loop could
                                                            // not get much bigger. (But N will have
                                                            // go get really large for that to happen,
                                                            // probably.
    acb_ptr invsqrts = _acb_vec_init(sqrt_lim);
    acb_set_ui(&invsqrts[1], 1u);
    for(int p = 2; p < sqrt_lim; p = next_prime(p)) {
        acb_reasonable_sqrt(&invsqrts[p], chi_values[p % level], prec);
        acb_conj(&invsqrts[p], &invsqrts[p]);
        for(int k = p; k < sqrt_lim; k += p) {
            acb_mul(&invsqrts[k], &invsqrts[k/p], &invsqrts[p], prec);
        }
    }


    arb_mat_t basis; // Will actually be a subset of the fourier coefficients
                     // of a basis. This will be a full rank symmetric matrix
                     // with entry (i,j) equal to T_{n_i, n_j}, where
                     // n_i is the ith entry in basis_data

    arb_mat_init(basis, dim, dim);
    for(int k = 0; k < dim; k++) {
        for(int j = k; j < dim; j++) {
            int n = basis_data[k];
            int m = basis_data[j];
            trace_TnTm(z1, n, m);
            acb_mul(z1, z1, &invsqrts[n], prec);
            acb_mul(z1, z1, &invsqrts[m], prec);
            arb_set(arb_mat_entry(basis, j, k), acb_realref(z1));
            arb_set(arb_mat_entry(basis, k, j), acb_realref(z1));
        }
    }


    if(verbose) cout << "computed basis matrix" << endl;

    arb_mat_t Tp_basis; // The basis after being acted on by Tp (or a sum of Tp)
    arb_mat_init(Tp_basis, dim, dim);
    int p = 2;
    bool found_unique_eigenvalues = false;

    acb_mat_t basis_transformation;
    arb_mat_t eigenvectors;
    arb_mat_t eigenvalues;

    acb_mat_init(basis_transformation, dim, dim);
    arb_mat_init(eigenvalues, dim, 1);
    arb_mat_init(eigenvectors, dim, dim);

    //rmatrix_t rr_Tp_basis(dim, dim);
    while(!found_unique_eigenvalues) {
        while(GCD(p, level) > 1) {p++;}
        if(verbose) cout << "computing hecke action for p = " << p << endl;
        compute_traces((p+1) * (basis_data[dim-1] + 1)*(basis_data[dim-1] + 1) + 5);

        for(int k = 0; k < dim; k++) {
            for(int j = k; j < dim; j++) {
                int n = basis_data[k];
                int m = basis_data[j];
                trace_TpTnTm(z1, p, n, m);
                acb_mul(z1, z1, &invsqrts[n], prec);
                acb_mul(z1, z1, &invsqrts[m], prec);
                acb_mul(z1, z1, &invsqrts[p], prec);
                arb_add(arb_mat_entry(Tp_basis, j, k), arb_mat_entry(Tp_basis, j, k), acb_realref(z1), prec);
                if(j != k)
                    arb_add(arb_mat_entry(Tp_basis, k, j), arb_mat_entry(Tp_basis, k, j), acb_realref(z1), prec);
            }
        }

        if(verbose)
            arb_mat_printd(Tp_basis, 10);
        if(verbose)
            cout << endl;

        int result = arb_mat_generalized_eigenproblem_symmetric_positive_definite(eigenvalues, eigenvectors, Tp_basis, basis, prec);

        if(verbose > 1) {
            cout << "eigenvalues:" << endl;
            arb_mat_print_sage_float(eigenvalues);
            cout << endl;
            arb_mat_printd(eigenvalues, prec/3);
            cout << endl;
            cout << "eigenvectors:" << endl;
            arb_mat_print_sage_float(eigenvectors);
            cout << endl;
            arb_mat_printd(eigenvectors, 10);
            cout << endl;
        }

        if(result == 0) {
            for(int j = 0; j < dim; j++) {
                for(int k = 0; k < dim; k++) {
                    acb_mul_arb(    acb_mat_entry(basis_transformation, j, k),
                                    &invsqrts[basis_data[j]],
                                    arb_mat_entry(eigenvectors, j, k),
                                    prec);
                }
            }
            found_unique_eigenvalues = true;
        }
        else {
            p++;
            if(verbose) {
                cout << "using more Hecke operators to get unique eigenvalues. next p = " << p << endl;
            }
        }

    }

    newspace_basis(out, 2);
    acb_mat_mul(out, out, basis_transformation, prec);
    for(int j = 0; j < dim; j++) {
        for(int k = 0; k < dim; k++) {
            acb_div(acb_mat_entry(basis_transformation, k,j), acb_mat_entry(basis_transformation,k,j), acb_mat_entry(out,1, j), prec);
        }
    }
    acb_mat_clear(out);
    newspace_basis(out, ncoeffs);
    acb_mat_mul(out, out, basis_transformation, prec);
}

const vector<int>& cuspforms_acb::newspace_basis_data() {
    if(basis_cols.size() > 0) return basis_cols;
    int d = new_dimension();
    if(d == 0) return basis_cols;
    basis_cols = modp_space->newspace_basis_data();
    return basis_cols;
}

void cuspforms_acb::compute_traces(int end) {
    // TODO
    if(verbose) cout << "cuspforms_acb: compute_traces called with end == " << end << endl;

    if(end < traces_computed) return;
    for(cuspforms_acb *subspace : subspaces) {
        subspace->compute_traces(end);
    }
    int start = traces_computed;
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
    if(verbose > 0) { cout << "A2:"; print_interval = round(2 * sqrt(4*end)/70); print_interval = std::max(print_interval, 1l); }

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
    fmpz_t polyterm;
    fmpz_init(polyterm);
    for(t = -t; (long)t*t < 4l*end; t++) {
        if(verbose > 0 && t % print_interval == 0)  {
            cout << '.';
            cout.flush();
        }

        for(long x = 0; x < level; x++) {
            //complex<double> chi_value = chi_values[x];
            //if(chi_value == 0.0) continue;
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
                        //complex<double> z = (complex<double>)(psi_table[level]/psi_table[level/g]);
                        acb_mul_arb(s1, chi_values[x], t1, prec);
                        acb_sub(traces[n], traces[n], s1, prec);
                        //z *= chi_value;
                        //z *= classnumbers[D/(f*f)];
                        //z /= 2.0;
                        //z *= polyterm;
                        //traces[n] -= z;
                        //z = nmod_mul(z, one_over_two, modp);
                        //traces[n] = nmod_sub(traces[n], z, modp);
                    }
                    k++;
                    f = square_divisors[k];
                } while(f != 1);
            }
        }
    }
    fmpz_clear(polyterm);
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
