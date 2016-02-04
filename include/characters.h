#ifndef __CHARACTERS_H__
#define __CHARACTERS_H__

#include "slint.h"

#include <complex>
#include <cmath>
#include <vector>
#include <set>
#include <fftw3.h>

#ifdef USE_MPFI
#include "mpfi_fft.h"
#endif

inline std::complex<double> e(double z) {
    std::complex<double>twopii(0,2*M_PI);
    return std::exp(twopii * z);
}

class DirichletGroup;


class DirichletCharacter {
public:
    long m;         // the label
    DirichletGroup * parent;

    DirichletCharacter() {}
    DirichletCharacter(DirichletGroup * parent_, long m_);

    ~DirichletCharacter() {}
    long exponent(long n);
    std::complex<double> max(long* index);
    std::complex<double> sum(long end);
    std::complex<double> value(long n);
    bool is_primitive();
    bool is_primitive_at_two();
    bool is_primitive_at_odd_part();
    bool is_even();
    void primitive_part_at_two(long * index, long * conductor);
    void primitive_part_at_known_p(long * index, long * conductor, long j);
    long conductor(long * index = NULL);
    std::complex<double> gauss_sum();
    void values_mod_p(long &p, long * chi_values);
    std::set<long> galois_orbit();
    std::complex<double> L_one();
};


class DirichletGroup {
    //
    //
public:
    long q;         // the modulus
    long q_odd;     // the odd part of the modulus
    long q_even;

    int k;          // the number of odd prime factors of q
 
    std::vector<long> * primes;
    std::vector<int> * exponents;
    std::vector<long> * generators;

    long ** A;      // exponent vectors:
                    //      for each m coprime to q we store an array
                    //      with the property that
                    //
                    //          m == g[j]**A[m][k] mod p[j]**e[j]
                    //
                    //      where g[j] is a primitive root mod p[j]**e[j],
                    //      and p[j] is the j-th prime factor of q.
                    //      (We don't actually store g[k], p[j], or e[j] anywhere.)
    
    long * B;

    long * PHI;     // PHI[j] = phi(q/phi(p[j]**e[j])). This will make it easier
                    // to compute the characters.

    long phi_q_odd;
    long phi_q;     // phi(q)


    std::complex<double> * zeta_powers_odd;     // zeta_powers[n] == e(n/phi(q))
    std::complex<double> * zeta_powers_even;     // zeta_powers[n] == e(n/phi(q))
    
    long index_from_primitive_character(long q0, long m);

#ifdef USE_MPFI
    mpfi_c_t * zeta_powers_odd_mpfi;
    mpfi_c_t * zeta_powers_even_mpfi;

    mpfi_c_t twopii_on_qodd;
    mpfi_c_t twopii_on_qeven;

    void DFTsum(mpfi_c_t * out, mpfi_c_t * in);
#endif

    void DFTsum(std::complex<double> * out, std::complex<double> * in);

    bool is_coprime_to_q(long n) {
        if(q_even > 1 && n % 2 == 0)
            return false;
        else if(q_odd == 1)
            return true;
        else
            return A[n % q_odd][0] != -1;
    }

    DirichletGroup() {
        q = 0;
    }
    DirichletGroup(long q_) : q(q_) {
        if(q == 0) return;
        q_even = 1;
        q_odd = q;
        while(q_odd % 2 == 0) {
            q_odd /= 2;
            q_even *= 2;
        }

        k = 0;
        phi_q = euler_phi(q);
        dft_translation_built = false;

        if(q_odd > 1) {
            primes = new std::vector<long>();
            exponents = new std::vector<int>();
            generators = new std::vector<long>();

            factors(q_odd, primes, exponents);
            k = primes->size();
            phi_q_odd = euler_phi(q_odd);
            
            PHI = new long[k];
            A = new long*[q_odd];
            for(int n = 0; n < q_odd; n++) {
                A[n] = new long[k];
            }
            zeta_powers_odd = new std::complex<double>[phi_q_odd];

#ifdef USE_MPFI
            zeta_powers_odd_mpfi = new mpfi_c_t[phi_q_odd];
#endif

            for(long j = 0; j < k; j++) {
                long x = pow(primes->at(j), exponents->at(j));
                long g = primitive_root(x);
                generators->push_back(g);
                long phi = (pow(primes->at(j), exponents->at(j) - 1) * (primes->at(j) - 1));
                PHI[j] = phi_q_odd/phi;
                long a = 1;
                for(long l = 0; l < phi; l++) {
                    for(long m = a; m < q_odd; m += x) {
                        A[m][j] = l;
                    }
                    a = (a * g) % x;
                }
            }
        
            //
            // for each m, 0 <= m < q, if (m,q) > 1, store
            // this as a flag in A[m][0]
            //
            for(long m = 0; m < q_odd; m++) {
                if(GCD(m,q_odd) > 1) {
                    A[m][0] = -1;
                }
            }

#ifdef USE_MPFI
            mpfi_c_init(twopii_on_qodd);
            mpfi_c_zero(twopii_on_qodd);
            mpfi_const_pi(twopii_on_qodd[0].im);
            mpfi_mul_ui(twopii_on_qodd[0].im, twopii_on_qodd[0].im, 2ul);
            mpfi_div_ui(twopii_on_qodd[0].im, twopii_on_qodd[0].im, phi_q_odd);
#endif

            for(unsigned long n = 0; n < phi_q_odd; n++) {
                zeta_powers_odd[n] = e(n/(double)phi_q_odd);

#ifdef USE_MPFI
                mpfi_c_init(zeta_powers_odd_mpfi[n]);
                mpfi_c_mul_ui(zeta_powers_odd_mpfi[n], twopii_on_qodd, n);
                mpfi_c_exp(zeta_powers_odd_mpfi[n], zeta_powers_odd_mpfi[n]);
#endif

            }
        } // end of initialization of everything having to do
          // with q_odd
        else {
            phi_q_odd = 1;
        }
        
        if(q_even > 4) {
            B = new long[q_even];
            zeta_powers_even = new std::complex<double>[q_even/4];

#ifdef USE_MPFI
            zeta_powers_even_mpfi = new mpfi_c_t[q_even/4];

            mpfi_c_init(twopii_on_qeven);
            mpfi_c_zero(twopii_on_qeven);
            mpfi_const_pi(twopii_on_qeven[0].im);
            mpfi_mul_ui(twopii_on_qeven[0].im, twopii_on_qeven[0].im, 2ul);
            mpfi_div_ui(twopii_on_qeven[0].im, twopii_on_qeven[0].im, q_even/4);
#endif

            for(unsigned long n = 0; n < q_even/4; n++) {
                zeta_powers_even[n] = e(4*n/double(q_even));

#ifdef USE_MPFI
                mpfi_c_init(zeta_powers_even_mpfi[n]);
                mpfi_c_mul_ui(zeta_powers_even_mpfi[n], twopii_on_qeven, n);
                mpfi_c_exp(zeta_powers_even_mpfi[n], zeta_powers_even_mpfi[n]);
#endif

            }
            long pow_five = 1;
            for(long e = 0; e < q_even/4; e++) {
                B[pow_five] = e;
                B[pow_five - 1] = 1;
                B[q_even - pow_five] = e;
                B[q_even - pow_five - 1] = -1;
                pow_five *= 5;
                pow_five %= q_even;
            }
        }

    }

    ~DirichletGroup() {
        if(q == 0) return;
        if(q_odd > 1) {
            delete [] zeta_powers_odd;
            for(int n = 0; n < q_odd; n++) {
                delete [] A[n];
            }

#ifdef USE_MPFI
            for(int n = 0; n < phi_q_odd; n++) {
                mpfi_c_clear(zeta_powers_odd_mpfi[n]);
            }
            delete [] zeta_powers_odd_mpfi;
            mpfi_c_clear(twopii_on_qodd);
#endif

            delete [] A;
            delete [] PHI;
            delete primes;
            delete exponents;
            delete generators;
        }
        if(q_even > 4) {
            delete [] B;
            delete [] zeta_powers_even;

#ifdef USE_MPFI
            for(int n = 0; n < q_even/4; n++) {
                mpfi_c_clear(zeta_powers_even_mpfi[n]);
            }
            delete [] zeta_powers_even_mpfi;
            mpfi_c_clear(twopii_on_qeven);
#endif
        }

        if(dft_translation_built) {
            delete [] dft_translation;
            delete [] idft_translation;
            delete [] dft_lengths;
        }
    }

    long chi_odd_exponent(long m, long n) {
        long x = 0;
        for(int j = 0; j < k; j++) {
            x += A[m][j]*A[n][j]*PHI[j];
            //x = x % phi_q_odd;
            if(x > 4294967296)
                x = x % phi_q_odd;
        }
        if(x >= phi_q_odd) 
            x = x % phi_q_odd;
        
        return x;
    }
    
    long chi_even_exponent(long m, long n) {
        long x = B[m]*B[n];
        if(B[m-1] == -1 && B[n-1] == -1)
            x += q_even/8;
        return x % (q_even/4);
    }

    std::complex<double> chi(long m, long n) {
        std::complex<double> even_part = 1.0;
        std::complex<double> odd_part = 1.0;
        if(q_even > 1) {
            if(m % 2 == 0 || n % 2 == 0) {
                return 0;
            }
            else if(q_even == 2) even_part = 1.0;
            else if(q_even == 4) {
                if(m % 4 == 3 && n % 4 == 3) even_part = -1.0;
                else even_part = 1.0;
            }
            else {
                even_part = zeta_powers_even[chi_even_exponent(m % q_even, n % q_even)];
            }
        }
        if(m >= q_odd)
            m %= q_odd;
        if(n >= q_odd);
            n %= q_odd;
        if(q_odd > 1) {
            if(A[m][0] == -1 || A[n][0] == -1)
                return 0;
            else
                odd_part = zeta_powers_odd[chi_odd_exponent(m, n)];
            if(q_even == 1)
                return odd_part;
        }
        return odd_part * even_part;
    }

    long exponent(long m, long n) {
        // return the a such that chi(m,n) == e(a/phi(q))
        long exponent;
        long odd_exponent = 0;
        if(q_odd > 1) {
            odd_exponent = chi_odd_exponent(m % q_odd, n % q_odd);
        }
        long even_exponent = 0;
        if(q_even > 4) {
            even_exponent = chi_even_exponent(m % q_even, n % q_even);
            even_exponent *= 2;
                                // the function just above computes the exponent of
                                // e(1/ (q_even/4) ), but we want the exponent of
                                // e(1/phi(q_even)) = e(1/(q_even/2))
        }
        else if(q_even == 4) {
            if(m % q_even == 3 && n % q_even == 3) {
                even_exponent = 1;
            }
            else {
                even_exponent = 0;
            }
        }
        
        if(q_even == 1) { // special case because phi(1) != 1/2.
            exponent = odd_exponent;
        }
        else {
            exponent = odd_exponent * q_even/2 + even_exponent * phi_q_odd;
        }
    
        // we now have the value of chi(m) as e(exponent/phi(q))

        // it could be equal to phi(q), though, and in that case we
        // want it to be zero...
        if(exponent == phi_q) {
            exponent -= phi_q;
        }
        return exponent;
    }

#ifdef USE_MPFI
    void chi(mpfi_c_t out, long m, long n) {
        std::complex<double> odd_part = 1.0;
        if(q_even > 1) {
            if(m % 2 == 0 || n % 2 == 0) {
                mpfi_c_zero(out);
                return;
            }
            else if(q_even == 2) {
                mpfi_c_set_ui(out, 1, 0);
            }
            else if(q_even == 4) {
                if(m % 4 == 3 && n % 4 == 3) {
                    mpfi_set_si(out[0].re, -1);
                    mpfi_set_si(out[0].im, 0);
                }
                else {
                    mpfi_c_set_ui(out, 1, 0);
                }
            }
            else {
                    mpfi_c_set(out, zeta_powers_even_mpfi[chi_even_exponent(m % q_even, n % q_even)]);
            }
        }
        else mpfi_c_set_ui(out, 1, 0);

        if(m >= q_odd)
            m %= q_odd;
        if(n >= q_odd);
            n %= q_odd;
        if(q_odd > 1) {
            if(A[m][0] == -1 || A[n][0] == -1)
                mpfi_c_zero(out);
            else
                mpfi_c_mul(out, out, zeta_powers_odd_mpfi[chi_odd_exponent(m, n)]);
        }
    }
#endif

    DirichletCharacter character(long m) {
        return DirichletCharacter(this, m);
    }

    void DFTsum_direct (std::complex<double> * out, std::complex<double> * in) {
        //
        // Set out[n] to the sum
        //
        //     sum_{k=1}^{q-1} in[k] chi(n,k)
        //
        // (computed naively, for testing purposes, or because
        //  something better hasn't been implemented yet.)

        
        for(int n = 0; n < q; n++) {
            std::complex<double> S = 0.0;
            for(int k = 0; k < q; k++) {
                S += in[k] * chi(n,k);
            }
            out[n] = S;
        }
    }

#ifdef USE_MPFI
    void DFTsum_direct (mpfi_c_t * out, unsigned long * in) {
        mpfi_c_t x, y;
        mpfi_c_init(x);
        mpfi_c_init(y);

        for(int n = 0; n < q; n++) {
            mpfi_c_zero(out[n]);
            if(!is_coprime_to_q(n)) continue;
            for(int k = 0; k < q; k++) {
                if(!is_coprime_to_q(k)) continue;
                chi(y, n, k);
                mpfi_c_mul_ui(x, y, in[k]);
                mpfi_c_add(out[n], out[n], x);
            }
        }

        mpfi_c_clear(x);
        mpfi_c_clear(y);
    }
#endif

    void all_sums(std::complex<double> * out, long end) {
        std::complex<double> * dft_in = new std::complex<double>[q]();
        for(long k = 0; k <= end; k++) {
            dft_in[k] = 1.0;
        }
        DFTsum(out, dft_in);
        delete [] dft_in;
    }

#ifdef USE_MPFI
    void all_sums(mpfi_c_t * out, long end) {
        mpfi_c_t * dft_in = new mpfi_c_t[q];
        for(long k = 0; k < q; k++) {
            mpfi_c_init(dft_in[k]);
            if(k <= end) {
                mpfi_c_set_ui(dft_in[k], 1ul, 0ul);
            }
            else {
                mpfi_c_zero(dft_in[k]);
            }
        }
        DFTsum(out, dft_in);
        for(int k = 0; k < q; k++) {
            mpfi_c_clear(dft_in[k]);
        }
        delete [] dft_in;
    }
#endif



    int * dft_translation;
    int * idft_translation;
    int * dft_lengths;
    int dft_dimension;
    bool dft_translation_built;

    void build_dft_translation(bool force=false) {
        if(dft_translation_built && !force) return;
        if(dft_translation_built && force) {
            delete [] dft_translation;
            delete [] idft_translation;
            delete [] dft_lengths;
        }

        int even_dimension;
        if(q_even <= 2)      even_dimension = 0;
        else if(q_even == 4) even_dimension = 1;
        else                 even_dimension = 2;

        dft_dimension = k + even_dimension;
        dft_lengths = new int[dft_dimension];
        
        if(even_dimension >= 1)
            dft_lengths[0] = 2;
        if(even_dimension == 2)
            dft_lengths[1] = q_even/4;

        for(int n = 0; n < k; n++) {
            long p = primes->at(n);
            long e = exponents->at(n);
            dft_lengths[n+even_dimension] = (p - 1)*pow(p, e-1);
        }

        dft_translation = new int[q];
        idft_translation = new int[q];
        int dlog;
        //          m == g[j]**A[m][k] mod p[j]**e[j]
        for(int n = 0; n < q; n++) {
            if((q_even > 1 && n % 2 == 0) || (q_odd > 1 && A[n % q_odd][0] == -1)) {
                dft_translation[n] = -1;
                idft_translation[n] = -1;
            }
            else {
                int dft_index = 0;
                int idft_index = 0;
                if(even_dimension >= 1) {
                    idft_index += (n % 4 == 3);
                    dft_index += (n % 4 == 3);
                    if(dft_dimension > 1) {
                        idft_index *= dft_lengths[1];
                        dft_index *= dft_lengths[1];
                    }
                }
                if(even_dimension == 2) {
                    dlog = B[n % q_even];
                    idft_index += dlog;
                    if(dlog > 0)
                        dft_index += (q_even/4 - B[n % q_even]);
                    if(dft_dimension > 2) {
                        idft_index *= dft_lengths[2];
                        dft_index *= dft_lengths[2];
                    }
                }
                for(int j = 0; j < k; j++) {
                    dlog = A[n % q_odd][j];
                    idft_index += dlog;
                    if(dlog > 0)
                        dft_index += (dft_lengths[j + even_dimension] - A[n % q_odd][j]);
                    if(j < k-1) {
                        idft_index *= dft_lengths[j + even_dimension + 1];
                        dft_index *= dft_lengths[j + even_dimension + 1];
                    }
                }
                dft_translation[n] = dft_index;
                idft_translation[n] = idft_index;
            }
        }

        dft_translation_built = true;
    }

    void print_dft_translation() {
        build_dft_translation();
        for(int n = 0; n < q; n++) {
            if(idft_translation[n] != -1) {
                std::cout << n << " " << idft_translation[n] << std::endl;
            }
        }
    }

    std::set< std::set<long> > galois_orbits() {
        std::set< std::set<long> > orbits;
        std::set<long> not_in_orbit_yet;
        for(long n = 1; n < q; n++) {
            if(GCD(n, q) == 1) not_in_orbit_yet.insert(n);
        }
        while(not_in_orbit_yet.size() > 0) {
            long next = *(not_in_orbit_yet.begin());
            auto o = character(next).galois_orbit();
            orbits.insert(o);
            for(long n : o) {
                not_in_orbit_yet.erase(n);
            }
        }
        return orbits;
    }
};

inline long DirichletGroup::index_from_primitive_character(long q0, long primitive_index) {
    //
    // Given a conductor q0 dividing q and an index of a character
    // mod q0, return the index of the character mod q which is the character
    // induced by chi_q0(m, -)
    //
    if(q0 == 1) return 1;
    if(primitive_index == 1) return 1;
    if(q % q0 != 0) return -1;

    std::vector<long> moduli;            // We're going to do this for each prime power
    std::vector<long> local_indices;     // dividing the modulus, and then CRT the results
                                    // together.
    for(int j = 0; j < k; j++) {
        long p = (*primes)[j];
        long modulus_at_p = ipow(p, (*exponents)[j]);
        moduli.push_back(modulus_at_p);

        long conductor_at_p = 1;
        long qq = q0;
        while(qq % p == 0) {
            conductor_at_p *= p;
            qq = qq/p;
        }
        long primitive_index_at_p = primitive_index % conductor_at_p;
        long index_at_p;
        if(conductor_at_p == 1) index_at_p = 1;
        else {
            index_at_p = PowerMod(primitive_index_at_p, modulus_at_p/conductor_at_p, modulus_at_p);
        }
        local_indices.push_back(index_at_p);
    }

    // now we need to deal with two
    // At 2, we can do almost the same thing, except that we also need to deal
    // with +- 1.

    long conductor_at_two = 1;
    long modulus_at_two = q_even;
    long qq = q0;
    while(qq % 2 == 0) {
        conductor_at_two *= 2;
        qq = qq/2;
    }
    long primitive_index_at_two = primitive_index % conductor_at_two;
    long index_at_two;
    if(conductor_at_two == 1) index_at_two = 1;
    else {
        index_at_two = PowerMod(primitive_index_at_two, modulus_at_two/conductor_at_two, modulus_at_two);
        // This is right up to a sign of +- 1. We need to adjust to make sure
        // that the indices are the same mod 4

        if(index_at_two % 4 != primitive_index_at_two % 4) {
            index_at_two = modulus_at_two - index_at_two;
        }
    }
    moduli.push_back(modulus_at_two);
    local_indices.push_back(index_at_two);
    return CRT(local_indices, moduli);
}

inline void DirichletGroup::DFTsum (std::complex<double> * out, std::complex<double> * in) {
    if(q < 4) {
        // just don't want to deal with problems with 2 right now.
        DFTsum_direct(out, in);
        return;
    }

    build_dft_translation();

    std::complex<double> *a, *X;
    a = new std::complex<double>[phi_q];
    X = new std::complex<double>[phi_q];


    for(int n = 0; n < q; n++) {
        if(dft_translation[n] == -1) continue;
        else a[dft_translation[n]] = in[n];
    }

    fftw_plan plan = fftw_plan_dft(dft_dimension, dft_lengths,
                                        (fftw_complex *)a,
                                        (fftw_complex *)X,
                                        FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    for(int n = 0; n < q; n++) {
        if(idft_translation[n] == -1) out[n] = 0.0;
        else out[n] = X[idft_translation[n]];
    }

    fftw_destroy_plan(plan);
    delete [] a;
    delete [] X;
}

#ifdef USE_MPFI
inline void DirichletGroup::DFTsum (mpfi_c_t * out, mpfi_c_t * in) {
    //if(q < 4) {
    //    // just don't want to deal with problems with 2 right now.
    //    DFTsum_direct(out, in);
    //    return;
    //}

    build_dft_translation();

    mpfi_c_t *a;//, *X;
    a = new mpfi_c_t[phi_q];
    //X = new mpfi_c_t[phi_q];

    for(int n = 0; n < phi_q; n++) {
        mpfi_c_init(a[n]);
        //mpfi_c_init(X[n]);
    }


    for(int n = 0; n < q; n++) {
        // This might be inefficient.
        // Since the mpfi_c array is really just an array of pointers,
        // we might be able to just copy the pointers without copying
        // the data. The fft wouldn't be operating on continuous memory, then,
        // but it still might be better than this copy.
        if(dft_translation[n] == -1) continue;
        else mpfi_c_set(a[dft_translation[n]], in[n]);
    }

    unsigned long * lengths = new unsigned long[dft_dimension];
    for(int n = 0; n < dft_dimension; n++) lengths[n] = dft_lengths[n];
    ndft(a, phi_q, dft_dimension, lengths);
    //fftw_plan plan = fftw_plan_dft(dft_dimension, dft_lengths,
    //                                    (fftw_complex *)a,
    //                                    (fftw_complex *)X,
    //                                    FFTW_BACKWARD, FFTW_ESTIMATE);

    for(int n = 0; n < q; n++) {
        if(idft_translation[n] == -1) mpfi_c_zero(out[n]);
        else mpfi_c_set(out[n], a[idft_translation[n]]);
    }

    for(int n = 0; n < phi_q; n++) {
        mpfi_c_clear(a[n]);
    }

    delete [] a;

}
#endif



inline DirichletCharacter::DirichletCharacter(DirichletGroup * parent_, long m_) : parent(parent_), m(m_) {
        //
        // This doesn't actually do anything. Computing characters
        // using DirichletCharacter is not going to be any faster
        // than using DirichletGroup::chi(), at least for now.
        //
}

inline std::complex<double> DirichletCharacter::gauss_sum() {
    std::complex<double> S = 0;
    long q = parent->q;
    std::complex<double> z = e(1.0/q);
    std::complex<double> x = z;
    for(long n = 1; n < q; n++) {
        S = S + value(n) * x;
        x = x * z;
    }
    return S;
}

inline std::complex<double> DirichletCharacter::max(long * index) {
    // return the max of the partial sums of the character,
    // and set *index to the spot where the max occurs,
    // unless index = 0
    //
    // note that *index is always <= (q-1)/2
    
    //if(parent->A[m][0] == -1) {
    //    if(index != 0){
    //        *index = -1;
    //    }
    //    return 0;
    //}

    std::complex<double> S(0,0);
    double absmax = 0.0;
    std::complex<double> current_max = 0.0;
    long max_location;

    for(long n = 0; n <= (parent->q-1)/2; n++) {
        S = S + value(n);
        if(abs(S) > absmax) {
            absmax = abs(S);
            current_max = S;
            max_location = n;
        }
    }

    if(index != 0) {
        *index = max_location;
    }
    return current_max;
}

inline std::complex<double> DirichletCharacter::sum(long end) {
    std::complex<double> S(0,0);

    for(long n = 0; n <= end; n++) {
        S = S + value(n);
    }

    return S;
}

// Is this right...?
//
//long DirichletCharacter::exponent(long n) {
//    long x = 0;
//    for(int j = 0; j < parent->k; j++) {
//        x += parent->A[m][j]*parent->A[n][j]*parent->PHI[j];
//        if(x > 4294967296)
//            x = x % parent->phi_q;
//    }
//    if(x >= parent->phi_q) {
//        x = x % parent->phi_q;
//    }
//    return x;
//}

inline long DirichletCharacter::exponent(long n) {
    return parent->exponent(m,n);
}

inline std::complex<double> DirichletCharacter::value(long n) {
    return parent->chi(m,n);
    //if(parent->A[m][0] == -1 || parent->A[n][0] == -1)
    //    return 0;
    //return parent->zeta_powers[exponent(n)];
}

inline bool DirichletCharacter::is_primitive() {
    // return whether or not this character is primitive

    return is_primitive_at_two() && is_primitive_at_odd_part();
}

inline bool DirichletCharacter::is_primitive_at_odd_part() {
    // this is computed one prime at a time, and the
    // character will be primitive if and only if it
    // is primitive at every prime.
    if(parent->q_odd == 1) return true;
    else {
        long n = m % parent->q_odd;
        if(parent->A[n][0] == -1) return false;
        for(int j = 0; j < parent->k; j++) {
            if(parent->A[n][j] % parent->primes->at(j) == 0)
                return false;
        }
        return true;
    }
}

inline bool DirichletCharacter::is_primitive_at_two() {
    long q_even = parent->q_even;
    long * B = parent->B;
    long n = m % q_even;
    if(q_even == 1) return true;
    else if(q_even == 2) return false;
    else if(q_even == 4) return n == 3;
    else {
        n = n % 8;
        if(n == 3) return true;
        return n == 5;
    }
}

inline bool DirichletCharacter::is_even() {
    // return whether or not this character is even
    //
    // We just figure out if the character is even by evaluating it at q-1.
    // The evaluation isn't going to be exact, but since the number is just
    // +-1 we can just check which of these it is closest to.
    //
    
    return abs(value(parent->q - 1) - 1.0) < .5;
}

inline void DirichletCharacter::primitive_part_at_known_p(long * index, long * conductor, long j) {
    //
    //    Return the conductor and the index of the primitive character
    //    associated to the p-part of the character for the j-th odd prime
    //    factor p of the modulus.
    //

    long p = parent->primes->at(j);
    long mm = m % parent->q_odd;
    //long dlog = parent->A[mm * parent->k + j];
    long dlog = parent->A[mm][j];
    if(dlog == 0) {
        *conductor = 1; *index = 1; return;
    }
    int e = 0;
    while(dlog % p == 0) {
        dlog /= p;
        e = e + 1;
    }
    *conductor = PowerMod(p, parent->exponents->at(j) - e, 1l << 62); // this is silly, but we're
                                                                     // not going to deal with
                                                                     // conductors larger than this

    *index = PowerMod(parent->generators->at(j), dlog, *conductor);
}

inline void DirichletCharacter::primitive_part_at_two(long * index, long * conductor) {
    // Return the conductor and the index of the primitive
    // character associated to the even part of the modulus.

    // We return a pair (n, q), where q is the even part of the
    // conductor of chi and n in the index of the even part of the
    // primitive character inducing chi.

    long q_even = parent->q_even;
    long mm = m % q_even;

    // This is an annoying computation.

    if(q_even == 1) {
        *index = 1;
        *conductor = 1;
    }
    else if(q_even == 2) {
        *index = 1;
        *conductor = 1;
    }
    else if(q_even == 4) {              // When q_even is 4, the conductor
        if(mm == 1) {*index = 1; *conductor = 1;}
        if(mm == 3) {*index = 3; *conductor = 4;}
    }
    else if(q_even == 8) {
        if(mm == 1) {*index = 1; *conductor = 1;}
        if(mm == 3) {*index = 3; *conductor = 8;}
        if(mm == 5) {*index = 5; *conductor = 8;}
        if(mm == 7) {*index = 3; *conductor = 4;}
    }
    else {
        long alpha   =   parent->B[mm];
        long epsilon = parent->B[mm-1];
        if(alpha == 0) {
            if(epsilon == 1) {*index = 1; *conductor = 1;}
            else {*index = 3; *conductor = 4;}
            return;
        }
        int f = 0;
        while(alpha % 2 == 0) {
            alpha /= 2;
            f = f + 1;
        }
        *conductor = q_even >> f;
        *index = PowerMod(5, alpha, *conductor);
        if(epsilon == -1) {
            *index = *conductor - *index;
        }
    }
}

inline long DirichletCharacter::conductor(long * index) {
    std::vector<long> moduli;
    std::vector<long> indices;
    long q = parent->q;
    long q1;
    long m1;
    for(int j = 0; j < parent->k; j++) {
        primitive_part_at_known_p(&m1, &q1, j);
        moduli.push_back(q1);
        indices.push_back(m1);
    }

    if(q % 2 == 0) {
        primitive_part_at_two(&m1, &q1);
        moduli.push_back(q1);
        indices.push_back(m1);
    }

    long q0 = 1;
    for(auto q1 : moduli) {
        q0 *= q1;
    }

    if(index) {
        *index = CRT(indices, moduli);
    }
    return q0;
}

inline void DirichletCharacter::values_mod_p(long & p, long * chi_values) {
    //
    // Given a prime p such that the order of chi divides p - 1, fill the array
    // chi_values with values of this character mod p.
    //
    // If p == 0, we choose an appropriate prime and set p to that prime.
    //

    int q = parent->q;
    int order;
    if(m == 0)
        order = 1;
    else
        order = order_mod(m, q);

    int phi_q = parent->phi_q;
    if(p == 0) {
        p = order + 1;
        while(p < 10000) p += order;
    }
    else {
        if( (p - 1) % order != 0) {
            p += order - (p - 1) % order;
        }
    }
    while(!is_prime(p) || phi_q % p == 0) p += order;
            // can p divide phi(N)?...
            // not going to think about
            // that right now.
    
    if( (p - 1) % order != 0) {
        std::cerr << "values_mod_p called with an invalid prime."
             << "p = " << p << std::endl
             << "order = " << order << std::endl;

        exit(1);
    }
    long g = primitive_root(p);
    int exponent_adjustment = phi_q/order;
    g = PowerMod(g, (p-1)/order, p);

    for(int k = 0; k < q; k++) {
        if(GCD(k, q) != 1) {
            chi_values[k] = 0;
            continue;
        }
        long e = exponent(k);
        if(e % exponent_adjustment != 0) {
            std::cerr << q << " " << m << " " << e << " " << order << " " << phi_q << " " <<
                exponent_adjustment << std::endl;
            std::cerr << "error" << std::endl;
            exit(1);
        }
        chi_values[k] = PowerMod(g, e/exponent_adjustment, p);
    }

}

inline std::set<long> DirichletCharacter::galois_orbit() {
    std::set<long> orbit;
    long z = 1;
    for(int l = 0; l < parent->k; l++) {
        long p = parent->primes->at(l);
        long e = parent->exponents->at(l);
        z = LCM(z, (p - 1)*ipow(p, e - 1));
    }
    if(parent->q_even > 1) {
        z = LCM(z, parent->q_even/2);
    }

    for(int a = 1; a < z; a++) {
        if(GCD(a,z) == 1) {
            orbit.insert(PowerMod(m, a, parent->q));
        }
    }

    return orbit;
}

inline std::complex<double> DirichletCharacter::L_one() {
    // we don't do this very efficiently. But then
    // since we already compute a full table of discrete logs,
    // we don't do anything very efficiently.

    // (We can easily compute all the values of L(1, chi), just about
    // optimally, using an FFT with the function DFTsum()

    // WARNING: We assume that chi is primitive.

    long q = parent->q;

    if(is_even()) {
        std::complex<double> S = 0.0;
        for(long k = 1; k < q; k++) {
            S = S + value(k)*std::log(std::sin(M_PI * k/(double)q));
        }
        S = std::conj(S);
        S *= -gauss_sum()/(double)q;
        return S;
    }
    else {
        std::complex<double> S = 0.0;
        for(long k = 1; k <= q/2; k++) {
            S += value(k);
        }
        S *= (std::complex<double>(0, 1) * M_PI) / ((2.0 - value(2)) * gauss_sum());
        S = std::conj(S);
        return S;
    }
}
#endif
