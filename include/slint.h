//
// Simple Library for Number Theory
//

#ifndef _SLINT_H_
#define _SLINT_H_

#include <vector>
#include <iostream>

long next_prime(long);

// taken from NTL, stripped of
// overflow checking
void XGCD(long& d, long& s, long& t, long a, long b)
{
   long  u, v, u0, v0, u1, v1, u2, v2, q, r;

   long aneg = 0, bneg = 0;

   if (a < 0) {
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      b = -b;
      bneg = 1;
   }

   u1=1; v1=0;
   u2=0; v2=1;
   u = a; v = b;

   while (v != 0) {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 =  u1 - q*u2;
      v2 = v1- q*v2;
      u1 = u0;
      v1 = v0;
   }

   if (aneg)
      u1 = -u1;

   if (bneg)
      v1 = -v1;

   d = u;
   s = u1;
   t = v1;
}
   
// taken from NTL:
long InvMod(long a, long n)
{
   long d, s, t;

   XGCD(d, s, t, a, n);
   if (d != 1) return -1;
   if (s < 0)
      return s + n;
   else
      return s;
}

// taken from NTL:
static inline long MulMod(long a, long b, long n)
{

    long q, res;

    q  = (long) ((((double) a) * ((double) b)) / ((double) n)); 
    res = a*b - q*n;
    if (res >= n)
        res -= n;
    else if (res < 0)
        res += n;
    return res;
}

// taken from NTL
long PowerMod(long a, long ee, long n)
{
    long x, y;

    unsigned long e;

    if (ee < 0)
        e = - ((unsigned long) ee);
    else
        e = ee;

    x = 1;
    y = a;
    while (e) {
        if (e & 1) x = MulMod(x, y, n);
        y = MulMod(y, y, n);
        e = e >> 1;
    }

    if (ee < 0) x = InvMod(x, n);

    return x;
}

// taken from NTL,
// stripped of overflow checking
long GCD(long a, long b)
{
   long u, v, t, x;

   if (a < 0) {
      a = -a;
   }

   if (b < 0) {
      b = -b;
   }


   if (b==0)
      x = a;
   else {
      u = a;
      v = b;
      do {
         t = u % v;
         u = v; 
         v = t;
      } while (v != 0);

      x = u;
   }

   return x;
}

long euler_phi(long n) {
    // yes, this is stupidly slow...
    //
    // I don't care.
    //
    long phi = 1;
    long p = 2;
    long p_power = 1;
    while(n > 1) {
        p_power = 1;
        while( (n % p) == 0 ) {
            n = n / p;
            p_power *= p;
        }
        phi = phi * ( p_power - p_power/p ); // Note: if p_power == 1, then p_power/p == 0
        if(p == 2)
            p = 3;
        else
            p = p + 2;
    }
    return phi;
}

void factors(long n, std::vector<long> * primes, std::vector<int> * exponents) {
    //
    // appends the prime factors of n to *primes,
    // and if exponents if not NULL, appends the exponents
    // of those factor to *exponents
    //
    // yes, this is stupidly slow.
    //
    // i don't care...
    //
    long p = 2;
    int a = 0;
    while(n > 1) {
        a = 0;
        while( (n % p) == 0 ) {
            n = n / p;
            a++;
        }
        if(a != 0) {
            primes->push_back(p);
            if(exponents != NULL) {
                exponents->push_back(a);
            }
        }
        if(p == 2)
            p = 3;
        else
            p = p + 2;
    }
}

bool is_squarefree(long n) {
    //
    // factors n by trial division, checking for squarefreeness along the
    // way.
    //
    // yes, this is stupidly slow.
    //
    // i don't care...
    //
    long p = 2;
    int a = 0;
    while(n > 1) {
        a = 0;
        while( (n % p) == 0 ) {
            n = n / p;
            a++;
        }
        if(a > 1) return false;
        if(p == 2)
            p = 3;
        else
            p = next_prime(p);
            //p = p + 2;
    }
    return true;
}

bool is_fundamental_discriminant(long n) {
    if(n % 4 == 1) {
        return is_squarefree(n);
    }
    else if(n % 4 == 0) {
        long m = n/4;
        if(m % 4 == 2 || m % 4 == 3) {
            return is_squarefree(m);
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}

long primitive_root(long n) {
    //
    // Return a primitive root mod n.
    //
    // If n == 1 or 2, returns 1
    // If n == 4, returns 3
    // If n is an odd prime power p^e with p < 3037000499, returns
    // the smallest primitive root mod p which is a primitive
    // root for all e.
    // If n in an odd prime larger than 3037000499, returns the
    // smallest primitive root mod n
    //
    if(n < 2) {
        return n;
    }
    if(n == 2)
        return 1;
    if(n == 4)
        return 3;
    std::vector<long> prime_factors;
    factors(n, &prime_factors, NULL);
    if(prime_factors.size() > 1)
        return -1;

    long p = prime_factors[0];
    if(p == 2)
        return -1;
    long p2;
    if(p > 3037000499) {
        p2 = p; // when p is too large, we still compute
                // a primitive root, but we don't verify
                // that it is a primitive root for all powers of p
    }
    else {
        p2 = p*p;
    }
    long phi = p - 1;
    std::vector<long> phi_prime_factors;
    factors(phi, &phi_prime_factors, NULL);
    long a = 1;
    while(a < n) {
        a++;
        if(GCD(a,n) > 1)
            continue;
        bool root = true;
        for(    std::vector<long>::iterator i = phi_prime_factors.begin();
                i != phi_prime_factors.end();
                i++     ) {
            //std::cout << p << " " << *i << std::endl;
            if(PowerMod(a, phi/(*i), p) == 1) {
                root = false;
                break;
            }
        }
        if(root) {
            if(p == p2)
                return a;
            else {
                long x = PowerMod(a, p, p2);
                if(x != a)
                    return a;
            }
        }
    }
    return -1;
}

bool is_prime_power(long q) {
    std::vector<long> primes;
    factors(q, &primes, NULL);
    if(primes.size() == 1)
        return true;
    else
        return false;
}

bool MR_test(long n, long a) {
    long d = n - 1;
    int s = 0;
    while(d % 2 == 0) {
        d = d/2;
        s = s + 1;
    }
    long x = PowerMod(a, d, n);
    if(x == 1 || x == n-1) {
        return true;
    }
    int r = 1;
    while(r < s) {
        x = x*x % n;
        if(x == 1) return false;
        if(x == n-1) return true;
        r++;
    }
    return false;
}

bool is_prime(long q) {
    if(q == 2 || q == 7 || q == 61) return true;
    if(q < 4759123141l) {
        return MR_test(q,2) && MR_test(q,7) && MR_test(q,61);
    }
    else {
        std::cerr << "You should be using FLINT." << std::endl;
        std::vector<long> primes;
        std::vector<int> exponents;
        factors(q, &primes, &exponents);
        if(primes.size() == 1 && exponents[0] == 1)
            return true;
        else
            return false;
    }
}

long next_prime(long n) {
    if(n < 2)
        return 2;
    if(n == 2)
        return 3;
    if(n % 2 == 0)
        n += 1;
    else
        n += 2;
    while(!is_prime(n)) {
        n += 2;
    }
    return n;
}

long odd_part(long n) {
    if(n == 0) {
        return 1;
    }
    while(n % 2 == 0) {
        n = n/2;
    }
    return n;
}

long kronecker(long n, long m) {
    if(GCD(n,m) != 1) {
        return 0;
    }
    long t = 1;
    while(m > 1) {
        long m_odd, m_even;
        m_odd = odd_part(m);
        m_even = m/m_odd;
        if(n % 8 == 3 || n % 8 == 5) {
            while(m_even % 2 == 0) {
                t = -t;
                m_even /= 2;
            }
        }
        //if(m_even == 2) {
        //    if(n % 8 == 3 || n % 8 == 5) t = -t;
        //}
        n = n % m_odd;
        if(odd_part(n) % 4 == 3 && m_odd % 4 == 3) t = -t;
        long x = n;
        n = m_odd;
        m = x;
    }
    if(m == 0) {
        if(n != 1) return 0;
    }
    return t;
}



long kronecker2(long n, long m) {
    long t = 1;
    //cout << n << " " << m << " " << t << endl;
    long m_even, m_odd;
    m_odd = odd_part(m);
    m_even = m/m_odd;
    while(m > 2) {
        if(m_even == 2) {
            if(n % 8 == 3 || n % 8 == 5) t = -t;
        }
        n = n % m_odd;
        if(odd_part(n) % 4 == 3 && m_odd % 4 == 3) t = -t;
        //if(odd_part(m) % 4 == 3) t = -t;
        long x = n;
        n = m_odd;
        m = x;
        //cout << n << " " << m << " " << t << endl;
        m_odd = odd_part(m);
        m_even = m/m_odd;
    }
    if(m == 2) {
        if(n % 2 == 0) return 0;
        if(n % 8 == 3 || n % 8 == 5) t = -t;
    }
    else if (m == 0) {
        if(n != 1) return 0;
    }
    return t;
}

void prime_range(std::vector<long> * primes, long end, long start = 2) {
    // fill primes with a list of prime numbers between
    // start and end (including start but not end)

    if(start > 2) {
        std::cerr << "that's not implemented." << std::endl;
        return;
    }

    // a really simple sieve...
    bool * sieve_range = new bool[end]();
    sieve_range[0] = 1;
    sieve_range[1] = 1;
    long p = 2;
    while(p < end) {
        primes->push_back(p);
        for(long k = 2*p; k < end; k += p) {
            sieve_range[k] = 1;
        }
        do {
            p++;
        } while (p < end && sieve_range[p] == 1);
    }

    delete [] sieve_range;
}

void squarefree_range(std::vector<long> * squarefrees, long end, long start = 1) {
    // fill primes with a list of prime numbers between
    // start and end (including start but not end)

    if(start > 1) {
        std::cerr << "that's not implemented." << std::endl;
        return;
    }

    // a really simple sieve...
    bool * sieve_range = new bool[end]();
    sieve_range[0] = 1;
    sieve_range[1] = 0;
    squarefrees->push_back(1);
    long n = 2;
    while(n < end) {
        squarefrees->push_back(n);
        for(long k = n*n; k < end; k += n*n) {
            sieve_range[k] = 1;
        }
        do {
            n++;
        } while (n < end && sieve_range[n] == 1);
    }

    delete [] sieve_range;
}

void discriminant_range(std::vector<long> * discs, long end, long start = 1) {
    std::vector<long> sqfrees;
    squarefree_range(&sqfrees, end, start);

    int index1 = 0;
    int index2 = 0;
    int size = sqfrees.size();
    long n1 = sqfrees[index1];
    long n2 = sqfrees[index2];
    while(n2 < end + 1) {
        while(4*n2 < n1) {
            if(n2 % 4 == 2 || n2 % 4 == 3) discs->push_back(4*n2);
            index2++;
            n2 = sqfrees[index2];
        }
        if(n1 == end) break;
        if(n1 % 4 == 1) discs->push_back(n1);
        index1++;
        if(index1 > size) n1 = end;
        else n1 = sqfrees[index1];
    }
}

int order_mod(int n, int q) {
    //
    // return the smallest e such that n^e == 1 mod q, or
    // -1 if there is no such e
    //

    if(GCD(n,q) != 1) return -1;
    int z = n;
    int e = 1;
    while(z != 1) {
        z = (z * n) % q;
        e++;
    }
    return e;
}

long CRT(long a, long b, long m, long n) {
    //
    // return x == a mod m and b mod n
    //

    long minv = InvMod(m, n);
    long ninv = InvMod(n, m);
    long q = m*n;

    a = a % m;
    b = b % n;

    long x1 = MulMod(n, MulMod(a,ninv,q), q);
    long x2 = MulMod(m, MulMod(b,minv,q), q);
    
    return (x1 + x2) % q;
}

long CRT(std::vector<long> a, std::vector<long> n) {
    long x = a[0];
    long q = n[0];
    for(int k = 1; k < a.size(); k++) {
        x = CRT(x, a[k], q, n[k]);
        q = q*n[k];
    }
    if(q == 1) return 0;
    return x;
}

inline long ipow(long a, long n) {
    long z = 1;
    for(int j = 0; j < n; j++) {
        z = z * a;
    }
    return z;
}

std::vector<long> divisors(long n) {
    std::vector<long> x;
    for(int d = 1; d <= n; d++) {
        if(n % d == 0) x.push_back(d);
    }
    return x;
}


struct int_factor_t {
    int p;
    int e;
    int f; // f == p^e
};

int_factor_t * factor_table;
int factor_table_size = 0;

void build_factor_table(int size) {
    if(factor_table_size > 0) {
        delete [] factor_table;
    }
    if(size == 0) {
        factor_table_size = 0;
        return;
    }
    factor_table = new int_factor_t[size];
    for(int k = 0; k < size; k++) {
        factor_table[k] = {1,1,1};
    }
    std::vector<long> primes;
    prime_range(&primes, size);
    for(auto p : primes) {
        int e = 1;
        long f = p;
        while(f < size) {
            for(int n = f; n < size; n += f) {
                factor_table[n] = {(int)p, e, (int)f};
            }
            f = f * p;
            e++;
        }
    }
    factor_table_size = size;
}

struct int_factorization_t {
    // factorization of 32 bit (signed) integers
    int_factor_t factors[9];
    int n;
    int nfactors;
};

void factor(int n, int_factorization_t &factorization) {
    int sign = 1;
    if(n < 0) {sign = -1; n = -n;}
    if(n >= factor_table_size) {
        std::cerr << "factorization of integers larger than the precomputed size not supported yet." << std::endl;
    }
    factorization.n = n;
    factorization.nfactors = 0;
    while(n > 1) {
        int_factor_t f = factor_table[n];
        factorization.factors[factorization.nfactors] = f;
        n = n/f.f;
        factorization.nfactors++;
    }
}

#endif
