//
// Simple Library for Number Theory
//

#ifndef _SLINT_H_
#define _SLINT_H_

#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

// XXX lazy...
// 0 terminated array
const int prime_powers_table[] = {
    2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43,
    47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97, 101, 103, 107, 109,
    113, 121, 125, 127, 128, 131, 137, 139, 149, 151, 157, 163, 167, 169, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 243, 251,
    256, 257, 263, 269, 271, 277, 281, 283, 289, 293, 307, 311, 313, 317, 331,
    337, 343, 347, 349, 353, 359, 361, 367, 373, 379, 383, 389, 397, 401, 409,
    419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
    503, 509, 512, 521, 523, 529, 541, 547, 557, 563, 569, 571, 577, 587, 593,
    599, 601, 607, 613, 617, 619, 625, 631, 641, 643, 647, 653, 659, 661, 673,
    677, 683, 691, 701, 709, 719, 727, 729, 733, 739, 743, 751, 757, 761, 769,
    773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 841, 853, 857, 859, 863,
    877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 961, 967, 971,
    977, 983, 991, 997, 1009, 1013, 1019, 1021, 1024, 1031, 1033, 1039, 1049,
    1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123,
    1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
    1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301,
    1303, 1307, 1319, 1321, 1327, 1331, 1361, 1367, 1369, 1373, 1381, 1399,
    1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481,
    1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559,
    1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627,
    1637, 1657, 1663, 1667, 1669, 1681, 1693, 1697, 1699, 1709, 1721, 1723,
    1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823,
    1831, 1847, 1849, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907,
    1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003,
    2011, 2017, 2027, 2029, 2039, 2048, 2053, 2063, 2069, 2081, 2083, 2087,
    2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179,
    2187, 2197, 2203, 2207, 2209, 2213, 2221, 2237, 2239, 2243, 2251, 2267,
    2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347,
    2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2401, 2411, 2417,
    2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539,
    2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647,
    2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711,
    2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797,
    2801, 2803, 2809, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887,
    2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999,
    3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089,
    3109, 3119, 3121, 3125, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203,
    3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307,
    3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389,
    3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3481, 3491,
    3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
    3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671,
    3673, 3677, 3691, 3697, 3701, 3709, 3719, 3721, 3727, 3733, 3739, 3761,
    3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853,
    3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943,
    3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051,
    4057, 4073, 4079, 4091, 4093, 4096, 4099, 4111, 4127, 4129, 4133, 4139,
    4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243,
    4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349,
    4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457,
    4463, 4481, 4483, 4489, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549,
    4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651,
    4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759,
    4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877,
    4889, 4903, 4909, 4913, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967,
    4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5041,
    5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153,
    5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273,
    5279, 5281, 5297, 5303, 5309, 5323, 5329, 5333, 5347, 5351, 5381, 5387,
    5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471,
    5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563,
    5569, 5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659,
    5669, 5683, 5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779,
    5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
    5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981,
    5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089,
    6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199,
    6203, 6211, 6217, 6221, 6229, 6241, 6247, 6257, 6263, 6269, 6271, 6277,
    6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361,
    6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481,
    6491, 6521, 6529, 6547, 6551, 6553, 6561, 6563, 6569, 6571, 6577, 6581,
    6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701,
    6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803,
    6823, 6827, 6829, 6833, 6841, 6857, 6859, 6863, 6869, 6871, 6883, 6889,
    6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983,
    6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103,
    7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213,
    7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307, 7309, 7321, 7331,
    7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477,
    7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549,
    7559, 7561, 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643,
    7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741,
    7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873,
    7877, 7879, 7883, 7901, 7907, 7919, 7921, 7927, 7933, 7937, 7949, 7951,
    7963, 7993, 8009, 8011, 8017, 8039, 8053, 8059, 8069, 8081, 8087, 8089,
    8093, 8101, 8111, 8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8192,
    8209, 8219, 8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
    8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, 8389, 8419,
    8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, 8513, 8521, 8527, 8537,
    8539, 8543, 8563, 8573, 8581, 8597, 8599, 8609, 8623, 8627, 8629, 8641,
    8647, 8663, 8669, 8677, 8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731,
    8737, 8741, 8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
    8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, 8933, 8941,
    8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, 9013, 9029, 9041, 9043,
    9049, 9059, 9067, 9091, 9103, 9109, 9127, 9133, 9137, 9151, 9157, 9161,
    9173, 9181, 9187, 9199, 9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277,
    9281, 9283, 9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
    9391, 9397, 9403, 9409, 9413, 9419, 9421, 9431, 9433, 9437, 9439, 9461,
    9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, 9539, 9547, 9551,
    9587, 9601, 9613, 9619, 9623, 9629, 9631, 9643, 9649, 9661, 9677, 9679,
    9689, 9697, 9719, 9721, 9733, 9739, 9743, 9749, 9767, 9769, 9781, 9787,
    9791, 9803, 9811, 9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883,
    9887, 9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973, 0};


static long next_prime(long);

// taken from NTL, stripped of
// overflow checking
static void XGCD(long& d, long& s, long& t, long a, long b)
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
static long InvMod(long a, long n)
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
static long PowerMod(long a, long ee, long n)
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
static long GCD(long a, long b)
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

static long LCM(long a, long b) {
    long g = GCD(a, b);
    return a * (b/g);
}

static long euler_phi(long n) {
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

static void factors(long n, std::vector<long> * primes, std::vector<int> * exponents) {
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

static bool is_squarefree(long n) {
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

static bool is_fundamental_discriminant(long n) {
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

static long primitive_root(long n) {
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
    if(n == 40487) return 10;
#ifdef FLINT_VERSION
    if(n_is_prime(n)) {
        return n_primitive_root_prime(n);
    }
#endif
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

static bool is_prime_power(long q) {
    std::vector<long> primes;
    factors(q, &primes, NULL);
    if(primes.size() == 1)
        return true;
    else
        return false;
}

static bool MR_test(long n, long a) {
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

static bool is_prime(long q) {
#ifdef FLINT_VERSION
    return n_is_prime(q);
#else
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
#endif
}

static long next_prime(long n) {
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

static long odd_part(long n) {
    if(n == 0) {
        return 1;
    }
    while(n % 2 == 0) {
        n = n/2;
    }
    return n;
}

static long kronecker(long n, long m) {
    if(GCD(n,m) != 1) {
        return 0;
    }
    if(n < 0) {
        if(m % 2 == 1) {
            n = n % m;
            if(n != 0) n = n + m;
        }
        else {
            n = n % (4*m);
            if(n != 0) n += 4*m;
        }
    }
    //if(n > m) {
    //    n = n % m;
    //}
    long t = 1;
    /*
    std::cout << std::endl;
    std::cout << n << " " << m << std::endl;
    while(m > 1) {
        long m_odd = odd_part(m);
        long n_odd = odd_part(n);
        if(m_odd % 4 == 3 && n_odd % 4 == 3) t = -t;
        long x = n;
        n = m % n;
        m = x;
    }
    std::cout << t << " " << n << " " << m << std::endl;
    */
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



static long kronecker2(long n, long m) {
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

static void prime_range(std::vector<long> * primes, long end, long start = 2) {
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

static std::vector<long> prime_range(long end, long start = 2) {
    // fill primes with a list of prime numbers between
    // start and end (including start but not end)

    std::vector<long> primes;
    if(start > 2) {
        std::cerr << "that's not implemented." << std::endl;
        return primes;
    }


    // a really simple sieve...
    bool * sieve_range = new bool[end]();
    sieve_range[0] = 1;
    sieve_range[1] = 1;
    long p = 2;
    while(p < end) {
        primes.push_back(p);
        for(long k = 2*p; k < end; k += p) {
            sieve_range[k] = 1;
        }
        do {
            p++;
        } while (p < end && sieve_range[p] == 1);
    }

    delete [] sieve_range;
    return primes;
}



static void squarefree_range(std::vector<long> * squarefrees, long end, long start = 1) {
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

static void discriminant_range(std::vector<long> * discs, long end, long start = 1) {
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

static int order_mod(int n, int q) {
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

static long CRT(long a, long b, long m, long n) {
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

static long CRT(std::vector<long> a, std::vector<long> n) {
    if(a.size() == 0) return 0;
    long x = a[0];
    long q = n[0];
    for(int k = 1; k < a.size(); k++) {
        x = CRT(x, a[k], q, n[k]);
        q = q*n[k];
    }
    if(q == 1) return 0;
    return x;
}

static inline long ipow(long a, long n) {
    long z = 1;
    for(int j = 0; j < n; j++) {
        z = z * a;
    }
    return z;
}

//static std::vector<long> divisors(long n) {
//    std::vector<long> x;
//    for(int d = 1; d <= n; d++) {
//        if(n % d == 0) x.push_back(d);
//    }
//    return x;
//}

struct int_factor_t {
    int p;
    int e;
    int f; // f == p^e
};

struct long_factor_t {
    long p;
    int e;
    long f; // f == p^e
};

static int_factor_t * factor_table;
static int factor_table_size = 0;
static int mmapped_table = 0;

static void build_factor_table(int size) {
    if(mmapped_table) {
        std::cerr << "increasing the size of an mmapped factorization table is currently unsupported." << std::endl;
        exit(1);
    }
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

static void load_factor_table() {
    std::string home = std::getenv("HOME");
    std::string filelocation = home + "/include/int-factorization-table";
    int fd = open(filelocation.c_str(), O_RDONLY);
    if(fd == -1) {
        std::cerr << "error opening factor table file" << std::endl;
        exit(1);
    }
    struct stat fs;
    if(fstat(fd, &fs) == -1) {
        std::cerr << "error in fstat opening factor table." << std::endl;
        exit(1);
    }

    size_t filesize = fs.st_size;
    factor_table_size = filesize/sizeof(int_factor_t);

    factor_table = (int_factor_t *)mmap(NULL, filesize, PROT_READ, MAP_SHARED, fd, 0);
    if(factor_table == MAP_FAILED) {
        std::cerr << "error opening mmapping factor table file" << std::endl;
        exit(1);
    }
    mmapped_table = 1;
}


static long_factor_t * long_factor_table;
static long long_factor_table_size = 0;

static void build_long_factor_table(long size) {
    if(long_factor_table_size > 0) {
        delete [] long_factor_table;
    }
    if(size == 0) {
        long_factor_table_size = 0;
        return;
    }
    long_factor_table = new long_factor_t[size];
    for(long k = 0; k < size; k++) {
        long_factor_table[k] = {1,1,1};
    }
    std::vector<long> primes;
    prime_range(&primes, size);
    for(auto p : primes) {
        int e = 1;
        long f = p;
        while(f < size) {
            for(long n = f; n < size; n += f) {
                long_factor_table[n] = {p, e, f};
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


struct long_factorization_t {
    // factorization of 64 bit (signed) integers
    long_factor_t factors[16];
    long n;
    int nfactors;
};



static void factor(int n, int_factorization_t &factorization) {
    const int max_table_size = 100000000;
    int sign = 1;
    if(n < 0) {sign = -1; n = -n;}
    if(n >= factor_table_size) {
        if(mmapped_table) {
            std::cerr << "The mmapped factorization table is not large enough. Quitting." << std::endl;
            exit(0);
        }
        if(n < max_table_size) {
            build_factor_table(std::max(100000, std::min(2*n, max_table_size)));
        }
        else {
            std::cerr << "Factorization of integers > " << max_table_size << " not currently supported." << std::endl;
            std::cerr << "You might try changing the default max size in slint.h if you have lots of ram available." << std::endl;
            exit(0);
        }
    }
    factorization.n = sign * n;
    factorization.nfactors = 0;
    while(n > 1) {
        int_factor_t f = factor_table[n];
        factorization.factors[factorization.nfactors] = f;
        n = n/f.f;
        factorization.nfactors++;
    }
}

static void factor_long(long n, long_factorization_t &factorization) {
    const long max_table_size = 100000000;
    long sign = 1;
    if(n < 0) {sign = -1; n = -n;}
    if(n >= factor_table_size) {
        if(n < max_table_size) {
            build_long_factor_table(std::max(100000l, std::min(2*n, max_table_size)));
        }
        else {
            std::cerr << "Factorization of integers > " << max_table_size << " not currently supported." << std::endl;
            std::cerr << "You might try changing the default max size in slint.h if you have lots of ram available." << std::endl;
            exit(0);
        }
    }
    factorization.n = sign * n;
    factorization.nfactors = 0;
    while(n > 1) {
        long_factor_t f = long_factor_table[n];
        factorization.factors[factorization.nfactors] = f;
        n = n/f.f;
        factorization.nfactors++;
    }
}



//std::vector<int> divisors(int n) {
//    std::vector<int> x;
//    for(int d = 1; d <= n; d++) {
//        if(n % d == 0) x.push_back(d);
//    }
//    return x;
//}

static std::vector<int> divisors(int n) {
    std::vector<int> x;
    int_factorization_t fac;
    factor(n, fac);
    int e[fac.nfactors];
    for(int k = 0; k < fac.nfactors; k++) {e[k] = 0;}   // We'll write f = prod(fac.factors[k].p^e[k])
    int f = 1;                                          // but we don't directly compute it this way.
    do {
        x.push_back(f);
        int j = 0;
        while(j < fac.nfactors && e[j] + 1 > fac.factors[j].e) {
            while(e[j] > 0) {
                e[j]--;
                f = f/fac.factors[j].p;
            }
            j++;
        }
        if(j == fac.nfactors) f = 0;
        else {
            e[j]++;
            f *= fac.factors[j].p;
        }
    } while(f != 0);

    std::sort(x.begin(), x.end());

    return x;
}


static std::vector<long> divisors(long n) {
    std::vector<long> x;
    long_factorization_t fac;
    factor_long(n, fac);
    int e[fac.nfactors];
    for(int k = 0; k < fac.nfactors; k++) {e[k] = 0;}   // We'll write f = prod(fac.factors[k].p^e[k])
    int f = 1;                                          // but we don't directly compute it this way.
    do {
        x.push_back(f);
        int j = 0;
        while(j < fac.nfactors && e[j] + 1 > fac.factors[j].e) {
            while(e[j] > 0) {
                e[j]--;
                f = f/fac.factors[j].p;
            }
            j++;
        }
        if(j == fac.nfactors) f = 0;
        else {
            e[j]++;
            f *= fac.factors[j].p;
        }
    } while(f != 0);

    std::sort(x.begin(), x.end());

    return x;
}




static int ndivisors(int n) {
    int_factorization_t fac;
    factor(n, fac);
    int x = 1;
    for(int k = 0; k < fac.nfactors; k++) {
        x = x * (fac.factors[k].e + 1);
    }
    return x;
}

static int euler_phi(int n) {
    int_factorization_t f;
    factor(n, f);
    int phi = 1;
    for(int k = 0; k < f.nfactors; k++) {
        int p = f.factors[k].p;
        int e = f.factors[k].e;
        int p_to_e = f.factors[k].f;
        if(e == 1) {
            phi *= (p-1);
        }
        else {
            phi *= p_to_e/p * (p - 1);
        }
    }

    return phi;
}

static int mobius(int n) {
    int_factorization_t fac;
    factor(n, fac);
    int mu = 1;
    for(int k = 0; k < fac.nfactors; k++) {
        if(fac.factors[k].e > 1) {
            return 0;
        }
        mu = -mu;
    }
    return mu;
}

static int squarefree_part(int n) {
    int_factorization_t f;
    factor(n, f);
    int s = 1;
    for(int k = 0; k < f.nfactors; k++) {
        if(f.factors[k].e % 2 == 1) s *= f.factors[k].p;
    }
    return s;
}

#endif
