#include "cuspforms_modp.h"

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
