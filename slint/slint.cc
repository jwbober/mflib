//
// Simple Library for Number Theory
//

#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <mutex>

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "slint.h"

static long next_prime(long);

static int_factor_t * factor_table;
static int factor_table_size = 0;
static int mmapped_table = 0;

void build_factor_table(int size) {
    static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    if(mmapped_table) {
        std::cerr << "increasing the size of an mmapped factorization table is currently unsupported." << std::endl;
        exit(1);
    }
    if(size <= factor_table_size) return;
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

void load_factor_table() {
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

void build_long_factor_table(long size) {
    static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    if(size <= long_factor_table_size) return;
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

void factor(int n, int_factorization_t &factorization) {
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

void factor_long(long n, long_factorization_t &factorization) {
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

std::vector<int> divisors(int n) {
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

std::vector<long> divisors(long n) {
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

int ndivisors(int n) {
    int_factorization_t fac;
    factor(n, fac);
    int x = 1;
    for(int k = 0; k < fac.nfactors; k++) {
        x = x * (fac.factors[k].e + 1);
    }
    return x;
}

int euler_phi(int n) {
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

int mobius(int n) {
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

int squarefree_part(int n) {
    int_factorization_t f;
    factor(n, f);
    int s = 1;
    for(int k = 0; k < f.nfactors; k++) {
        if(f.factors[k].e % 2 == 1) s *= f.factors[k].p;
    }
    return s;
}
