import sys
from dirichlet_conrey import *

level = int(sys.argv[1])
chi_number = int(sys.argv[2])
#start = int(sys.argv[3])
#end = int(sys.argv[4])
NN = int(sys.argv[3])
MM = int(sys.argv[4])
p = int(sys.argv[5])

def print_traces(chi, start, end):
    L = Newforms(chi, names='a')
    def trace(n):
        if len(L) == 0:
            return 0
        return sum(f[n].trace() for f in L)

    q = chi.level()
    def trace_product(m, n):
        if m*n == 0:
            return 0
        T = 0
        g = gcd(m,n)
        for d in divisors(g):
            if gcd(d, q) == 1:
                T = T + trace(m*n/d^2) * chi(d) * d

        return T
    sys.stdout.write(str(chi.level()) + " ")
    for m in range(start, end + 1):
        sys.stdout.write(str(trace(m) % p))
        if(m < end):
            sys.stdout.write('\t')
    print
    sys.stdout.flush()

def print_TnTm(chi, NN, MM):
    L = Newforms(chi, names='a')
    def trace(n):
        if len(L) == 0:
            return 0
        return sum(f[n].trace() for f in L)

    q = chi.level()
    def trace_product(m, n):
        if m*n == 0:
            return 0
        T = 0
        g = gcd(m,n)
        for d in divisors(g):
            if gcd(d, q) == 1:
                T = T + trace(m*n/d^2) * chi(d) * d
        return T

    sys.stdout.write(str(chi.level()) + "\n")
    for n in range(NN + 1):
        for m in range(MM + 1):
            z = trace_product(n,m)
            sys.stdout.write(str(ZZ(trace_product(n, m)) % p))
            if m < MM:
                sys.stdout.write('\t')
        print
        sys.stdout.flush()
    print
            
    sys.stdout.flush()



chi = DirichletGroup_conrey(level)[chi_number].sage_character()
q = chi.conductor()
print p
for M in divisors(level):
    if M == 1:
        continue
    if M % q == 0:
        print_TnTm(chi.restrict(M), NN, MM)
