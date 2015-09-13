import sys
from dirichlet_conrey import *

dimensions = {}
max_level = 0
for line in open('weight1_dimensions'):
    q, chi, p, d = line.split()
    q = int(q)
    chi = int(chi)
    d = int(d)
    dimensions[(q, chi)] = d
    if q > max_level:
        max_level = q

dihedral_dimensions = [0] * 1001
for line in open('dihedral.txt'):
    q, something, something_else, blank = line.split(',')
    q = int(q)
    dihedral_dimensions[q] += 1

total_dimensions = [0] * 1001

for q in range(1, max_level + 1):
    G = DirichletGroup_conrey(q)
    for orbit in G.galois_orbits():
        chi = orbit[0]
        if (q, chi.number()) not in dimensions:
            continue
        d = dimensions[(q, chi.number())]
        q0 = chi.conductor()
        for M1 in divisors(q):
            if M1 == q:
                continue
            if M1 % q0 == 0:
                G2 = DirichletGroup_conrey(M1)
                psi_orbit = G2.from_sage_character(chi.sage_character().restrict(M1)).galois_orbit()
                for psi in psi_orbit:
                    if (M1, psi.number()) in dimensions:
                        d -= number_of_divisors(q/M1) * dimensions[(M1, psi.number())]
        dimensions[(q, chi.number())] = d
        if d != 0:
            print q, chi.number(), d, d * len(orbit)
            total_dimensions[q] += d * len(orbit)
            sys.stdout.flush()

print
for q in range(1, max_level + 1):
    print q, total_dimensions[q], dihedral_dimensions[q], total_dimensions[q] - dihedral_dimensions[q]
