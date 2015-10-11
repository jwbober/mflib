import sys
import subprocess
from dirichlet_conrey import *
prime = sys.argv[3]
def do_q(q):
    G = DirichletGroup_conrey(q)
    for orbit in G.galois_orbits():
        chi = orbit[0]
        if not chi.is_odd():
            continue
        correct_orders = True
        #for psi in chi.decomposition():
        #    if is_squarefree(psi.modulus()) and psi.multiplicative_order() not in [2,3,4,5]:
        #        correct_orders = False
        #        break
        #if not correct_orders:
        #    continue
        command = 'weight1 {} {} {} 100'.format(q, chi.number(), prime)
        print command
        sys.stdout.flush()
        #output = subprocess.check_output(command, shell=True)
        #level, character, prime, dimension = output.split()
        #dimension = int(dimension)
        #if dimension != 0:
        #    print level, character, prime, dimension, dimension * len(orbit)
        #    sys.stdout.flush()

if __name__ == "__main__":
    import sys
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    for q in range(start, end):
        do_q(q)
