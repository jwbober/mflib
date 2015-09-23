#include <iostream>
#include <cstdlib>

#include "flint/nmod_mat.h"
#include "flint/nmod_poly.h"

#include "characters.h"
#include "S2dimensions.h"
#include "trace-formula.h"
#include "classnumbers.h"

using std::vector;
using std::cerr;
using std::cout;
using std::endl;

void print_nmod_mat_t(nmod_mat_t M) {
    int nrows = nmod_mat_nrows(M);
    int ncols = nmod_mat_ncols(M);
    cout << "matrix(Integers(" << M->mod.n << "),[" << endl;
    for(int i = 0; i < nrows; i++) {
        cout << " [";
        for(int j = 0; j < ncols; j++) {
            cout << nmod_mat_entry(M, i, j);
            if(j + 1 < ncols) cout << ", ";
        }
        cout << "]";
        if(i + 1 < nrows) cout << ",";
        cout << endl;
    }
    cout << "])" << endl;
}

int main(int argc, char ** argv) {
    if(argc < 6) {
        const char * usage = "./test-tf level chi p0 ncoeffs extra_rows\n";
        cout << usage;
        return 0;
    }
    init_classnumbers();
    int level = atoi(argv[1]);
    int chi_number = atoi(argv[2]);
    long p0 = atol(argv[3]);
    int ncoeffs = atoi(argv[4]);
    int extra_rows = atoi(argv[5]);
    DirichletGroup G(level);
    DirichletCharacter chi = G.character(chi_number);

    nmod_mat_t * bases = new nmod_mat_t[level + 1];
    newspace_bases_weight2_modp(bases, ncoeffs, level, p0, chi, extra_rows);
    int q = chi.conductor();
    for(int M : divisors(level)) {
        if(M % q == 0 && M > 1) {
            int dimension;
            if(nmod_mat_nrows(bases[M]) == 0) {
                dimension = 0;
            }
            else dimension = nmod_mat_entry(bases[M], 0, 1);
            int rank = nmod_mat_rank(bases[M]);
            cout << M << " " << dimension << " " << rank << endl;
            print_nmod_mat_t(bases[M]);
            cout << endl;
        }
    }
    cout << endl;
    nmod_mat_t fullspace_basis;
    cuspform_basis_weight2_modp(fullspace_basis, ncoeffs, level, p0, chi);
    
    int dimension = nmod_mat_nrows(fullspace_basis);
    print_nmod_mat_t(fullspace_basis);
    for(int k = 0; k < dimension; k++) {
        int rank = nmod_mat_rank(fullspace_basis);
        cout << dimension - k << " " << rank << endl;
        for(int j = 0; j < nmod_mat_ncols(fullspace_basis); j++) {
            nmod_mat_entry(fullspace_basis, dimension - k - 1, j) = 0;
        }
    }
    return 0;
}
