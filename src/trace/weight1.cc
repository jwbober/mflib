#include <iostream>
#include <cstdlib>

#include "flint/nmod_mat.h"
#include "flint/nmod_poly.h"

#include "characters.h"
#include "S2dimensions.h"
#include "trace-formula.h"

using std::vector;
using std::cerr;
using std::cout;
using std::endl;


//#undef nmod_mat_entry
//static mp_limb_t& nmod_mat_entry(nmod_mat_t mat, int i, int j) {
//    if(i >= mat->r) { cout << "nmod_mat_entry() called with i out of range."; exit(0); }
//    if(j >= mat->c) { cout << "nmod_mat_entry() called with j out of range."; exit(0); }
//    return mat->rows[i][j];
//}

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

void compute_more_eisenstein_coefficients(int number_of_coefficients);

vector<long> E3Z;
vector<long> E4Z;
vector<long> E7Z;

void eisenstein3(nmod_poly_t E3, unsigned degree) {
    long mod = nmod_poly_modulus(E3);
    if(degree >= E3Z.size()) {
        compute_more_eisenstein_coefficients(degree + 1);
    }
    for(unsigned int k = 0; k < degree; k++) {
        nmod_poly_set_coeff_ui(E3, k, E3Z[k] % mod);
    }
}
void eisenstein4(nmod_poly_t E4, unsigned int degree) {
    long mod = nmod_poly_modulus(E4);
    if(degree >= E4Z.size()) {
        compute_more_eisenstein_coefficients(degree + 1);
    }
    for(unsigned int k = 0; k < degree; k++) {
        nmod_poly_set_coeff_ui(E4, k, E4Z[k] % mod);
    }
}
void eisenstein7(nmod_poly_t E7, unsigned int degree) {
    long mod = nmod_poly_modulus(E7);
    if(degree >= E7Z.size()) {
        compute_more_eisenstein_coefficients(degree + 1);
    }
    for(unsigned int k = 0; k < degree; k++) {
        nmod_poly_set_coeff_ui(E7, k, E7Z[k] % mod);
    }
}

void compute_more_eisenstein_coefficients(int number_of_coefficients) {
    for(int n = E3Z.size(); n < number_of_coefficients; n++) {
        if(n == 0) {
            E3Z.push_back(1);
            E4Z.push_back(1);
            E7Z.push_back(1);
        }
        else {
            long a3 = 0;
            long a4 = 0;
            long a7 = 0;
            for(long d : divisors(n)) {
                switch(d % 3) {
                    case 0:
                        break;
                    case 1:
                        a3 += 1;
                        break;
                    case 2:
                        a3 -= 1;
                        break;
                }
                switch(d % 4) {
                    case 0:
                    case 2:
                        break;
                    case 1:
                        a4 += 1;
                        break;
                    case 3:
                        a4 -= 1;
                        break;
                }
                switch(d % 7) {
                    case 1:
                    case 2:
                    case 4:
                        a7 += 1;
                        break;
                    case 3:
                    case 5:
                    case 6:
                        a7 -= 1;
                        break;
                    case 0:
                        break;
                }

            }
            E3Z.push_back(6 * a3);
            E4Z.push_back(4 * a4);
            E7Z.push_back(2 * a7);
        }
    }
}

int weight1_dimension_bound(int level, DirichletCharacter& chi, long& p, int extra_coeffs, int verbose) {
    //
    // Compute (an upper bound for?) the dimension of the space of weight 1
    // cusp forms mod p with character chi, where p is an appropriate prime.
    //
    // If p is 0, we choose the first appropriate prime > 10000. Otherwise
    // we reset p to the first appropriate prime >= p.
    //
    // "Appropriate" means that p - 1 is divisible by both the orders of
    // chi * chi3 and chi * chi4, where chi3 is the quadratic character
    // mod 3 and chi4 is the quadratic character mod 4.
    //
    if(chi.is_even()) return 0;

    long q3 = level;
    if(q3 % 3 != 0) q3 *= 3;  // q3 = LCM(3, level)

    long q4 = level;
    if(q4 % 2 != 0) q4 *= 2;
    if(q4 % 4 != 0) q4 *= 2;

    long primitive_index;
    int conductor = chi.conductor(&primitive_index);
    DirichletGroup G3(q3);
    DirichletGroup G4(q4);

    DirichletCharacter chi3 = G3.character(G3.index_from_primitive_character(3, 2)                                  //
                                           * G3.index_from_primitive_character(conductor, primitive_index) % q3);   // What kind of person can
    DirichletCharacter chi4 = G4.character(G4.index_from_primitive_character(4, 3)                                  // write such ugly code?
                                           * G4.index_from_primitive_character(conductor, primitive_index) % q4);   //

    int order3 = order_mod(chi3.m, q3);
    int order4 = order_mod(chi4.m, q4);

    if(order3 != order4 && order3 != 2 * order4 && 2*order3 != order4) {
        cerr << "Something went wrong. Error1." << endl;
        exit(0);
    }

    int order = std::max(order3, order4);

    if(p == 0) {
        p = order + 1;
        while(p < 10000) p += order;
    }
    else {
        if( (p - 1) % order != 0) {
            p += order - (p - 1) % order;
        }
    }
    int phiN = euler_phi(level);
    while(!n_is_prime(p) || phiN % p == 0) p += order; // can p divide phi(N)?...
                                                      // not going to think about
                                                      // that right now.

    nmod_mat_t basis3;
    nmod_mat_t basis4;

    int d3 = S2_max_dimension[q3];
    int d4 = S2_max_dimension[q4];
    int Q = d3 + d4 + 10;
    Q += extra_coeffs;

    //cout << q3 << " " << chi3.m << endl;
    //cout << q4 << " " << chi4.m << endl;
    cuspform_basis_weight2_modp(basis3, Q, q3, p, chi3, verbose);
    cuspform_basis_weight2_modp(basis4, Q, q4, p, chi4, verbose);

    //print_nmod_mat_t(basis3);
    //print_nmod_mat_t(basis4);


    d3 = nmod_mat_nrows(basis3);
    d4 = nmod_mat_nrows(basis4);
    
    //cout << nmod_mat_ncols(basis3) << endl;
    //cout << nmod_mat_ncols(basis4) << endl;

    //Q = d3 + d4 + 10;

    if(d3 != nmod_mat_rank(basis3) || d4 != nmod_mat_rank(basis4)) {
        cerr << "didn't get a matrix of full rank." << endl;
        cerr << "d3 = " << d3 << endl;
        cerr << "for chi3 = (" << q3 << ", " << chi3.m << ") rank(basis3) = " << nmod_mat_rank(basis3) << endl;
        cerr << "d4 = " << d4 << endl;
        cerr << "for chi4 = (" << q3 << ", " << chi4.m << ") rank(basis4) = " << nmod_mat_rank(basis4) << endl;
        nmod_mat_clear(basis3);
        nmod_mat_clear(basis4);
        return -1;
    }

    nmod_poly_t E3;  nmod_poly_init(E3, p);
    nmod_poly_t E4;  nmod_poly_init(E4, p);
    nmod_poly_t f;   nmod_poly_init(f, p);
    nmod_poly_t g;   nmod_poly_init(g, p);

    eisenstein3(E3, Q);
    eisenstein4(E4, Q);

    nmod_mat_t M;
    nmod_mat_init(M, d3 + d4, Q, p);
 
    for(int j = 0; j < d3; j++) {
        nmod_poly_zero(g);
        for(int k = 0; k < Q; k++) {
            nmod_poly_set_coeff_ui(g, k, nmod_mat_entry(basis3, j, k));
        }
        nmod_poly_mul(f, g, E4);
        for(int k = 0; k < Q; k++) {
            nmod_mat_entry(M, j, k) = nmod_poly_get_coeff_ui(f, k);
        }
    }
    
    for(int j = d3; j < d3 + d4; j++) {
        nmod_poly_zero(g);
        for(int k = 0; k < Q; k++) {
            nmod_poly_set_coeff_ui(g, k, nmod_mat_entry(basis4, j - d3, k));
        }
        nmod_poly_mul(f, g, E3);
        for(int k = 0; k < Q; k++) {
            nmod_mat_entry(M, j, k) = nmod_poly_get_coeff_ui(f, k);
        }
    }

    int weight2_dimension = d3 + d4 - nmod_mat_rank(M);

    nmod_mat_clear(M);
    nmod_mat_clear(basis3);
    nmod_mat_clear(basis4);
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(E3);
    nmod_poly_clear(E4);

    return weight2_dimension;
}

int main(int argc, char ** argv) {
    if(argc < 2) {
        const char * usage = 
            "./weight1 level [chi] [p0] [extracoeffs] [verbose]\n"
            "\n"
            "Compute (an upper bound for) the dimension of S_1(level, chi) mod p,\n"
            "where p is the first prime >= p0 which is 1 mod the order of chi*chi3\n"
            "and chi*chi4.\n"
            "\n"
            "The output has the form\n"
            "\n"
            "level chi p dimension\n"
            "\n"
            "where p is the prime actually used. If p0 is 0, or not set we choose\n"
            "p0 = 10000. If extracoeffs is nonzero, we that many more coefficients\n"
            "than we normally would. If chi is zero, or unset, we comput this\n"
            "dimension for every chi mod level.\n";
        cout << usage;
        return 0;
    }
    int level = atoi(argv[1]);
    int chi_number = 0;
    long p0 = 0;
    int extra_coeffs = 0;
    int verbose = 0;
    if(argc > 2) chi_number = atoi(argv[2]);
    if(argc > 3) p0 = atol(argv[3]);
    if(argc > 4) extra_coeffs = atoi(argv[4]);
    if(argc > 5) verbose = atoi(argv[5]);

    DirichletGroup G(level);
    if(chi_number != 0) {
        long p = p0;
        DirichletCharacter chi = G.character(chi_number);
        int bound = weight1_dimension_bound(level, chi, p, extra_coeffs, verbose);
        if(bound < 0)
            cout << level << " " << chi_number << " " << p << " " << '?' << endl;
        else
            cout << level << " " << chi_number << " " << p << " " << bound << endl;
    }
    else {
        for(int k = 1; k < level; k++) {
            if(GCD(k, level) != 1) continue;
            DirichletCharacter chi = G.character(k);
            long p = p0;
            int bound = weight1_dimension_bound(level, chi, p, extra_coeffs, verbose);
            if(bound < 0)
                cout << level << " " << chi_number << " " << p << " " << '?' << endl;
            else
                cout << level << " " << chi_number << " " << p << " " << bound << endl;
        }
    }
    return 0;
}
