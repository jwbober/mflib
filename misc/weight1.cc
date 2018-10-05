#include <iostream>

#include "flint/nmod_poly.h"
#include "flint/nmod_mat.h"

#include "cuspforms_modp.h"
#include "classnumbers.h"

using namespace std;

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

void compute_more_eisenstein_coefficients(int number_of_coefficients) {
    for(int n = E3Z.size(); n < number_of_coefficients; n++) {
        if(n == 0) {
            E3Z.push_back(1);
            E4Z.push_back(1);
            //E7Z.push_back(1);
        }
        else {
            long a3 = 0;
            long a4 = 0;
            //long a7 = 0;
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
                /*
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
                */
            }
            E3Z.push_back(6 * a3);
            E4Z.push_back(4 * a4);
            //E7Z.push_back(2 * a7);
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

    // What kind of person can write such horrible code?
    // Some of this Dirichlet character stuff really need to be rewritten...
    DirichletCharacter chi3 =
        G3.character(G3.index_from_primitive_character(3, 2)
        * G3.index_from_primitive_character(conductor, primitive_index) % q3);
    DirichletCharacter chi4 =
        G4.character(G4.index_from_primitive_character(4, 3)
        * G4.index_from_primitive_character(conductor, primitive_index) % q4);

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

    cuspforms_modp * S3 = get_cuspforms_modp(chi3, 2, p);
    cuspforms_modp * S4 = get_cuspforms_modp(chi4, 2, p);

    int d3 = S3->dimension();
    int d4 = S4->dimension();
    int Q = d3 + d4 + 10;

    Q += extra_coeffs;

    S3->basis(basis3, Q);
    S4->basis(basis4, Q);

    if(d3 != nmod_mat_rank(basis3) || d4 != nmod_mat_rank(basis4)) {
        cerr << "didn't get a matrix of full rank." << endl;
        cerr << "d3 = " << d3 << endl;
        cerr << "for chi3 = (" << q3 << ", " << chi3.m << ") rank(basis3) = " << nmod_mat_rank(basis3) << endl;
        cerr << "d4 = " << d4 << endl;
        cerr << "for chi4 = (" << q3 << ", " << chi4.m << ") rank(basis4) = " << nmod_mat_rank(basis4) << endl;
        cerr << "I think maybe this shouldn't happen? If it does, see extra_coefficients larger." << endl;
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
    //nmod_mat_init(M, d3 + d4, Q, p);
    nmod_mat_init(M, Q, d3 + d4, p);

    for(int j = 0; j < d3; j++) {
        nmod_poly_zero(g);
        for(int k = 0; k < Q; k++) {
            nmod_poly_set_coeff_ui(g, k, nmod_mat_entry(basis3, j, k));
        }
        nmod_poly_mul(f, g, E4);
        for(int k = 0; k < Q; k++) {
            //nmod_mat_entry(M, j, k) = nmod_poly_get_coeff_ui(f, k);
            nmod_mat_entry(M, k, j) = nmod_poly_get_coeff_ui(f, k);
        }
    }

    for(int j = d3; j < d3 + d4; j++) {
        nmod_poly_zero(g);
        for(int k = 0; k < Q; k++) {
            nmod_poly_set_coeff_ui(g, k, nmod_mat_entry(basis4, j - d3, k));
        }
        nmod_poly_mul(f, g, E3);
        for(int k = 0; k < Q; k++) {
            //nmod_mat_entry(M, j, k) = nmod_poly_get_coeff_ui(f, k);
            nmod_mat_entry(M, k, j) = nmod_poly_get_coeff_ui(f, k);
        }
    }

    int weight2_dimension = d3 + d4 - nmod_mat_rank(M);
    nmod_mat_t nullspace;
    nmod_mat_init(nullspace, d3 + d4, weight2_dimension, p);
    nmod_mat_nullspace(nullspace, M);

    nmod_mat_t half_nullspace;
    nmod_mat_t transpose_basis3;
    nmod_mat_init(transpose_basis3, Q + 1, d3, p);

    nmod_mat_transpose(transpose_basis3, basis3);

    nmod_mat_t weight1_basis;
    nmod_mat_init(weight1_basis, Q + 1, weight2_dimension, p);
    nmod_mat_window_init(half_nullspace, nullspace, 0, 0, d3, weight2_dimension);
    nmod_mat_mul(weight1_basis, transpose_basis3, half_nullspace);

    // we should now have some weight 2 forms (in the matrix weight1_basis)
    // which should give us weight 1 forms when divided by E3.

    nmod_poly_t E3_inverse;
    nmod_poly_init(E3_inverse, p);
    nmod_poly_inv_series(E3_inverse, E3, Q);
    cout << p << endl;
    for(int k = 0; k < weight2_dimension; k++) {
        for(int j = 0; j < Q; j++) {
            nmod_poly_set_coeff_ui(g, j, nmod_mat_entry(weight1_basis, j, k));
        }
        nmod_poly_mul(f, g, E3_inverse);
        unsigned long z = 0;
        int j = 1;
        while(z == 0) {
            z = nmod_poly_get_coeff_ui(f, j);
            j++;
        }
        z = InvMod(z, p);
        nmod_poly_scalar_mul_nmod(f, f, z);
        for(int j = 0; j < Q/2; j++) {
            cout << nmod_poly_get_coeff_ui(f, j) << " ";
        }
        cout << endl;
        cout << endl;
    }
    //nmod_mat_print_pretty(weight1_basis);

    nmod_mat_clear(M);
    nmod_mat_clear(basis3);
    nmod_mat_clear(basis4);
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(E3);
    nmod_poly_clear(E4);

    return weight2_dimension;
}

int weight1_basis_data(nmod_mat_t basis_transformation, cuspforms_modp ** S3, int level, DirichletCharacter& chi, long& p, int extra_coeffs, int verbose) {
    //
    // Compute a matrix basis_transformation such that a basis for the
    // weight 1 space mod p can be computed as
    // (*S3)->basis(basis3, ncoeffs);
    // nmod_mat_transpose(basis3_transpose, basis3);
    // nmod_mat_mul(weight1_basis, basis3_transpose, basis_transformation);
    //
    // and then multiplying the resulting basis by the inverse of the eisenstein series E3.
    //
    // Also returns the dimension of the weight 1 space.

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

    // What kind of person can write such horrible code?
    // Some of this Dirichlet character stuff really need to be rewritten...
    DirichletCharacter chi3 =
        G3.character(G3.index_from_primitive_character(3, 2)
        * G3.index_from_primitive_character(conductor, primitive_index) % q3);
    DirichletCharacter chi4 =
        G4.character(G4.index_from_primitive_character(4, 3)
        * G4.index_from_primitive_character(conductor, primitive_index) % q4);

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

    *S3 = get_cuspforms_modp(chi3, 2, p);
    cuspforms_modp * S4 = get_cuspforms_modp(chi4, 2, p);

    int d3 = (*S3)->dimension();
    int d4 = S4->dimension();
    int Q = d3 + d4 + 10;

    Q += extra_coeffs;

    (*S3)->basis(basis3, Q);
    S4->basis(basis4, Q);

    if(d3 != nmod_mat_rank(basis3) || d4 != nmod_mat_rank(basis4)) {
        cerr << "didn't get a matrix of full rank." << endl;
        cerr << "d3 = " << d3 << endl;
        cerr << "for chi3 = (" << q3 << ", " << chi3.m << ") rank(basis3) = " << nmod_mat_rank(basis3) << endl;
        cerr << "d4 = " << d4 << endl;
        cerr << "for chi4 = (" << q3 << ", " << chi4.m << ") rank(basis4) = " << nmod_mat_rank(basis4) << endl;
        cerr << "I think maybe this shouldn't happen? If it does, see extra_coefficients larger." << endl;
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
    //nmod_mat_init(M, d3 + d4, Q, p);
    nmod_mat_init(M, Q, d3 + d4, p);

    for(int j = 0; j < d3; j++) {
        nmod_poly_zero(g);
        for(int k = 0; k < Q; k++) {
            nmod_poly_set_coeff_ui(g, k, nmod_mat_entry(basis3, j, k));
        }
        nmod_poly_mul(f, g, E4);
        for(int k = 0; k < Q; k++) {
            //nmod_mat_entry(M, j, k) = nmod_poly_get_coeff_ui(f, k);
            nmod_mat_entry(M, k, j) = nmod_poly_get_coeff_ui(f, k);
        }
    }

    for(int j = d3; j < d3 + d4; j++) {
        nmod_poly_zero(g);
        for(int k = 0; k < Q; k++) {
            nmod_poly_set_coeff_ui(g, k, nmod_mat_entry(basis4, j - d3, k));
        }
        nmod_poly_mul(f, g, E3);
        for(int k = 0; k < Q; k++) {
            //nmod_mat_entry(M, j, k) = nmod_poly_get_coeff_ui(f, k);
            nmod_mat_entry(M, k, j) = nmod_poly_get_coeff_ui(f, k);
        }
    }

    int weight2_dimension = d3 + d4 - nmod_mat_rank(M);
    nmod_mat_t nullspace;
    nmod_mat_init(nullspace, d3 + d4, weight2_dimension, p);
    nmod_mat_nullspace(nullspace, M);

    nmod_mat_t half_nullspace;

    nmod_mat_window_init(half_nullspace, nullspace, 0, 0, d3, weight2_dimension);
    nmod_mat_init_set(basis_transformation, half_nullspace);

    nmod_mat_window_clear(half_nullspace);

    nmod_mat_clear(M);
    nmod_mat_clear(basis3);
    nmod_mat_clear(basis4);
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(E3);
    nmod_poly_clear(E4);
    return weight2_dimension;

    //nmod_mat_mul(weight1_basis, transpose_basis3, half_nullspace);

    // we should now have some weight 2 forms (in the matrix weight1_basis)
    // which should give us weight 1 forms when divided by E3.

    /*
    nmod_poly_t E3_inverse;
    nmod_poly_init(E3_inverse, p);
    nmod_poly_inv_series(E3_inverse, E3, Q);
    cout << p << endl;
    for(int k = 0; k < weight2_dimension; k++) {
        for(int j = 0; j < Q; j++) {
            nmod_poly_set_coeff_ui(g, j, nmod_mat_entry(weight1_basis, j, k));
        }
        nmod_poly_mul(f, g, E3_inverse);
        unsigned long z = 0;
        int j = 1;
        while(z == 0) {
            z = nmod_poly_get_coeff_ui(f, j);
            j++;
        }
        z = InvMod(z, p);
        nmod_poly_scalar_mul_nmod(f, f, z);
        for(int j = 0; j < Q/2; j++) {
            cout << nmod_poly_get_coeff_ui(f, j) << " ";
        }
        cout << endl;
        cout << endl;
    }
    //nmod_mat_print_pretty(weight1_basis);

    nmod_mat_clear(M);
    nmod_mat_clear(basis3);
    nmod_mat_clear(basis4);
    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(E3);
    nmod_poly_clear(E4);

    return weight2_dimension;
    */
}


void weight1_basis(nmod_mat_t basis, nmod_mat_t basis_transformation, cuspforms_modp * S3, long p, int ncoeffs, int verbose) {
    nmod_mat_t basis3, basis3_transpose;
    S3->basis(basis3, ncoeffs);

    nmod_mat_init(basis3_transpose, nmod_mat_ncols(basis3), nmod_mat_nrows(basis3), p);
    nmod_mat_init(basis, nmod_mat_nrows(basis3_transpose), nmod_mat_ncols(basis_transformation), p);

    nmod_mat_transpose(basis3_transpose, basis3);
    nmod_mat_mul(basis, basis3_transpose, basis_transformation);

    int weight1_dimension = nmod_mat_ncols(basis_transformation);

    nmod_poly_t E3, f, g;
    nmod_poly_t E3_inverse;

    nmod_poly_init(E3, p);
    nmod_poly_init(f, p);
    nmod_poly_init(g, p);
    nmod_poly_init(E3_inverse, p);

    eisenstein3(E3, ncoeffs);
    nmod_poly_inv_series(E3_inverse, E3, ncoeffs);

    for(int k = 0; k < weight1_dimension; k++) {
        for(int j = 0; j < ncoeffs; j++) {
            nmod_poly_set_coeff_ui(g, j, nmod_mat_entry(basis, j, k));
        }
        nmod_poly_mul(f, g, E3_inverse);
        unsigned long z = 0;
        int j = 1;
        while(z == 0) {
            z = nmod_poly_get_coeff_ui(f, j);
            j++;
        }
        z = InvMod(z, p);
        nmod_poly_scalar_mul_nmod(f, f, z);
        for(int j = 0; j < ncoeffs; j++) {
            cout << nmod_poly_get_coeff_ui(f, j) << " ";
        }
        cout << endl;
        cout << endl;
    }

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
    init_classnumbers();
    int level = atoi(argv[1]);
    int chi_number = 0;
    long p0 = 0;
    int ncoeffs = 0;
    int verbose = 0;
    if(argc > 2) chi_number = atoi(argv[2]);
    if(argc > 3) p0 = atol(argv[3]);
    if(argc > 4) ncoeffs = atoi(argv[4]);
    if(argc > 5) verbose = atoi(argv[5]);

    DirichletGroup G(level);
    if(chi_number != 0) {
        long p = p0;
        DirichletCharacter chi = G.character(chi_number);
        nmod_mat_t basis_transformation;
        nmod_mat_t basis;
        cuspforms_modp * S3;
        //int dimension = weight1_basis_data(basis_transformation, &S3, level, chi, p, 20, verbose);
        weight1_basis_data(basis_transformation, &S3, level, chi, p, 20, verbose);
        weight1_basis(basis, basis_transformation, S3, p, ncoeffs, verbose);
    }
    else {
        for(int k = 1; k < level; k++) {
            if(GCD(k, level) != 1) continue;
            DirichletCharacter chi = G.character(k);
            long p = p0;
            int bound = weight1_dimension_bound(level, chi, p, 20, verbose);
            if(bound < 0)
                cout << level << " " << chi.m << " " << p << " " << '?' << endl;
            else
                cout << level << " " << chi.m << " " << p << " " << bound << endl;
        }
    }
    return 0;
}
