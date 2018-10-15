#include <vector>
#include "flint/nmod_poly.h"

#include "slint.h"

//
// XXX not thread safe
//

using namespace std;

static vector<long> E3Z;
static vector<long> E4Z;
static void compute_more_eisenstein_coefficients(int number_of_coefficients) {
    for(int n = E3Z.size(); n < number_of_coefficients; n++) {
        if(n == 0) {
            E3Z.push_back(1);
            E4Z.push_back(1);
        }
        else {
            long a3 = 0;
            long a4 = 0;
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
            }
            E3Z.push_back(6 * a3);
            E4Z.push_back(4 * a4);
        }
    }
}

void eisenstein3_weight1_modp(nmod_poly_t E3, unsigned int degree) {
    long mod = nmod_poly_modulus(E3);
    if(degree >= E3Z.size()) {
        compute_more_eisenstein_coefficients(degree + 1);
    }
    for(unsigned int k = 0; k < degree; k++) {
        nmod_poly_set_coeff_ui(E3, k, E3Z[k] % mod);
    }
}
void eisenstein4_weight1_modp(nmod_poly_t E4, unsigned int degree) {
    long mod = nmod_poly_modulus(E4);
    if(degree >= E4Z.size()) {
        compute_more_eisenstein_coefficients(degree + 1);
    }
    for(unsigned int k = 0; k < degree; k++) {
        nmod_poly_set_coeff_ui(E4, k, E4Z[k] % mod);
    }
}


