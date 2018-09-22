#include <ostream>

#include <cstring>

#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"

using namespace std;

ostream& operator << (ostream& out, fmpz_t in) {
    char * outstr = fmpz_get_str(NULL, 10, in);
    out.write(outstr, strlen(outstr));
    flint_free(outstr);
    return out;
}

ostream& operator << (ostream& out, fmpz_poly_t in) {
    for(int k = 0; k <= fmpz_poly_degree(in); k++) {
        if(k != 0) out << " ";
        out << fmpz_poly_get_coeff_ptr(in, k);
    }
    return out;
}
