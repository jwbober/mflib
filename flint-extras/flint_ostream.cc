#include <ostream>

#include <cstring>

#include "flint/fmpz.h"

using namespace std;

ostream& operator << (ostream& out, fmpz_t in) {
    char * outstr = fmpz_get_str(NULL, 10, in);
    out.write(outstr, strlen(outstr));
    flint_free(outstr);
    return out;
}
