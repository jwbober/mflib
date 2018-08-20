#include "acb.h"

#include <ostream>
#include <cstdio>

using namespace std;

ostream& operator << (ostream& out, arf_t in) {
    char * outstr;
    size_t length;
    FILE * pseudofile = open_memstream(&outstr, &length);
    arf_fprintd(pseudofile, in, out.precision());
    fclose(pseudofile);
    out.write(outstr, length);
    free(outstr);
    return out;
}

ostream& operator << (ostream& out, mag_t in) {
    char * outstr;
    size_t length;
    FILE * pseudofile = open_memstream(&outstr, &length);
    mag_fprintd(pseudofile, in, out.precision());
    fclose(pseudofile);
    out.write(outstr, length);
    free(outstr);
    return out;
}

ostream& operator << (ostream& out, arb_t in) {
    return out << '(' << arb_midref(in) << " +/- " << arb_radref(in) << ")";
}

ostream& operator << (ostream& out, acb_t in) {
    return out << '[' << acb_realref(in) << ", " << acb_imagref(in) << "]";
}
