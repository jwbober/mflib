#include <iostream>
#include <string>

#include <cstdio>

#include "mfformat.h"

#include "acb.h"

using namespace std;

int count = 0;

string outpath;
int callback(struct mfheader * header, int coeff_datasize, const void * coeff_data, acb_ptr coeffs) {
    string outfilename = outpath + "/"
                            + to_string(header->level) + "."
                            + to_string(header->weight) + "."
                            + to_string(header->chi) + "."
                            + to_string(header->j);
    FILE * outfile = fopen(outfilename.c_str(), "w");
    if(!outfile) return 1;
    cout << "writing to " << outfilename << endl;
    if(!write_mfheader(outfile, header)) return 1;
    if(fwrite(coeff_data, coeff_datasize, 1, outfile) != 1) return 1;
    fclose(outfile);
    return 0;
}

int main(int argc, char ** argv) {
    outpath = argv[2];
    iterate_through_sqlitefile(argv[1], callback, 0);
    return 0;
}
