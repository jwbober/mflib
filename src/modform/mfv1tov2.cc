#include <iostream>
#include <string>
#include "acb.h"
#include "mfformat.h"

#include "slint.h"

using namespace std;

int main(int argc, char ** argv) {
    FILE * infile = fopen(argv[1], "r");

    string outfilename = argv[2];
    string::size_type pos = outfilename.rfind("/");
    if(pos != string::npos) {
        string outpath = outfilename.substr(0, pos + 1);
        string command = "mkdir -p " + outpath;
        system(command.c_str());
    }

    FILE * outfile = fopen(argv[2], "w");

    struct mfheader header;
    read_mfheader(infile, &header);
    if(header.version != MFV1) {
        cout << "not in version 1 format" << endl;
        return -1;
    }

    header.version = MFV2;
    if(!write_mfheader(outfile, &header)) {
        cout << "error writing outfile" << endl;
        return -1;
    }
    fmpz_t x;
    fmpz_t y;
    fmpz_init(x);
    fmpz_init(y);
    int k = 0;
    int p = prime_powers_table[k];
    for(int j = 0; j < header.ncoeffs; j++) {
        fmpz_inp_raw(x, infile);
        if(header.chi != 1) fmpz_inp_raw(y, infile);
        if(j + 1 == p) {
            fmpz_out_raw(outfile, x);
            if(header.chi != 1) fmpz_out_raw(outfile, y);
            k = k + 1;
            p = prime_powers_table[k];
        }
    }

    fclose(outfile);
    return 0;
}
