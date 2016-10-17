#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <cstdint>

#include "acb.h"

using namespace std;


// Taken from http://stackoverflow.com/questions/236129/split-a-string-in-c/7408245#7408245
std::vector<std::string> split(const std::string &text, char sep) {
  std::vector<std::string> tokens;
  std::size_t start = 0, end = 0;
  while ((end = text.find(sep, start)) != std::string::npos) {
    tokens.push_back(text.substr(start, end - start));
    start = end + 1;
  }
  tokens.push_back(text.substr(start));
  return tokens;
}


void acb_set_zzzzzzzz(acb_t out, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d,
                                 fmpz_t e, fmpz_t f, fmpz_t g, fmpz_t h) {
    arf_set_fmpz_2exp(arb_midref(acb_realref(out)), a, b);
    mag_set_fmpz_2exp_fmpz(arb_radref(acb_realref(out)), c, d);
    arf_set_fmpz_2exp(arb_midref(acb_imagref(out)), e, f);
    mag_set_fmpz_2exp_fmpz(arb_radref(acb_imagref(out)), g, h);
}

int parse_textfile(acb_ptr out, string filename, int maxlines) {
    ifstream infile(filename);
    if(!infile) return 0;

    fmpz_t a, b, c, d, e, f, g, h;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(d);
    fmpz_init(e);
    fmpz_init(f);
    fmpz_init(g);
    fmpz_init(h);

    int count = 0;
    while(infile && count < maxlines) {
        string x;
        infile >> x;
        if(!infile) break;
        fmpz_set_str(a, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(b, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(c, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(d, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(e, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(f, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(g, x.c_str(), 10);
        infile >> x;
        fmpz_set_str(h, x.c_str(), 10);

        acb_set_zzzzzzzz(out + count, a, b, c, d, e, f, g, h);
        count++;
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(d);
    fmpz_clear(e);
    fmpz_clear(f);
    fmpz_clear(g);
    fmpz_clear(h);
    return count;
}


const int lines_expected = 4000;

// 3229261
struct mfheader {
    uint32_t level;
    uint32_t weight;
    uint32_t chi;
    uint32_t orbit;
    uint32_t j;
    uint32_t prec;
    uint32_t ncoeffs;
};

int write_mfheader(FILE * outfile, struct mfheader * header) {
    if(!fwrite((char*)&header->level, sizeof(header->level), 1, outfile)) return 0;
    if(!fwrite((char*)&header->weight, sizeof(header->weight), 1, outfile)) return 0;
    if(!fwrite((char*)&header->chi, sizeof(header->chi), 1, outfile)) return 0;
    if(!fwrite((char*)&header->orbit, sizeof(header->orbit), 1, outfile)) return 0;
    if(!fwrite((char*)&header->j, sizeof(header->j), 1, outfile)) return 0;
    if(!fwrite((char*)&header->prec, sizeof(header->prec), 1, outfile)) return 0;
    if(!fwrite((char*)&header->ncoeffs, sizeof(header->ncoeffs), 1, outfile)) return 0;
    return 1;
}

int main(int argc, char ** argv) {
    const char * usage =
        "convert-to-binary infile outfile [prec] [ncoeffs]\n"
        "\n"
        "Convert the old-style text format for modular form coefficient data\n"
        "to a more efficient binary format.\n"
        "\n"
        "Reads from infile, writes to outfile.\n"
        "In the output, the coeffcients will all have an implicit error interval\n"
        "with radius 2^{-prec}. If prec is not specified, of is set to 4294967295, use the largest value possible.\n"
        "If the input isn't precise enough to acheive the desired output, print an error\n"
        "message and quit. If ncoeffs is specified, the output will be truncated to the\n"
        "specified number of coefficients.\n";

    if(argc < 3) {
        cout << usage;
        return 0;
    }
    slong prec = 4294967295;
    if(argc > 3) {
        prec = atol(argv[3]);
    }

    string infilename = argv[1];
    string outfilename = argv[2];
    ifstream filelist(argv[1]);

    acb_ptr coeffs = _acb_vec_init(lines_expected);
    int lines = parse_textfile(coeffs, infilename, lines_expected);
    if(lines == 0) {
        cout << "Error reading input file." << endl;
        return -1;
    }

    if(argc > 4) {
        int z = atoi(argv[4]);
        if(z < lines) lines = z;
    }


    arf_t radius;
    arf_init(radius);
    slong max_error = -ARF_PREC_EXACT;
    for(int k = 0; k < lines; k++) {
        arf_set_mag(radius, arb_radref(acb_realref(coeffs + k)));
        slong error_bound = arf_abs_bound_lt_2exp_si(radius);
        if(error_bound > max_error)
            max_error = error_bound;

        //cout << error_bound << " ";
        arf_set_mag(radius, arb_radref(acb_imagref(coeffs + k)));
        error_bound = arf_abs_bound_lt_2exp_si(radius);
        if(error_bound > max_error)
            max_error = error_bound;
        //cout << error_bound << endl;
    }


    if(prec == 4294967295) {
        prec = -max_error - 1;
    }

    cout << prec << endl;
    cout << max_error << endl;

    if(prec > -max_error - 1) {
        cout << "Input file not precise enough." << endl;
        return -1;
    }

    if(prec > 100000) {
        cout << "Something went wrong or the precision is too large. You don't want to do this." << endl;
        return -1;
    }

    if(prec < 0) {
        cout << "negative precision not supported." << endl;
        return -1;
    }

    int a = infilename.rfind("/");
    string basefilename = infilename.substr(a+1);

    cout  << basefilename << endl;
    vector<string> z = split(basefilename, '.');
    mfheader header;
    header.level = atoi(z[0].c_str());
    header.weight = atoi(z[1].c_str());
    header.chi = atoi(z[2].c_str());
    header.orbit = 0;
    header.j = atoi(z[3].c_str());
    header.prec = prec;
    header.ncoeffs = lines;

    fmpz_t x;
    fmpz_t y;
    fmpz_init(x);
    fmpz_init(y);

    FILE * outfile = fopen(outfilename.c_str(), "w");

    cout << header.level << " " << header.weight << " " << header.chi << " " << header.j << endl;
    unsigned int magic = 3229261;
    fwrite( (void*)&magic, sizeof(magic), 1, outfile);
    write_mfheader(outfile, &header);

    for(int k = 0; k < lines; k++) {
        acb_mul_2exp_si(coeffs + k, coeffs + k, prec);
        if( !arb_get_unique_fmpz(x, acb_realref(coeffs + k)) ) {
            arb_floor(acb_realref(coeffs + k), acb_realref(coeffs + k), prec + 100);
            if( !arb_get_unique_fmpz(x, acb_realref(coeffs + k)) ) {
                cout << "error" << endl;
                exit(0);
            }
        }
        if( !arb_get_unique_fmpz(y, acb_imagref(coeffs + k)) ) {
            arb_floor(acb_imagref(coeffs + k), acb_imagref(coeffs + k), prec + 100);
            if( !arb_get_unique_fmpz(y, acb_imagref(coeffs + k)) ) {
                cout << "error" << endl;
                exit(0);
            }
        }
        fmpz_out_raw(outfile, x);
        if(header.chi != 1)
            fmpz_out_raw(outfile, y);
    }
    fclose(outfile);
    _acb_vec_clear(coeffs, lines_expected);
    return 0;
}
