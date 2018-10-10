#include <iostream>
#include <fstream>

using namespace std;

#include <pari/pari.h>


int classnumber(int D) {
    int ret;
    pari_sp av = avma; ret=itos(classno(stoi(D))); avma = av;
    return ret;
}

int main(int argc, char ** argv) {
    pari_init(80000000, 6000000);
    int start = atoi(argv[1]);
    int end = atoi(argv[2]);
    ofstream outfile(argv[3]);

    int zero = 0;
    for(int k = start; k < end; k++) {
        //if(k % 10000 == 0) cout << k << endl;
        if(k > 0 && (k % 4 == 3 || k % 4 == 0)) {
            int h = classnumber(-k);
            outfile.write((char*)&h, sizeof(h));
        }
        else {
            outfile.write((char*)&zero, sizeof(zero));
        }
    }
    outfile.close();
    return 0;
}
