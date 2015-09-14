#include <iostream>

using namespace std;

#include <pari/pari.h>


int classnumber(int D) {
    int ret;
    pari_sp av = avma; ret=itos(classno(stoi(D))); avma = av;
    return ret;
}

int main(int argc, char ** argv) {
    pari_init(80000000, 6000000);
    int end = atoi(argv[1]);
    cout << "const int classnumbers [] = {" << endl;
    cout << 0;
    for(int k = 1; k < end; k++) {
        if(k % 10000 == 0) cerr << k << endl;
        cout << ',' << endl;
        if(k % 4 == 3 || k % 4 == 0) {
            cout << classnumber(-k);
        }
        else {
            cout << -1;
        }
    }
    cout << "};" << endl;
    return 0;
}
