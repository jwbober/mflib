/*
 * chimax.cc -- Program to compute the maximum of character sums
 *              and the location of the maximum.
 *
 *   The output will be a table of numbers, ascii (text) formatted by default.
 *   Each line is:
 *
 *      index parity max_index abs(max) arg(max)
 *
 *  The parity is 1 (True) is chi is odd, and 0 if it is even.
 *
 *  If binary output is requested on the command line, then this same output
 *  is written as (unsigned int) (unsigned char) (unsigned int) (double) (double).
 *  See the typedef struct outdata in the code.
 */


#include "slint.h"
#include "characters.h"


#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>

using namespace std;

pthread_mutex_t * output_mutex;     // We run multiple threads at once on a
pthread_mutex_t * next_mutex;       // single q, so we have to coordinate a
                                    // little bit.
                                    
DirichletGroup * character_group;   // We are only ever going to operate on a
long q;                             // single q in any run of this program,
                                    // hence the global variables.

bool binary = false;                // Whether to print the output as binary.
                                    // (If false, we print ascii output.)

long get_next_character_index() {
    //
    // Worker threads call this function to find the next
    // character that they should work on.
    // 
    // We start with character number 2 and proceed up to
    // q-1, except that to avoid computing the complex
    // conjugate of a character that we've already computed,
    // we compute the inverse mod q of an index we might
    // return. If it is smaller, then we already computed
    // that character.

    static long n = 2;

    long next;
    pthread_mutex_lock(next_mutex);
    
    if(n < q) {
        next = n;
        n = n + 1;
        while(GCD(next,q) != 1 || InvMod(next,q) < next) {
            next = n;
            n = n + 1;
        }
    }
    else {
        next = -1;
    }
    
    pthread_mutex_unlock(next_mutex);
    return next;
}

typedef struct {
    unsigned int m;
    unsigned char parity;
    unsigned int index;
    double absmax;
    double argmax;
} outdata;

void output_max(DirichletCharacter * chi, long index, complex<double> max) {
    if(index == -1)
        return;
    pthread_mutex_lock(output_mutex);

    if(binary) {
        outdata o;
        o.m = chi->m;
        o.parity = !(chi->is_even());
        o.index = index;
        o.absmax = abs(max);
        o.argmax = arg(max);

        cout.write((char*)&o, sizeof(o));
    }
    else {
        cout << chi->m << " ";
        cout << !(chi->is_even()) << " ";
        cout << index << " ";
        cout << abs(max) << " ";
        cout << arg(max) << endl;
    }

    pthread_mutex_unlock(output_mutex);
}

void * do_stuff(void * data) {
    long n = get_next_character_index();
    if(n < 0){
        pthread_exit(NULL);
    }

    DirichletCharacter chi;

    chi = character_group->character(n);

    long index;
    complex<double> max;
    while(n >= 0) {
        max = chi.max(&index);
        output_max(&chi, index, max);
        n = get_next_character_index();
        if(n >= 0) {
            chi = character_group->character(n);
        }
    }

}

void usage() {
    cout << "usage: max q number_of_threads" << endl;
}


int main(int argc, char ** argv) {
    cout << setprecision(17);
    if(argc < 3) {
        usage();
        return 1;
    }
    if(argc >= 4) {
        binary = atoi(argv[3]);
    }

    q = atol(argv[1]);

    int number_of_threads = atoi(argv[2]);

    character_group = new DirichletGroup(q);

    output_mutex = new pthread_mutex_t;
    next_mutex = new pthread_mutex_t;

    pthread_mutex_init(output_mutex, NULL);
    pthread_mutex_init(next_mutex, NULL);

    pthread_t threads[number_of_threads];

    for(int n = 0; n < number_of_threads; n++) {
        pthread_create(&threads[n], NULL, do_stuff, NULL);
    }
    for(int n = 0; n < number_of_threads; n++) {
        pthread_join(threads[n], NULL);
    }
    delete character_group;
        
    return 0;
}
