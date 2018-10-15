#include <iostream>
#include <map>
#include <tuple>
#include <mutex>

#include "cuspforms_weight1_modp.h"

using std::map;
using std::pair;
using std::tuple;
using std::cout;
using std::cerr;
using std::endl;
using std::recursive_mutex;
using std::lock_guard;

typedef tuple<int,int,long> space_desc_t;
static map<space_desc_t, cuspforms_weight1_modp*> cache;

static recursive_mutex cache_mutex;

cuspforms_weight1_modp * get_cuspforms_weight1_modp(DirichletCharacter &chi,
                                                    long p,
                                                    int verbose) {
    lock_guard<recursive_mutex> lock(cache_mutex);
    long order = order_mod(chi.m, chi.parent->q);
    if(order % 2 != 0) order *= 2; // this may not be the order anymore, of course,
                                   // but we need the main point is that p needs to
                                   // split in all appropriate fields
                                   // (We may want to do more work here later to take
                                   //  into account the dihedral forms and other conditions...)
    if(p == 0) {
        p = order + 1;
        while(p < 10000) p += order;
    }
    else {
        if( (p - 1) % order != 0) {
            p += order - (p - 1) % order;
        }
    }
    while(!is_prime(p)) {p += order;}

    if(verbose > 2) cout << "getting space cuspforms_weight1_modp("
                         << chi.parent->q << ", "
                         << chi.m << ")" << " with p == " << p << endl;

    auto result = cache.find( space_desc_t(chi.parent->q, chi.m, p) );
    if(result != cache.end()) return result->second;

    cuspforms_weight1_modp * S = new cuspforms_weight1_modp(chi, p, verbose);
    cache[space_desc_t(chi.parent->q, chi.m, p)] = S;
    return S;
}

void clear_cuspforms_weight1_modp() {
    lock_guard<recursive_mutex> lock(cache_mutex);
    for(auto x : cache) {
        delete x.second;
    }
    cache.clear();
}
