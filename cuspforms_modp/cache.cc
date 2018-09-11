#include <iostream>
#include <map>
#include <tuple>
#include <mutex>

#include "cuspforms_modp.h"

using std::map;
using std::pair;
using std::tuple;
using std::cout;
using std::cerr;
using std::endl;
using std::recursive_mutex;
using std::lock_guard;

typedef tuple<int,int,int,long> space_desc_t;
static map<space_desc_t, cuspforms_modp*> cache;

static recursive_mutex cache_mutex;

cuspforms_modp * get_cuspforms_modp(DirichletCharacter &chi, int weight, long p, int verbose) {
    lock_guard<recursive_mutex> lock(cache_mutex);
    long order = order_mod(chi.m, chi.parent->q);
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

    if(verbose > 2) cout << "getting space cuspforms_modp(" << chi.parent->q << ", " << weight << ", " << chi.m << ")" << " with p == " << p << endl;
    auto result = cache.find( space_desc_t(chi.parent->q, chi.m, weight, p) );
    if(result != cache.end()) return result->second;

    cuspforms_modp * S = new cuspforms_modp(chi, weight, p, verbose);
    cache[space_desc_t(chi.parent->q, chi.m, weight, p)] = S;
    return S;
}

void clear_cuspforms_modp() {
    lock_guard<recursive_mutex> lock(cache_mutex);
    for(auto x : cache) {
        delete x.second;
    }
    cache.clear();
}
