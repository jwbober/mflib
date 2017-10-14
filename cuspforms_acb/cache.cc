#include <iostream>
#include <map>
#include <tuple>

#include "cuspforms_acb.h"

using std::map;
using std::pair;
using std::tuple;
using std::cout;
using std::cerr;
using std::endl;

typedef tuple<int,int,int> space_desc_t;
static map<space_desc_t, cuspforms_acb*> cache;

cuspforms_acb * get_cuspforms_acb(DirichletCharacter &chi, int weight, int nthreads, int verbose) {
    if(verbose > 2) cout << "getting space cuspforms_acb(" << chi.parent->q << ", " << weight << ", " << chi.m << ")" << endl;
    auto result = cache.find( space_desc_t(chi.parent->q, chi.m, weight) );
    if(result != cache.end()) return result->second;

    cuspforms_acb * S = new cuspforms_acb(chi, weight, nthreads, verbose);
    cache[space_desc_t(chi.parent->q, chi.m, weight)] = S;
    return S;
}

void clear_cuspform_cache() {
    for(auto S : cache) {
        delete S.second;
    }
    cache.clear();
}
