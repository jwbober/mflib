#include "arb_poly.h"

// compute gcd(a,b) given its degree
// returns 0 on success, -1 if some leading coeff contains 0
int arb_poly_gcd(arb_poly_t res,const arb_poly_t a,const arb_poly_t b,slong d, slong prec) {
	static arb_poly_t r0,r1,q,temp;
	static int init;
	int i;
	arb_ptr p;
        
	if (!init) {
		arb_poly_init(r0); arb_poly_init(r1);
		arb_poly_init(q);  arb_poly_init(temp);
		init = 1;
	}
        printf("arb_poly_gcd %ld %ld %ld\n", arb_poly_degree(a), arb_poly_degree(b), d);
	if (arb_poly_degree(a) >= arb_poly_degree(b)) {
		arb_poly_set(r0,a);
		arb_poly_set(r1,b);
	} else {
		arb_poly_set(r1,a);
		arb_poly_set(r0,b);
	}
        if(arb_poly_degree(r1) == -1) {
            if(arb_poly_degree(r0) == d) {
                arb_poly_set(res,r0);
                return 0;
            }
            else return -1;
        }
	for (;;) {
                printf("in arb_poly_gcd degree is %ld.\n", arb_poly_degree(r1));
		arb_poly_set(temp,r0);
		arb_poly_set(r0,r1);
		arb_poly_divrem(q,r1,temp,r1,prec);
		if ((i=arb_poly_degree(r1)) < d) break;
		p = arb_poly_get_coeff_ptr(r1,i);
		if (arb_contains_zero(p)) return -1;
		arb_poly_scalar_div(r1,r1,p,prec);
		arb_one(arb_poly_get_coeff_ptr(r1,i));
	}

	arb_poly_set(res,r0);
	return 0;
}
