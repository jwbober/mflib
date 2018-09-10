#include "arb_poly.h"

// sgn(f(a)*f(b)), where x=[a,b]
int arb_poly_sign_change(arb_poly_t f,arb_srcptr x,slong pr) {
	static arb_t a,b;
	static int init;

	if (!init) {
		arb_init(a); arb_init(b);
		init = 1;
	}
	arb_get_interval_arf(arb_midref(a),arb_midref(b),x,pr);
	mag_zero(arb_radref(a)); mag_zero(arb_radref(b));
	arb_poly_evaluate(a,f,a,pr);
	arb_poly_evaluate(b,f,b,pr);
	arb_mul(a,a,b,pr);
	if (arb_is_positive(a)) return 1;
	if (arb_is_negative(a)) return -1;
	return 0;
}
