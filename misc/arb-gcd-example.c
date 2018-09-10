#include <string.h>

#include "arb-extras.h"

static slong prec;

// (f(t*z)+f(t/z))/2 and (f(t*z)-f(t/z))/(2*I), where z=e(a/q)
static void twist_poly(arb_poly_t r1,arb_poly_t r2,const fmpz_poly_t f,slong a,slong q) {
	static arb_t c,s,t;
	static int init;
	slong d,n;

	if (!init) {
		arb_init(c); arb_init(s);
		arb_init(t);
		init = 1;
	}
	d = fmpz_poly_degree(f);
	arb_poly_truncate(r1,d+1);
	arb_poly_truncate(r2,d+1);
	for (;d>=0;d--) {
		n = a*d % q;
		if (n < 0) n += q;
		arb_set_si(t,2*n);
		arb_div_si(t,t,q,prec);
		arb_sin_cos_pi(s,c,t,prec);

		arb_mul_fmpz(t,c,fmpz_poly_get_coeff_ptr(f,d),prec);
		if (4*n == q || 4*n == 3*q) arb_zero(t);
		arb_poly_set_coeff_arb(r1,d,t);

		arb_mul_fmpz(t,s,fmpz_poly_get_coeff_ptr(f,d),prec);
		if (!n || 2*n == q) arb_zero(t);
		arb_poly_set_coeff_arb(r2,d,t);
	}
}

static int phase(acb_srcptr a,int order) {
	static arb_t t,pi;
	static fmpz_t z;
	static int init;
	int res;

	if (order == 1) return 0;
	if (!init) {
		arb_init(t); arb_init(pi);
		fmpz_init(z);
		init = 1;
	}
	acb_arg(t,a,prec);
	arb_const_pi(pi,prec);
	arb_div(t,t,pi,prec);
	arb_mul_si(t,t,order,prec);
	if (!arb_get_unique_fmpz(z,t)) {
		fprintf(stderr,"cannot determine phase in arbgcd\n");
		exit(0);
	}
	res = fmpz_get_si(z) % order;
	if (res < 0) res += order;
	return res;
}

static unsigned int gcd(unsigned int x,unsigned int y) {
	int k,l;

	if (!x) return y;
	if (!y) return x;

	// __builtin_ctz() counts bottom 0 bits
	// compile with -march=corei7
	k = __builtin_ctzl(x); x >>= k;
	l = __builtin_ctzl(y); y >>= l;
	while (x != y)
	if (x < y)
		y -= x, y >>= __builtin_ctzl(y);
	else
		x -= y, x >>= __builtin_ctzl(x);
	return x<<(k<l?k:l);
}

// returns 0 on success, -1 on error (in which case res is undefined)
int bober(int *res,const fmpz_poly_t *poly,int npolys,const acb_t *a,int ncoeffs,int order,slong pr) {
	static arb_poly_t g;
	static arb_t t;
	static acb_t x,z;
	static int init;
	arb_poly_t *f;
	int i,j,u,phi,retval=0,s;

	if (!init) {
		arb_poly_init(g); arb_init(t);
		acb_init(x); acb_init(z);
		init = 1;
	}

	prec = pr;
	f = (arb_poly_t *)malloc(npolys*sizeof(f[0]));
	for (i=0;i<npolys;i++)
		arb_poly_init(f[i]);

	for (u=0,phi=0;u<order;u++)
		if (gcd(u,order) == 1) phi++;
	for (j=0;j<ncoeffs;j++)
		res[j] = npolys;
	for (u=0;u<order;u++)
	if (gcd(u,order) == 1) {
		acb_set_si(z,-u);
		acb_div_si(z,z,order,prec);
		acb_exp_pi_i(z,z,prec);
		for (i=0;i<npolys;i++) {
			twist_poly(f[i],g,poly[i],u,2*order);
			if (arb_poly_gcd(f[i],f[i],g,fmpz_poly_degree(poly[i])/phi, prec) < 0) {
				retval = -1;
				goto cleanup;
			}
		}
		for (j=0;j<ncoeffs;j++)
		if (phase(a[j],order) == u) {
			acb_mul(x,a[j],z,prec);
			for (i=0;i<npolys;i++) {
				s = arb_poly_sign_change(f[i],acb_realref(x),prec);
				if (!s || s < 0 && res[j] < npolys) {
					retval = -1;
					goto cleanup;
				}
				if (s < 0) res[j] = i;
			}
#if 1
printf("%d ",res[j]);
acb_printd(a[j],15);
printf("\n");
#endif
		}
	}

	for (j=0;j<ncoeffs;j++)
		if (res[j] == npolys) retval = -1;

cleanup:
	for (i=0;i<npolys;i++)
		arb_poly_clear(f[i]);
	free(f);
	return retval;
}

int main(void) {
	int i,res[1152];
	char *s,buf[1024];
	FILE *fp;
	fmpz_poly_t poly[2];
	acb_t a[1152];
	arb_t e,t;

	fp = fopen("381.4.13.T2","r");
	for (i=0;i<2;i++) {
		fmpz_poly_init(poly[i]);
		fmpz_poly_fread_pretty(fp,poly[i],&s);
	}
	fclose(fp);

	arb_init(e);
	arb_one(e);
	arb_mul_2exp_si(e,e,-190);
	arb_init(t);
	arb_neg(t,e);
	arb_union(e,e,t,200);

	fp = fopen("coefficients","r");
	for (i=0;i<1152;i++) {
		fgets(buf,sizeof(buf),fp);
		acb_init(a[i]);
		s = strtok(buf," ");
		arb_set_str(t,s,200);
		arb_add(acb_realref(a[i]),t,e,200);
		s = strtok(NULL,"\n");
		arb_set_str(t,s,200);
		arb_add(acb_imagref(a[i]),t,e,200);
	}
	fclose(fp);

	printf("%d\n",bober(res,poly,2,a,1152,7,20000));
	return 0;
}
