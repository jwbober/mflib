#include <iostream>
#include "characters.h"

double * bernoulli;

const double log_2Pi = 1.83787706640935;

using namespace std;



// "stolen" from lcalc...

//compute the nth derivative of log(GAMMA(z))
//n=0 gives log_GAMMA(z).
//Only bothers to handle Re(z)<0 if n=0.
complex<double> log_gamma (complex<double> z, int n = 0)
{
    int M;
    complex<double> log_G, r, r2, y;
    double xx = z.real();
    double yy = z.imag();

    if(abs(z) < 1e-10) return 0.0/0.0;
    //if(xx<0) xx=-xx;

    //else if Re(z)<0 return log(Gamma(z+1)/z)
    if(xx < 0){
       if(n == 0){
         return(log_gamma(z+1.0)-log(z));
       }
       else{
           cout << "error in log_GAMMA: derivative called with Re(z)<0" << endl;
       }
    }

    double x;
    int i,m;

    //assume the remainder stopping at the mth term is roughly bernoulli(2m)/(z+M)^(2m).
    //Now bernoulli(2m) is roughly (2m)!/(2Pi)^(2m). So remainder is more or less
    //(2m/(2ePi(z+M))^(2m). Setting 2m = Digits, we see that we need m/(ePi(z+M))
    //to be roughly 1/10 in size, i.e. |z+M| should be roughly 10*m/(ePi)=10*Digits/(2ePi).
    //squaring, gives us how large M should be.

    //n==0 leads to |z+M| >10*2m/(2*Pi*e) with 2m=Digits. However, when we differentiate
    //we end up with extra powers of (z+M) in the denominators, but extra (2m-1)...(2m+n-2)
    //in the numerators. So assume further that |z+M| > 2m+n = Digits+n

    if(n == 0){
        //.343 = 100./(4*Pi*Pi*exp(2.))
        if((xx*xx+yy*yy)> .343*15*15) M=0;
        else M=int(ceil(sqrt((15*15*.343-yy*yy))-xx+1));
    }
    else{
        if((xx*xx+yy*yy)> .343*(15+n)*(15+n)) M=0;
        else M=int(ceil(sqrt(((15+n)*(15+n)*.343-yy*yy))-xx+1));
    }

    if(n==0){
       log_G=log_2Pi/2+(z+(double)M-.5)*log(z+(double)M)-(z+(double)M);
    }
    else if(n==1)
       log_G=log(z+(double)M)-.5/(z+(double)M);
    else{
       r = 1.0;
       for(m=1;m<=n-1;m++){
          r=-r*(double)m/(z+(double)M);
       }
       log_G=log_G-r/(n-1.0)-.5*r/(z+(double)M);
    }

    r = 1.0;
    for(m=1;m<=n;m++){
       r=-r*(double)m/(z+(double)M);
    }
    r=r/(z+(double)M);

    r2=1.0/((z+(double)M)*(z+(double)M));
    m=2;
    x=abs(r);
    do{
        y=bernoulli[m]*r;
        log_G+=y/(m*(m-1.0));

        r*=(m+n-1.0)*((double)m+(double)n)*r2/((m-1.0)*m);
        m+=2;
    }while(m <= 15 && abs(y)*x > 1e-10);
    //}while(m<=DIGITS);

    for(m=0;m<=M-1;m++){
       if(n==0){
           log_G-=log(z+(double)m); //XXX might be faster to multiply and take one log,
                                 //but careful about overflow errors

       }
       else{
           r=double(1);
           for(i=1;i<=n;i++){
              r = -r*(double)i/(z+(double)m);
           }
           log_G += r/(double)n;
       }
    }

    return log_G;
}

int main(int argc, char ** argv) {
    if(argc < 3) {
        cout << "usage: dirichletL q n" << endl;
        return 0;
    }
    int q = atoi(argv[1]);
    int n = atoi(argv[2]);
    if(GCD(q, n) != 1) {
        cout << n << " is not coprime to " << q << endl;
        return 0;
    }

    bernoulli = new double[16];
    bernoulli[0] = 1;

    for(int k = 1; k <= 15; k++){
        double r = k + 1;
        double x = 0.;
        for(int j = 1; j <= k; j++){
            r = r*(k + 1 -j)*1./(j + 1);
            x = x - bernoulli[k - j]*r;
        }
        bernoulli[k] = x/(k + 1);
    }

    DirichletGroup G(q);
    DirichletCharacter chi = G.character(n);

    if(!chi.is_primitive()) {
        cout << "character number " << n << " is not a primitive character mod " << q << endl;
        return 0;
    }

    for(int k = 0; k < 10; k++) {
        complex<double> s(.5, k);
        cout << exp(log_gamma(s)) << endl;
    }
    return 0;
}
