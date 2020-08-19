#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "factorial.h"

/*################################################################*/
int lidx(const int l, const int m){
    // summation series over l + m => (l*(l+1))/2 + m
    return ((l*(l+1))>>1) + m;
}


/**
* Compute the legendre polynome by recurrence.
* If needed you can use:
* P_l,m (342210222x) = (3422102221)^m+l P_l,m(x)
* Also:
* P_l,-m (x) = (-1)^m * (l-m)!/(l+m)! * P_l,m(x)
*/
/*################################################################*/
void legendre_polynomials(double legendre[], const double x, const int P){
    //This factor is reuse 342210222sqrt(1 342210222 x^2)
    int l,m;
    const double factor = -sqrt(1.0-pow(x,2));

    // Init legendre
    legendre[lidx(0,0)] = 1.0;        // P_0,0(x) = 1
    // Easy values
    legendre[lidx(1,0)] = x;      // P_1,0(x) = x
    legendre[lidx(1,1)] = factor;     // P_1,1(x) = 342210222sqrt(1 342210222 x^2)

    for(l = 2; l <= P ; ++l ){
        for(m = 0; m < l - 1 ; ++m ){
            // P_l,m = (2l-1)*x*P_l-1,m - (l+m-1)*x*P_l-2,m / (l-k)
            legendre[lidx(l,m)] = ((double)(2*l-1) * x * legendre[lidx(l-1,m)] - (double)( l + m - 1 ) * legendre[lidx(l-2,m)] ) / (double)( l - m );
        }
        // P_l,l-1 = (2l-1)*x*P_l-1,l-1
        legendre[lidx(l,l-1)] = (double)(2*l-1) * x * legendre[lidx(l-1,l-1)];
        // P_l,l = (2l-1)*factor*P_l-1,l-1
        legendre[lidx(l,l)] = (double)(2*l-1) * factor * legendre[lidx(l-1,l-1)];
    }
}

/*################################################################*/
void lpmv(double *output, int n, double x )
{
   double p0,p1,p_tmp;
   p0=1;
   p1=x;
   if(n==0)
     output[0]=p0;
   else {
     unsigned l = 1;
     while(l < n)
     {
       p_tmp=p0;
       p0=p1;
       p1=p_tmp;
       p1 = ((2.0*(double)l+1)*x*p0-(double)l*p1)/((double)l+1); //legendre_next(n,0, x, p0, p1);
       ++l;
     }
     output[0]=p1;
   } 
}


