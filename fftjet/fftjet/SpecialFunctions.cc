// Most of the code below is cherry-picked from the Cephes library
// written by Stephen L. Moshier. This code is not thread-safe!

#include <cassert>
#include <cmath>
#include <cfloat>

#include "fftjet/SpecialFunctions.hh"

#define SQTPI 2.50662827463100050242
#define LOGPI 1.14472988584940017414

#define MAXGAM 171.624376956302725
#define MAXSTIR 143.01608
#define MAXLOG 709.782712893383996732224
#define MINLOG (-708.396418532264106224)
#define MACHEP 1.2e-16
#define MAXLGM 2.556348e305
#define MAXNUM 1.79769313486231570815e308
#define LS2PI 0.91893853320467274178

#define GAUSS_MAX_SIGMA 40.0

static double igam(double a, double x);
static double invgauss(double x);

static double polevl(const double x, const double *p, int i)
{
    double ans = *p++;
    do {
	ans = ans * x + *p++;
    } while (--i);
    return ans;
}

static double p1evl(double x, const double *p, const int N)
{
    double ans = x + *p++;
    int i = N - 1;
    do {
	ans = ans * x  + *p++;
    } while (--i);
    return ans;
}

/* Gamma function computed by Stirling's formula.
 * The polynomial STIR is valid for 33 <= x <= 172.
 */
static double stirf(const double x)
{
    static const double STIR[5] = {
        7.87311395793093628397E-4,
        -2.29549961613378126380E-4,
        -2.68132617805781232825E-3,
        3.47222221605458667310E-3,
        8.33333333333482257126E-2,
    };

    double y, w;
    w = 1.0/x;
    w = 1.0 + w * polevl( w, STIR, 4 );
    y = exp(x);
    if( x > MAXSTIR )
    {
        // Avoid overflow in pow()
	const double v = pow( x, 0.5 * x - 0.25 );
	y = v * (v / y);
    }
    else
    {
	y = pow( x, x - 0.5 ) / y;
    }
    return SQTPI * y * w;
}

// Logarithm of the gamma function
static double lgam(double x)
{
    static int sgngam = 1;

    double p, q, u, w, z;
    int i;

    static const double A[] = {
        8.11614167470508450300E-4,
        -5.95061904284301438324E-4,
        7.93650340457716943945E-4,
        -2.77777777730099687205E-3,
        8.33333333333331927722E-2
    };
    static const double B[] = {
        -1.37825152569120859100E3,
        -3.88016315134637840924E4,
        -3.31612992738871184744E5,
        -1.16237097492762307383E6,
        -1.72173700820839662146E6,
        -8.53555664245765465627E5
    };
    static const double C[] = {
        -3.51815701436523470549E2,
        -1.70642106651881159223E4,
        -2.20528590553854454839E5,
        -1.13933444367982507207E6,
        -2.53252307177582951285E6,
        -2.01889141433532773231E6
    };

    sgngam = 1;

    if( x < -34.0 )
    {
	q = -x;
	w = lgam(q); /* note this modifies sgngam! */
	p = floor(q);
	if( p == q )
            goto loverf;
	i = static_cast<int>(p);
	if( (i & 1) == 0 )
            sgngam = -1;
	else
            sgngam = 1;
	z = q - p;
	if( z > 0.5 )
        {
            p += 1.0;
            z = p - q;
        }
	z = q * sin( M_PI * z );
	if( z == 0.0 )
            goto loverf;
	// z = log(PI) - log( z ) - w;
	z = LOGPI - log( z ) - w;
	return( z );
    }

    if( x < 13.0 )
    {
	z = 1.0;
	p = 0.0;
	u = x;
	while( u >= 3.0 )
        {
            p -= 1.0;
            u = x + p;
            z *= u;
        }
	while( u < 2.0 )
        {
            if( u == 0.0 )
                goto loverf;
            z /= u;
            p += 1.0;
            u = x + p;
        }
	if( z < 0.0 )
        {
            sgngam = -1;
            z = -z;
        }
	else
            sgngam = 1;
	if( u == 2.0 )
            return( log(z) );
	p -= 2.0;
	x = x + p;
	p = x * polevl( x, B, 5 ) / p1evl( x, C, 6);
	return( log(z) + p );
    }

    if( x > MAXLGM )
    {
    loverf:
	assert(!"Overflow in lgam");
	return( sgngam * MAXNUM );
    }

    q = ( x - 0.5 ) * log(x) - x + LS2PI;
    if( x > 1.0e8 )
	return( q );

    p = 1.0/(x*x);
    if( x >= 1000.0 )
	q += ((   7.9365079365079365079365e-4 * p
                  - 2.7777777777777777777778e-3) *p
              + 0.0833333333333333333333) / x;
    else
	q += polevl( p, A, 4 ) / x;
    return( q );
}

// Complementary incomplete gamma
static double igamc(double a, double x)
{
    const double big = 4.503599627370496e15;
    const double biginv =  2.22044604925031308085e-16;

    double ans, ax, c, yc, r, t, y, z;
    double pk, pkm1, pkm2, qk, qkm1, qkm2;

    if( (x <= 0) || ( a <= 0) )
	return( 1.0 );

    if( (x < 1.0) || (x < a) )
	return( 1.0 - igam(a,x) );

    ax = a * log(x) - x - lgam(a);
    if( ax < -MAXLOG )
    {
	// mtherr( "igamc", UNDERFLOW );
        assert(!"Underflow in igamc");
	return( 0.0 );
    }
    ax = exp(ax);

    // continued fraction
    y = 1.0 - a;
    z = x + y + 1.0;
    c = 0.0;
    pkm2 = 1.0;
    qkm2 = x;
    pkm1 = x + 1.0;
    qkm1 = z * x;
    ans = pkm1/qkm1;

    do
    {
	c += 1.0;
	y += 1.0;
	z += 2.0;
	yc = y * c;
	pk = pkm1 * z  -  pkm2 * yc;
	qk = qkm1 * z  -  qkm2 * yc;
	if( qk != 0 )
        {
            r = pk/qk;
            t = fabs( (ans - r)/r );
            ans = r;
        }
	else
            t = 1.0;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;
	if( fabs(pk) > big )
        {
            pkm2 *= biginv;
            pkm1 *= biginv;
            qkm2 *= biginv;
            qkm1 *= biginv;
        }
    }
    while( t > MACHEP );

    return( ans * ax );
}

// Incomplete gamma
static double igam(double a, double x)
{
    double ans, ax, c, r;

    if( (x <= 0) || ( a <= 0) )
	return( 0.0 );

    if( (x > 1.0) && (x > a ) )
	return( 1.0 - igamc(a,x) );

    // Compute  x**a * exp(-x) / gamma(a)
    ax = a * log(x) - x - lgam(a);
    if( ax < -MAXLOG )
    {
	// mtherr( "igam", UNDERFLOW );
        assert(!"Underflow in igam");
	return( 0.0 );
    }
    ax = exp(ax);

    // power series
    r = a;
    c = 1.0;
    ans = 1.0;

    do
    {
	r += 1.0;
	c *= x/r;
	ans += c;
    }
    while( c/ans > MACHEP );

    return( ans * ax/a );
}

// Inverse of complemented incomplete gamma integral.
// Work only for y0 <= 0.5.
static double igami(double a, double y0)
{
    double x0, x1, x, yl, yh, y, d, lgm, dithresh;
    int i, dir;

    assert(y0 <= 0.5);

    /* bound the solution */
    x0 = MAXNUM;
    yl = 0;
    x1 = 0;
    yh = 1.0;
    dithresh = 5.0 * MACHEP;

    /* approximation to inverse function */
    d = 1.0/(9.0*a);
    y = ( 1.0 - d - invgauss(y0) * sqrt(d) );
    x = a * y * y * y;

    lgm = lgam(a);

    for( i=0; i<10; i++ )
    {
	if( x > x0 || x < x1 )
            goto ihalve;
	y = igamc(a,x);
	if( y < yl || y > yh )
            goto ihalve;
	if( y < y0 )
        {
            x0 = x;
            yl = y;
        }
	else
        {
            x1 = x;
            yh = y;
        }
        /* compute the derivative of the function at this point */
	d = (a - 1.0) * log(x) - x - lgm;
	if( d < -MAXLOG )
            goto ihalve;
	d = -exp(d);
        /* compute the step to the next approximation of x */
	d = (y - y0)/d;
	if( fabs(d/x) < MACHEP )
            goto done;
	x = x - d;
    }

    /* Resort to interval halving if Newton iteration did not converge. */
ihalve:

    d = 0.0625;
    if( x0 == MAXNUM )
    {
	if( x <= 0.0 )
            x = 1.0;
	while( x0 == MAXNUM )
        {
            x = (1.0 + d) * x;
            y = igamc( a, x );
            if( y < y0 )
            {
                x0 = x;
                yl = y;
                break;
            }
            d = d + d;
        }
    }
    d = 0.5;
    dir = 0;

    for( i=0; i<400; i++ )
    {
	x = x1  +  d * (x0 - x1);
	y = igamc( a, x );
	lgm = (x0 - x1)/(x1 + x0);
	if( fabs(lgm) < dithresh )
            break;
	lgm = (y - y0)/y0;
	if( fabs(lgm) < dithresh )
            break;
	if( x <= 0.0 )
            break;
	if( y >= y0 )
        {
            x1 = x;
            yh = y;
            if( dir < 0 )
            {
                dir = 0;
                d = 0.5;
            }
            else if( dir > 1 )
                d = 0.5 * d + 0.5; 
            else
                d = (y0 - yl)/(yh - yl);
            dir += 1;
        }
	else
        {
            x0 = x;
            yl = y;
            if( dir > 0 )
            {
                dir = 0;
                d = 0.5;
            }
            else if( dir < -1 )
                d = 0.5 * d;
            else
                d = (y0 - yl)/(yh - yl);
            dir -= 1;
        }
    }
    if( x == 0.0 )
	// mtherr( "igami", UNDERFLOW );
        assert(!"Underflow in igami");

done:
    return( x );
}

static double invgauss(const double P)
{
    assert(P > 0.0);
    assert(P < 1.0);

    // Translated from PPND16 algorithm of Wichura (originally in Fortran).
    // This was not taken from Cephes.
    static const double ZERO = 0., ONE = 1., HALF = 0.5, 
        SPLIT1 = 0.425, SPLIT2 = 5., 
        CONST1 = 0.180625, CONST2 = 1.6;

    static const double A0 = 3.3871328727963666080, 
        A1 = 1.3314166789178437745E+2, 
        A2 = 1.9715909503065514427E+3, 
        A3 = 1.3731693765509461125E+4, 
        A4 = 4.5921953931549871457E+4, 
        A5 = 6.7265770927008700853E+4, 
        A6 = 3.3430575583588128105E+4, 
        A7 = 2.5090809287301226727E+3, 
        B1 = 4.2313330701600911252E+1, 
        B2 = 6.8718700749205790830E+2, 
        B3 = 5.3941960214247511077E+3, 
        B4 = 2.1213794301586595867E+4, 
        B5 = 3.9307895800092710610E+4, 
        B6 = 2.8729085735721942674E+4, 
        B7 = 5.2264952788528545610E+3;

    static const double C0 = 1.42343711074968357734, 
        C1 = 4.63033784615654529590, 
        C2 = 5.76949722146069140550, 
        C3 = 3.64784832476320460504, 
        C4 = 1.27045825245236838258, 
        C5 = 2.41780725177450611770E-1, 
        C6 = 2.27238449892691845833E-2, 
        C7 = 7.74545014278341407640E-4, 
        D1 = 2.05319162663775882187, 
        D2 = 1.67638483018380384940, 
        D3 = 6.89767334985100004550E-1, 
        D4 = 1.48103976427480074590E-1, 
        D5 = 1.51986665636164571966E-2, 
        D6 = 5.47593808499534494600E-4, 
        D7 = 1.05075007164441684324E-9;

    static const double E0 = 6.65790464350110377720, 
        E1 = 5.46378491116411436990, 
        E2 = 1.78482653991729133580, 
        E3 = 2.96560571828504891230E-1, 
        E4 = 2.65321895265761230930E-2, 
        E5 = 1.24266094738807843860E-3, 
        E6 = 2.71155556874348757815E-5, 
        E7 = 2.01033439929228813265E-7, 
        F1 = 5.99832206555887937690E-1, 
        F2 = 1.36929880922735805310E-1, 
        F3 = 1.48753612908506148525E-2, 
        F4 = 7.86869131145613259100E-4, 
        F5 = 1.84631831751005468180E-5, 
        F6 = 1.42151175831644588870E-7, 
        F7 = 2.04426310338993978564E-15;

    const double Q = P - HALF;

    double R, PPND16;
    if (fabs(Q) <= SPLIT1)
    {
        R = CONST1 - Q * Q;
        PPND16 = Q * (((((((A7 * R + A6) * R + A5) * R + A4) * R + A3) 
                        * R + A2) * R + A1) * R + A0) / 
            (((((((B7 * R + B6) * R + B5) * R + B4) * R + B3) 
               * R + B2) * R + B1) * R + ONE);
    }
    else
    {
        if (Q < ZERO)
            R = P;
        else
            R = ONE - P;
        R = sqrt(-log(R));
        if (R <= SPLIT2)
        {
            R = R - CONST2;
            PPND16 = (((((((C7 * R + C6) * R + C5) * R + C4) * R + C3) 
                        * R + C2) * R + C1) * R + C0) / 
                (((((((D7 * R + D6) * R + D5) * R + D4) * R + D3) 
                   * R + D2) * R + D1) * R + ONE);
        }
        else
        {
            R = R - SPLIT2 ;
            PPND16 = (((((((E7 * R + E6) * R + E5) * R + E4) * R + E3) 
                        * R + E2) * R + E1) * R + E0) / 
                (((((((F7 * R + F6) * R + F5) * R + F4) * R + F3) 
                   * R + F2) * R + F1) * R + ONE);
        }
        if (Q < ZERO)
            PPND16 = -PPND16;
    }
    return PPND16;
}

/* Continued fraction expansion #1
 * for incomplete beta integral
 */
static double incbcf( double a, double b, double x )
{
    static const double big = 4.503599627370496e15;
    static const double biginv = 2.22044604925031308085e-16;

    double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
    double k1, k2, k3, k4, k5, k6, k7, k8;
    double r, t, ans, thresh;
    int n;

    k1 = a;
    k2 = a + b;
    k3 = a;
    k4 = a + 1.0;
    k5 = 1.0;
    k6 = b - 1.0;
    k7 = k4;
    k8 = a + 2.0;

    pkm2 = 0.0;
    qkm2 = 1.0;
    pkm1 = 1.0;
    qkm1 = 1.0;
    ans = 1.0;
    r = 1.0;
    n = 0;
    thresh = 3.0 * MACHEP;
    do
    {
	
	xk = -( x * k1 * k2 )/( k3 * k4 );
	pk = pkm1 +  pkm2 * xk;
	qk = qkm1 +  qkm2 * xk;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;

	xk = ( x * k5 * k6 )/( k7 * k8 );
	pk = pkm1 +  pkm2 * xk;
	qk = qkm1 +  qkm2 * xk;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;

	if( qk != 0 )
            r = pk/qk;
	if( r != 0 )
        {
            t = fabs( (ans - r)/r );
            ans = r;
        }
	else
            t = 1.0;

	if( t < thresh )
            goto cdone;

	k1 += 1.0;
	k2 += 1.0;
	k3 += 2.0;
	k4 += 2.0;
	k5 += 1.0;
	k6 -= 1.0;
	k7 += 2.0;
	k8 += 2.0;

	if( (fabs(qk) + fabs(pk)) > big )
        {
            pkm2 *= biginv;
            pkm1 *= biginv;
            qkm2 *= biginv;
            qkm1 *= biginv;
        }
	if( (fabs(qk) < biginv) || (fabs(pk) < biginv) )
        {
            pkm2 *= big;
            pkm1 *= big;
            qkm2 *= big;
            qkm1 *= big;
        }
    }
    while( ++n < 300 );

cdone:
    return(ans);
}

/* Continued fraction expansion #2
 * for incomplete beta integral
 */
static double incbd( double a, double b, double x )
{
    static const double big = 4.503599627370496e15;
    static const double biginv = 2.22044604925031308085e-16;

    double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
    double k1, k2, k3, k4, k5, k6, k7, k8;
    double r, t, ans, z, thresh;
    int n;

    k1 = a;
    k2 = b - 1.0;
    k3 = a;
    k4 = a + 1.0;
    k5 = 1.0;
    k6 = a + b;
    k7 = a + 1.0;;
    k8 = a + 2.0;

    pkm2 = 0.0;
    qkm2 = 1.0;
    pkm1 = 1.0;
    qkm1 = 1.0;
    z = x / (1.0-x);
    ans = 1.0;
    r = 1.0;
    n = 0;
    thresh = 3.0 * MACHEP;
    do
    {
	
	xk = -( z * k1 * k2 )/( k3 * k4 );
	pk = pkm1 +  pkm2 * xk;
	qk = qkm1 +  qkm2 * xk;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;

	xk = ( z * k5 * k6 )/( k7 * k8 );
	pk = pkm1 +  pkm2 * xk;
	qk = qkm1 +  qkm2 * xk;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;

	if( qk != 0 )
            r = pk/qk;
	if( r != 0 )
        {
            t = fabs( (ans - r)/r );
            ans = r;
        }
	else
            t = 1.0;

	if( t < thresh )
            goto cdone;

	k1 += 1.0;
	k2 -= 1.0;
	k3 += 2.0;
	k4 += 2.0;
	k5 += 1.0;
	k6 += 1.0;
	k7 += 2.0;
	k8 += 2.0;

	if( (fabs(qk) + fabs(pk)) > big )
        {
            pkm2 *= biginv;
            pkm1 *= biginv;
            qkm2 *= biginv;
            qkm1 *= biginv;
        }
	if( (fabs(qk) < biginv) || (fabs(pk) < biginv) )
        {
            pkm2 *= big;
            pkm1 *= big;
            qkm2 *= big;
            qkm1 *= big;
        }
    }
    while( ++n < 300 );
cdone:
    return(ans);
}

/* Power series for incomplete beta integral.
   Use when b*x is small and x not too close to 1.  */
static double pseries( double a, double b, double x )
{
    double s, t, u, v, n, t1, z, ai;

    ai = 1.0 / a;
    u = (1.0 - b) * x;
    v = u / (a + 1.0);
    t1 = v;
    t = u;
    n = 2.0;
    s = 0.0;
    z = MACHEP * ai;
    while( fabs(v) > z )
    {
	u = (n - b) * x / n;
	t *= u;
	v = t / (a + n);
	s += v; 
	n += 1.0;
    }
    s += t1;
    s += ai;

    u = a * log(x);
    if( (a+b) < MAXGAM && fabs(u) < MAXLOG )
    {
	t = fftjet::Gamma(a+b)/(fftjet::Gamma(a)*fftjet::Gamma(b));
	s = s * t * pow(x,a);
    }
    else
    {
	t = lgam(a+b) - lgam(a) - lgam(b) + u + log(s);
	if( t < MINLOG )
            s = 0.0;
	else
            s = exp(t);
    }
    return(s);
}

/* Incomplete beta function */
static double incbet( double aa, double bb, double xx )
{
    double a, b, t, x, xc, w, y;
    int flag;

    if( aa <= 0.0 || bb <= 0.0 )
	goto domerr;

    if( (xx <= 0.0) || ( xx >= 1.0) )
    {
	if( xx == 0.0 )
            return(0.0);
	if( xx == 1.0 )
            return( 1.0 );
    domerr:
	// mtherr( "incbet", DOMAIN );
        assert(!"Domain error in incbet");
	return( 0.0 );
    }

    flag = 0;
    if( (bb * xx) <= 1.0 && xx <= 0.95)
    {
	t = pseries(aa, bb, xx);
        goto done;
    }

    w = 1.0 - xx;

    /* Reverse a and b if x is greater than the mean. */
    if( xx > (aa/(aa+bb)) )
    {
	flag = 1;
	a = bb;
	b = aa;
	xc = xx;
	x = w;
    }
    else
    {
	a = aa;
	b = bb;
	xc = w;
	x = xx;
    }

    if( flag == 1 && (b * x) <= 1.0 && x <= 0.95)
    {
	t = pseries(a, b, x);
	goto done;
    }

    /* Choose expansion for better convergence. */
    y = x * (a+b-2.0) - (a-1.0);
    if( y < 0.0 )
	w = incbcf( a, b, x );
    else
	w = incbd( a, b, x ) / xc;

   /* Multiply w by the factor
   a      b   _             _     _
   x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   */

    y = a * log(x);
    t = b * log(xc);
    if( (a+b) < MAXGAM && fabs(y) < MAXLOG && fabs(t) < MAXLOG )
    {
	t = pow(xc,b);
	t *= pow(x,a);
	t /= a;
	t *= w;
	t *= fftjet::Gamma(a+b) / (fftjet::Gamma(a) * fftjet::Gamma(b));
	goto done;
    }

    /* Resort to logarithms.  */
    y += t + lgam(a+b) - lgam(a) - lgam(b);
    y += log(w/a);
    if( y < MINLOG )
	t = 0.0;
    else
	t = exp(y);

done:

    if( flag == 1 )
    {
	if( t <= MACHEP )
            t = 1.0 - MACHEP;
	else
            t = 1.0 - t;
    }
    return( t );
}

/* Inverse incomplete beta function */
static double incbi( double aa, double bb, double yy0 )
{
    double a, b, y0, d, y, x, x0, x1, lgm, yp, di, dithresh, yl, yh, xt;
    int i, rflg, dir, nflg;

    i = 0;
    if( yy0 <= 0 )
	return(0.0);
    if( yy0 >= 1.0 )
	return(1.0);
    x0 = 0.0;
    yl = 0.0;
    x1 = 1.0;
    yh = 1.0;
    nflg = 0;

    if( aa <= 1.0 || bb <= 1.0 )
    {
	dithresh = 1.0e-6;
	rflg = 0;
	a = aa;
	b = bb;
	y0 = yy0;
	x = a/(a+b);
	y = incbet( a, b, x );
	goto ihalve;
    }
    else
    {
	dithresh = 1.0e-4;
    }
    /* approximation to inverse function */

    yp = -fftjet::inverseGaussCdf(yy0);

    if( yy0 > 0.5 )
    {
	rflg = 1;
	a = bb;
	b = aa;
	y0 = 1.0 - yy0;
	yp = -yp;
    }
    else
    {
	rflg = 0;
	a = aa;
	b = bb;
	y0 = yy0;
    }

    lgm = (yp * yp - 3.0)/6.0;
    x = 2.0/( 1.0/(2.0*a-1.0)  +  1.0/(2.0*b-1.0) );
    d = yp * sqrt( x + lgm ) / x
	- ( 1.0/(2.0*b-1.0) - 1.0/(2.0*a-1.0) )
	* (lgm + 5.0/6.0 - 2.0/(3.0*x));
    d = 2.0 * d;
    if( d < MINLOG )
    {
	x = 1.0;
	goto under;
    }
    x = a/( a + b * exp(d) );
    y = incbet( a, b, x );
    yp = (y - y0)/y0;
    if( fabs(yp) < 0.2 )
	goto newt;

    /* Resort to interval halving if not close enough. */
ihalve:

    dir = 0;
    di = 0.5;
    for( i=0; i<100; i++ )
    {
	if( i != 0 )
        {
            x = x0  +  di * (x1 - x0);
            if( x == 1.0 )
                x = 1.0 - MACHEP;
            if( x == 0.0 )
            {
                di = 0.5;
                x = x0  +  di * (x1 - x0);
                if( x == 0.0 )
                    goto under;
            }
            y = incbet( a, b, x );
            yp = (x1 - x0)/(x1 + x0);
            if( fabs(yp) < dithresh )
                goto newt;
            yp = (y-y0)/y0;
            if( fabs(yp) < dithresh )
                goto newt;
        }
	if( y < y0 )
        {
            x0 = x;
            yl = y;
            if( dir < 0 )
            {
                dir = 0;
                di = 0.5;
            }
            else if( dir > 3 )
                di = 1.0 - (1.0 - di) * (1.0 - di);
            else if( dir > 1 )
                di = 0.5 * di + 0.5; 
            else
                di = (y0 - y)/(yh - yl);
            dir += 1;
            if( x0 > 0.75 )
            {
                if( rflg == 1 )
                {
                    rflg = 0;
                    a = aa;
                    b = bb;
                    y0 = yy0;
                }
                else
                {
                    rflg = 1;
                    a = bb;
                    b = aa;
                    y0 = 1.0 - yy0;
                }
                x = 1.0 - x;
                y = incbet( a, b, x );
                x0 = 0.0;
                yl = 0.0;
                x1 = 1.0;
                yh = 1.0;
                goto ihalve;
            }
        }
	else
        {
            x1 = x;
            if( rflg == 1 && x1 < MACHEP )
            {
                x = 0.0;
                goto done;
            }
            yh = y;
            if( dir > 0 )
            {
                dir = 0;
                di = 0.5;
            }
            else if( dir < -3 )
                di = di * di;
            else if( dir < -1 )
                di = 0.5 * di;
            else
                di = (y - y0)/(yh - yl);
            dir -= 1;
        }
    }

    // Partial loss of precision. Hm... Will ignore for now - igv
    // mtherr( "incbi", PLOSS );
    if( x0 >= 1.0 )
    {
	x = 1.0 - MACHEP;
	goto done;
    }
    if( x <= 0.0 )
    {
    under:
	// mtherr( "incbi", UNDERFLOW );
        assert(!"Underflow in incbi");
	x = 0.0;
	goto done;
    }

newt:

    if( nflg )
	goto done;
    nflg = 1;
    lgm = lgam(a+b) - lgam(a) - lgam(b);

    for( i=0; i<8; i++ )
    {
	/* Compute the function at this point. */
	if( i != 0 )
            y = incbet(a,b,x);
	if( y < yl )
        {
            x = x0;
            y = yl;
        }
	else if( y > yh )
        {
            x = x1;
            y = yh;
        }
	else if( y < y0 )
        {
            x0 = x;
            yl = y;
        }
	else
        {
            x1 = x;
            yh = y;
        }
	if( x == 1.0 || x == 0.0 )
            break;
	/* Compute the derivative of the function at this point. */
	d = (a - 1.0) * log(x) + (b - 1.0) * log(1.0-x) + lgm;
	if( d < MINLOG )
            goto done;
	if( d > MAXLOG )
            break;
	d = exp(d);
	/* Compute the step to the next approximation of x. */
	d = (y - y0)/d;
	xt = x - d;
	if( xt <= x0 )
        {
            y = (x - x0) / (x1 - x0);
            xt = x0 + 0.5 * y * (x - x0);
            if( xt <= 0.0 )
                break;
        }
	if( xt >= x1 )
        {
            y = (x1 - x) / (x1 - x0);
            xt = x1 - 0.5 * y * (x1 - x);
            if( xt >= 1.0 )
                break;
        }
	x = xt;
	if( fabs(d/x) < 128.0 * MACHEP )
            goto done;
    }

    /* Did not converge.  */
    dithresh = 256.0 * MACHEP;
    goto ihalve;

done:

    if( rflg )
    {
	if( x <= MACHEP )
            x = 1.0 - MACHEP;
	else
            x = 1.0 - x;
    }
    return( x );
}

namespace fftjet {
    double inverseGaussCdf(const double x)
    {
        assert(x >= 0.0);
        assert(x <= 1.0);
        if (x == 1.0)
            return GAUSS_MAX_SIGMA;
        else if (x == 0.0)
            return -GAUSS_MAX_SIGMA;
        else
            return invgauss(x);
    }

    double incompleteBeta(const double a, const double b,
                          const double x)
    {
        assert(a > 0.0);
        assert(b > 0.0);
        assert(x >= 0.0);
        assert(x <= 1.0);
        return incbet(a, b, x);
    }

    double inverseIncompleteBeta(const double a, const double b,
                                 const double x)
    {
        assert(a > 0.0);
        assert(b > 0.0);
        assert(x >= 0.0);
        assert(x <= 1.0);
        return incbi(a, b, x);
    }

    double Gamma(double x)
    {
        assert(x > 0.0);

        // The Stirling formula overflows for x >= MAXGAM
        assert(x < MAXGAM);

        static const double P[] = {
            1.60119522476751861407E-4,
            1.19135147006586384913E-3,
            1.04213797561761569935E-2,
            4.76367800457137231464E-2,
            2.07448227648435975150E-1,
            4.94214826801497100753E-1,
            9.99999999999999996796E-1
        };
        static const double Q[] = {
            -2.31581873324120129819E-5,
            5.39605580493303397842E-4,
            -4.45641913851797240494E-3,
            1.18139785222060435552E-2,
            3.58236398605498653373E-2,
            -2.34591795718243348568E-1,
            7.14304917030273074085E-2,
            1.00000000000000000320E0
        };

        double p, z = 1.0, q = x;

        if (q > 33.0)
            return stirf(q);

        while( x >= 3.0 )
	{
            x -= 1.0;
            z *= x;
	}

        while( x < 2.0 )
	{
            if( x < 1.e-9 )
		return( z/((1.0 + 0.5772156649015329 * x) * x) );
            z /= x;
            x += 1.0;
	}

        if( x == 2.0 )
            return(z);

        x -= 2.0;
        p = polevl( x, P, 6 );
        q = polevl( x, Q, 7 );
        return( z * p / q );
    }

    double incompleteGamma(const double a, const double x)
    {
        return igam(a, x);
    }

    double incompleteGammaC(const double a, const double x)
    {
        return igamc(a, x);
    }

    double inverseIncompleteGamma(const double a, const double x)
    {
        if (x >= 0.5)
            return igami(a, 1.0 - x);
        else
        {
            const double targetEps = 8.0*MACHEP;
            double xmin = 0.0;
            double xmax = igami(a, 0.5);
            while ((xmax - xmin)/(xmax + xmin + 1.0) > targetEps)
            {
                const double xtry = (xmax + xmin)/2.0;
                if (igam(a, xtry) > x)
                    xmax = xtry;
                else
                    xmin = xtry;
            }
            return (xmax + xmin)/2.0;
        }
    }
}
