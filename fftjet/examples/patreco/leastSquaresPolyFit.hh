//========================================================================
// leastSquaresPolyFit.hh
//
// Least squares polynomial fit to a bunch of points. In the current
// implementation, all the heavy lifting is done by LAPACK.
//
// I. Volobouev
// April 2009
//========================================================================

#ifndef LEASTSQUARESPOLYFIT_HH_
#define LEASTSQUARESPOLYFIT_HH_

// The number of data points "npoints" should be at least polyDegree+1.
// The dimension of the "coeffs" array should be at least polyDegree+1.
// The coefficients are returned in the order of decreasing degree.
// The function returns the sum of squared residuals.
// The "constTermError" will be set to the fit error of the constant
// term which equals to the fit error of the polynomial when x = 0
// (this is what we may need in FFTJet). The "constTermError" pointer
// can be NULL which can save some CPU time in case you do not need
// to estimate this error.
double leastSquaresPolyFit(unsigned polyDegree,
                           const double *x, const double *y,
                           unsigned npoints, double *coeffs,
                           double *constTermError=0);

// A convenience function which returns the value of the polynomial
// using a compatible format for the coefficients
double polyValue(const double *coeffs, unsigned polyDegree, double x);

#endif // LEASTSQUARESPOLYFIT_HH_
