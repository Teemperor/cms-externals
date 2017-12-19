//=========================================================================
// SpecialFunctions.hh
//
// A few mathematical special functions needed by this package.
// Note that these are not optimized in any way, and may be slow.
// They are not used in speed-critical parts of the package code.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_SPECIALFUNCTIONS_HH_
#define FFTJET_SPECIALFUNCTIONS_HH_

namespace fftjet {
    // Inverse cumulative distribition function for 1d Gaussian
    double inverseGaussCdf(double cdf);

    // Regularized incomplete beta function
    double incompleteBeta(double a, double b, double x);

    // Inverse regularized incomplete beta function
    double inverseIncompleteBeta(double a, double b, double x);

    // The gamma function for positive real arguments
    double Gamma(double x);

    // Incomplete gamma ratio and its complemented function
    double incompleteGamma(double a, double x);
    double incompleteGammaC(double a, double x);

    // Inverse incomplete gamma ratio
    double inverseIncompleteGamma(double a, double x);
}

#endif // FFTJET_SPECIALFUNCTIONS_HH_
