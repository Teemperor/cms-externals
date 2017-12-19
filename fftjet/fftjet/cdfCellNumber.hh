//=========================================================================
// cdfCellNumber.hh
//
// Utility functions for building inverse cumulative distributions
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_CDFCELLNUMBER_HH_
#define FFTJET_CDFCELLNUMBER_HH_

namespace fftjet {
    // A utility function for looking up the cumulative distribution
    // function cell number given the cdf value
    unsigned cdfCellNumber(double cdf, const double* cdfdata,
                           unsigned ncells, unsigned stride=1);

    // Look for the solution of the equation a x^2 + b x == c
    // on the x interval from 0 to 1
    double invertQuadraticCdf(double a, double b, double c);
}

#endif // FFTJET_CDFCELLNUMBER_HH_
