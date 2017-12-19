//=========================================================================
// invertJetResponse.hh
//
// The intended typical use of these functions is figuring out jet
// corrections from the jet response in the form of median (mean, etc)
// ratio between reconstructed and actual jet energy.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_INVERTJETRESPONSE_HH_
#define FFTJET_INVERTJETRESPONSE_HH_

namespace fftjet {
    // The following function solves numerically the equation x*f(x) = y
    // for unknown x. It is assumed that x*f(x) is a monotonously increasing
    // function, x*f(x) -> 0 when x -> 0, f(x) is positive, and that both
    // x and y are non-negative.
    template <class Functor>
    double invertJetResponse(const Functor& f, double y);

    // The following function solves numerically the equation x*f(a,x) = y
    // for unknown x, where "a" is some fixed parameter. It is assumed that
    // x*f(a,x) is a monotonously increasing function of x, x*f(a,x) -> 0 when
    // x -> 0, f(a,x) is positive, and that both x and y are non-negative.
    template <class Functor2d>
    double invertJetResponse2d(const Functor2d& f, double a, double y);
}

#include "fftjet/invertJetResponse.icc"

#endif // FFTJET_INVERTJETRESPONSE_HH_
