//========================================================================
// interpolate.hh
//
// A set of simple interpolating functions
//
// I. Volobouev
// March 2009
//========================================================================

#ifndef FFTJET_INTERPOLATE_HH_
#define FFTJET_INTERPOLATE_HH_

#include <cmath>

namespace fftjet {
    inline double interpolate_quadratic(
        const double xmin, const double xmax,
        const double v0, const double v1, const double v2,
        const double x)
    {
        const double h = (xmax - xmin)/2.0;
        const double d1 = (v2 - v0)/(xmax - xmin);
        const double d2 = (v2 + v0 - 2.0*v1)/h/h;
        const double dx = x - (xmax + xmin)/2.0;
        return v1 + dx*(d1 + d2*dx/2.0);
    }

    inline double lin_interpolate_1d(const double x0, const double x1,
                                     const double z0, const double z1,
                                     const double x)
    {
        return z0 + (z1 - z0)*((x - x0)/(x1 - x0));
    }

    inline double log_interpolate_1d(const double x0, const double x1,
                                     const double z0, const double z1,
                                     const double x)
    {
        return z0 + (z1 - z0)*(log(x/x0)/log(x1/x0));
    }
}

#endif // FFTJET_INTERPOLATE_HH_
