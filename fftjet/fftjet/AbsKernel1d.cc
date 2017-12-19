#include <cmath>

#include "fftjet/AbsKernel1d.hh"

namespace fftjet {
    double AbsKernel1d::intervalAverage(const double x, const double scale,
                                        const double dx) const
    {
        // 3-point Gaussian integration. Precise for 5th order polynomials.
        const double legendreRoot = sqrt(0.6);
        const double delta = legendreRoot*dx/2.0;
        return (5.0*operator()(x - delta, scale) +
                8.0*operator()(x, scale) +
                5.0*operator()(x + delta, scale))/18.0;
    }
}
