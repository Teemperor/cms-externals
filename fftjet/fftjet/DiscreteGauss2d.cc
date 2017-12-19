#include <cmath>
#include <cassert>

#include "fftjet/DiscreteGauss2d.hh"

namespace fftjet {
    DiscreteGauss2d::DiscreteGauss2d(
        const double sx, const double sy,
        const unsigned nEta, const unsigned nPhi)
        : AbsFrequencyKernel(),
          sEta_(sx),
          sPhi_(sy),
          nEta_(nEta),
          nPhi_(nPhi)
    {
        // Check that the inputs are reasonable
        assert(sEta_ > 0.0);
        assert(sPhi_ > 0.0);
        assert(nEta_);
        assert(nPhi_);
    }

    std::complex<double> DiscreteGauss2d::operator()(
        const int ix, const int iy, const double scale) const
    {
        const double cosum1 = cos(ix*2.0*M_PI/nEta_) - 1.0;
        const double cosvm1 = cos(iy*2.0*M_PI/nPhi_) - 1.0;
        const double detasq = 4.0*M_PI*M_PI/nEta_/nEta_;
        const double dphisq = 4.0*M_PI*M_PI/nPhi_/nPhi_;
        const double sxsq = sEta_*sEta_*scale*scale/detasq;
        const double sysq = sPhi_*sPhi_*scale*scale/dphisq;
        const double fourier = exp(sxsq*cosum1 + sysq*cosvm1);
        return std::complex<double>(fourier, 0.0);
    }
}
