#include <cmath>
#include <cassert>

#include "fftjet/DiscreteGauss1d.hh"

namespace fftjet {
    DiscreteGauss1d::DiscreteGauss1d(
        const double sx, const unsigned nPhi)
        : AbsFrequencyKernel1d(),
          sPhi_(sx),
          nPhi_(nPhi)
    {
        // Check that the inputs are reasonable
        assert(sPhi_ > 0.0);
        assert(nPhi_);
    }

    std::complex<double> DiscreteGauss1d::operator()(
        const int iy, const double scale) const
    {
        const double cosvm1 = cos(iy*2.0*M_PI/nPhi_) - 1.0;
        const double dphisq = 4.0*M_PI*M_PI/nPhi_/nPhi_;
        const double sysq   = sPhi_*sPhi_*scale*scale/dphisq;
        const double fourier = exp(sysq*cosvm1);
        return std::complex<double>(fourier, 0.0);
    }
}
