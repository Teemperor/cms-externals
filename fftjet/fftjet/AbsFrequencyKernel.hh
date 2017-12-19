//=========================================================================
// AbsFrequencyKernel.hh
//
// Interface class for 2d functions in the FFT transform space.
// Can be used to represent filters.
//
// This interface is much simpler than the interface used for real
// kernels. For the intended range of applications, we do not anticipate
// the need to represent delta functions in the frequency space and to
// generate random harmonics.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_ABSFREQUENCYKERNEL_HH_
#define FFTJET_ABSFREQUENCYKERNEL_HH_

#include <complex>

namespace fftjet {
    struct AbsFrequencyKernel
    {
        virtual ~AbsFrequencyKernel() {}

        // ix and iy are harmonics. The convention used here
        // is that they can be both positive and negative.
        virtual std::complex<double> operator()(int ix, int iy,
                                                double scale) const = 0;
    };
}

#endif // FFTJET_ABSFREQUENCYKERNEL_HH_
