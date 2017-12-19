//=========================================================================
// AbsFrequencyKernel1d.hh
//
// Interface class for 1d functions in the FFT transform space.
// Can be used to represent filters.
//
// This interface is much simpler than the interface used for real
// kernels. For the intended range of applications, we do not anticipate
// the need to represent delta functions in the frequency space and to
// generate random harmonics.
//
// I. Volobouev
// November 2009
//=========================================================================

#ifndef FFTJET_ABSFREQUENCYKERNEL1D_HH_
#define FFTJET_ABSFREQUENCYKERNEL1D_HH_

#include <complex>

namespace fftjet {
    struct AbsFrequencyKernel1d
    {
        virtual ~AbsFrequencyKernel1d() {}

        // i is the harmonic number (can be both positive and negative)
        virtual std::complex<double> operator()(int i, double scale) const = 0;
    };
}

#endif // FFTJET_ABSFREQUENCYKERNEL1D_HH_
