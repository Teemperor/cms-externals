//=========================================================================
// RealFrequencyKernel1d.hh
//
// 1d functions in the Fourier transform space whose transform
// is pure real. Implemented via real kernels.
//
// I. Volobouev
// November 2009
//=========================================================================

#ifndef FFTJET_REALFREQUENCYKERNEL1D_HH_
#define FFTJET_REALFREQUENCYKERNEL1D_HH_

#include <cassert>

#include "fftjet/AbsFrequencyKernel1d.hh"
#include "fftjet/AbsKernel1d.hh"

namespace fftjet {
    class RealFrequencyKernel1d : public AbsFrequencyKernel1d
    {
    public:
        inline RealFrequencyKernel1d(const AbsKernel1d* realKernel,
                                     const bool assumePointerOwnership=false)
            : kern_(realKernel), owns_(assumePointerOwnership) {assert(kern_);}
        inline virtual ~RealFrequencyKernel1d()
        {
            if (owns_) delete const_cast<AbsKernel1d*>(kern_);
        }

        inline bool ownsPointer() const {return owns_;}

        inline std::complex<double> operator()(
            const int ix, const double scale) const
        {
            // It is natural to invert the scale in the frequency space
            return std::complex<double>((*kern_)(ix, 1.0/scale), 0.0);
        }

    private:
        RealFrequencyKernel1d();
        RealFrequencyKernel1d(const RealFrequencyKernel1d&);
        RealFrequencyKernel1d& operator=(const RealFrequencyKernel1d&);

        const AbsKernel1d* const kern_;
        const bool owns_;
    };
}

#endif // FFTJET_REALFREQUENCYKERNEL1D_HH_
