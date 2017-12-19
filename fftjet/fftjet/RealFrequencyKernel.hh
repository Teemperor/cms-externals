//=========================================================================
// RealFrequencyKernel.hh
//
// 2d functions in the Fourier transform space whose transform
// is pure real. Implemented via real kernels.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_REALFREQUENCYKERNEL_HH_
#define FFTJET_REALFREQUENCYKERNEL_HH_

#include <cassert>

#include "fftjet/AbsFrequencyKernel.hh"
#include "fftjet/AbsKernel2d.hh"

namespace fftjet {
    class RealFrequencyKernel : public AbsFrequencyKernel
    {
    public:
        inline RealFrequencyKernel(const AbsKernel2d* realKernel,
                                   const bool assumePointerOwnership=false)
            : kern_(realKernel), owns_(assumePointerOwnership) {assert(kern_);}
        inline virtual ~RealFrequencyKernel()
        {
            if (owns_) delete const_cast<AbsKernel2d*>(kern_);
        }

        inline bool ownsPointer() const {return owns_;}

        inline std::complex<double> operator()(
            const int ix, const int iy, const double scale) const
        {
            // It is natural to invert the scale in the frequency space
            return std::complex<double>((*kern_)(ix, iy, 1.0/scale), 0.0);
        }

    private:
        RealFrequencyKernel();
        RealFrequencyKernel(const RealFrequencyKernel&);
        RealFrequencyKernel& operator=(const RealFrequencyKernel&);

        const AbsKernel2d* const kern_;
        const bool owns_;
    };
}

#endif // FFTJET_REALFREQUENCYKERNEL_HH_
