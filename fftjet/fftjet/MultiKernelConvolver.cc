#include <cassert>

#include "fftjet/MultiKernelConvolver.hh"

namespace fftjet {
    KernelSet::~KernelSet()
    {
        if (ownsPointers_)
        {
            unsigned n = numerator.size();
            for (unsigned i=0; i<n; ++i)
                delete numerator[i];

            n = denominator.size();
            for (unsigned i=0; i<n; ++i)
                delete denominator[i];

            n = filter.size();
            for (unsigned i=0; i<n; ++i)
                delete filter[i];

            n = denoiser.size();
            for (unsigned i=0; i<n; ++i)
                delete denoiser[i];
        }
    }

    KernelSet::KernelSet(const bool ownsPointers,
                         const double regularizationFraction)
        : ownsPointers_(ownsPointers)
    {
        setRegularizationFraction(regularizationFraction);
    }

    bool KernelSet::isEmpty() const
    {
        return numerator.size() == 0 &&
               denominator.size() == 0 &&
               filter.size() == 0 &&
               denoiser.size() == 0;
    }

    void KernelSet::setRegularizationFraction(const double fraction)
    {
        assert(fraction >= 0.0);
        assert(fraction <= 1.0);
        regFraction_ = fraction;
    }
}
