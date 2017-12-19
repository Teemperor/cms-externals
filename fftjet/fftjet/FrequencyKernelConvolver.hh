//=========================================================================
// FrequencyKernelConvolver.hh
//
// Class for performing kernel convolutions by FFT using a single kernel
// whose representation is given in the frequency domain
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_FREQUENCYKERNELCONVOLVER_HH_
#define FFTJET_FREQUENCYKERNELCONVOLVER_HH_

#include "fftjet/AbsKernelConvolver.hh"
#include "fftjet/AbsFrequencyKernel.hh"

namespace fftjet {
    template<typename Real, typename Complex>
    class FrequencyKernelConvolver : public AbsKernelConvolver<Real, Complex>
    {
    public:
        // This FrequencyKernelConvolver will not own the AbsFFTEngine
        // and AbsFrequencyKernel objects. It is the responsibility
        // of the user of this class to ensure that these objects
        // are not destroyed before the KernelConvolver itself.
        //
        // If the "minFixBin" and "maxFixBin" arguments are
        // changed from their default values, it means that
        // the convolver will try to fix the bump reconstruction
        // efficiency (will apply weights to the data to simulate
        // the effects of boundary kernel).
        FrequencyKernelConvolver(const AbsFFTEngine<Real,Complex>* fftEngine,
                                 const AbsFrequencyKernel* kernel,
                                 unsigned minFixBin=0, unsigned maxFixBin=0);
        inline virtual ~FrequencyKernelConvolver() {}

    protected:
        virtual KernelData<Real,Complex> buildKernelImage(double scale);

        const AbsFrequencyKernel* const kernel;
    };
}

#include "fftjet/FrequencyKernelConvolver.icc"

#endif // FFTJET_FREQUENCYKERNELCONVOLVER_HH_
