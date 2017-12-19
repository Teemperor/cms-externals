//=========================================================================
// KernelConvolver.hh
//
// Class for performing kernel convolutions by FFT using a single kernel
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_KERNELCONVOLVER_HH_
#define FFTJET_KERNELCONVOLVER_HH_

#include "fftjet/AbsKernelConvolver.hh"
#include "fftjet/AbsKernel2d.hh"

namespace fftjet {
    template<typename Real, typename Complex>
    class KernelConvolver : public AbsKernelConvolver<Real, Complex>
    {
    public:
        // This KernelConvolver will not own the AbsFFTEngine
        // and AbsKernel2d objects. It is the responsibility of
        // the user of this class to ensure that these objects
        // are not destroyed before the KernelConvolver itself.
        //
        // If the "minFixBin" and "maxFixBin" arguments are
        // changed from their default values, it means that
        // the convolver will try to fix the bump reconstruction
        // efficiency (will apply weights to the data to simulate
        // the effects of boundary kernel)
        KernelConvolver(const AbsFFTEngine<Real,Complex>* fftEngine,
                        const AbsKernel2d* kernel,
                        unsigned minFixBin=0, unsigned maxFixBin=0);
        inline virtual ~KernelConvolver() {}

    protected:
        virtual KernelData<Real,Complex> buildKernelImage(double scale);

        const AbsKernel2d* const kernel;
    };
}

#include "fftjet/KernelConvolver.icc"

#endif // FFTJET_KERNELCONVOLVER_HH_
