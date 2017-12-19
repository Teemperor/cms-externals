//=========================================================================
// SequentialConvolver.hh
//
// Class for performing kernel convolutions by FFT using separate kernels
// in eta and phi
//
// I. Volobouev
// November 2009
//=========================================================================

#ifndef FFTJET_SEQUENTIALCONVOLVER_HH_
#define FFTJET_SEQUENTIALCONVOLVER_HH_

#include <climits>

#include "fftjet/AbsSequentialConvolver.hh"
#include "fftjet/AbsKernel1d.hh"

namespace fftjet {
    template<typename Real, typename Complex>
    class SequentialConvolver : public AbsSequentialConvolver<Real, Complex>
    {
    public:
        // This SequentialConvolver will not own the AbsFFTEngine
        // and AbsKernel1d objects. It is the responsibility of
        // the user of this class to ensure that these objects
        // are not destroyed before the SequentialConvolver itself.
        //
        // The convolver will try to fix the bump reconstruction
        // efficiency (will apply weights to the data to simulate
        // the effects of boundary kernel) in case "fixEfficiency"
        // parameter is "true".
        //
        SequentialConvolver(const AbsFFTEngine<Real,Complex>* etaEngine,
                            const AbsFFTEngine<Real,Complex>* phiEngine,
                            const AbsKernel1d* etaKernel,
                            const AbsKernel1d* phiKernel,
                            const std::vector<double>& phiScales,
                            unsigned minEtaBin=0, unsigned maxEtaBin=UINT_MAX,
                            bool fixEfficiency=false);
        inline virtual ~SequentialConvolver() {}

    protected:
        virtual KernelData<Real,Complex> buildKernelImageEta(
            double scale, const AbsFFTEngine<Real,Complex>* engine,
            unsigned minFixBin, unsigned maxFixBin);
        virtual KernelData<Real,Complex> buildKernelImagePhi(
            double scale, const AbsFFTEngine<Real,Complex>* engine);

    private:
        KernelData<Real,Complex> buildImage(
            double scale, const AbsFFTEngine<Real,Complex>* engine,
            const AbsKernel1d* kernel,
            unsigned minFixBin, unsigned maxFixBin);

        const AbsKernel1d* const etaKernel;
        const AbsKernel1d* const phiKernel;
    };
}

#include "fftjet/SequentialConvolver.icc"

#endif // FFTJET_SEQUENTIALCONVOLVER_HH_
