//=========================================================================
// FrequencySequentialConvolver.hh
//
// Class for performing sequential convolutions by FFT using 1d kernels
// whose representations are given in the frequency domain
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_FREQUENCYSEQUENTIALCONVOLVER_HH_
#define FFTJET_FREQUENCYSEQUENTIALCONVOLVER_HH_

#include <climits>

#include "fftjet/AbsSequentialConvolver.hh"
#include "fftjet/AbsFrequencyKernel1d.hh"

namespace fftjet {
    template<typename Real, typename Complex>
    class FrequencySequentialConvolver :
        public AbsSequentialConvolver<Real, Complex>
    {
    public:
        FrequencySequentialConvolver(
            const AbsFFTEngine<Real,Complex>* etaEngine,
            const AbsFFTEngine<Real,Complex>* phiEngine,
            const AbsFrequencyKernel1d* etaKernel,
            const AbsFrequencyKernel1d* phiKernel,
            const std::vector<double>& phiScales,
            unsigned minEtaBin=0, unsigned maxEtaBin=UINT_MAX,
            bool fixEfficiency=false);
        inline virtual ~FrequencySequentialConvolver() {}

    protected:
        virtual KernelData<Real,Complex> buildKernelImageEta(
            double scale, const AbsFFTEngine<Real,Complex>* engine,
            unsigned minFixBin, unsigned maxFixBin);
        virtual KernelData<Real,Complex> buildKernelImagePhi(
            double scale, const AbsFFTEngine<Real,Complex>* engine);

    private:
        KernelData<Real,Complex> buildImage(
            double scale, const AbsFFTEngine<Real,Complex>* engine,
            const AbsFrequencyKernel1d* kernel,
            unsigned minFixBin, unsigned maxFixBin);

        const AbsFrequencyKernel1d* const etaKernel;
        const AbsFrequencyKernel1d* const phiKernel;
    };
}

#include "fftjet/FrequencySequentialConvolver.icc"

#endif // FFTJET_FREQUENCYSEQUENTIALCONVOLVER_HH_
