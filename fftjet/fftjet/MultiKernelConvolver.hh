//=========================================================================
// MultiKernelConvolver.hh
//
// Class for performing kernel convolutions by FFT using a set of kernels
//
// I. Volobouev
// April 2008
//=========================================================================

#ifndef FFTJET_MULTIKERNELCONVOLVER_HH_
#define FFTJET_MULTIKERNELCONVOLVER_HH_

#include <vector>

#include "fftjet/AbsKernelConvolver.hh"
#include "fftjet/AbsKernel2d.hh"
#include "fftjet/AbsFrequencyKernel.hh"

namespace fftjet {
    class KernelSet
    {
    public:
        explicit KernelSet(bool ownsPointers,
                           double regularizationFraction=0.0);
        ~KernelSet();

        inline bool ownsPointers() const {return ownsPointers_;}
        inline double regularizationFraction() const {return regFraction_;}
        bool isEmpty() const;

        // The regularization fraction (fraction of frequencies to zero out)
        // should be between 0 and 1.
        void setRegularizationFraction(double fraction);

        std::vector<const AbsFrequencyKernel*> filter;
        std::vector<const AbsKernel2d*> numerator;
        std::vector<const AbsKernel2d*> denominator;
        std::vector<const AbsKernel2d*> denoiser;

    private:
        KernelSet();
        KernelSet(const KernelSet&);
        KernelSet& operator=(const KernelSet&);

        double regFraction_;
        const bool ownsPointers_;
    };

    template<typename Real, typename Complex>
    class MultiKernelConvolver : public AbsKernelConvolver<Real, Complex>
    {
    public:
        // This MultiKernelConvolver will not own the AbsFFTEngine
        // and KernelSet objects. It is the responsibility of
        // the user of this class to ensure that these objects
        // are not destroyed before the MultiKernelConvolver itself.
        MultiKernelConvolver(const AbsFFTEngine<Real,Complex>* fftEngine,
                             const KernelSet* kernelSequence,
                             unsigned minFixBin=0, unsigned maxFixBin=0);
        virtual ~MultiKernelConvolver();

    protected:
        virtual KernelData<Real,Complex> buildKernelImage(double scale);

        const KernelSet* kernelSequence;

    private:
        void regularizeDenominator(Complex* imageBuf) const;

        Complex* denomBuf;
    };
}

#include "fftjet/MultiKernelConvolver.icc"

#endif // FFTJET_MULTIKERNELCONVOLVER_HH_
