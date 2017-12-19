//=========================================================================
// CUFFTFloatEngine.hh
//
// Single precision FFT engine for jet reconstruction which uses
// CUDA-based FFT library from NVIDIA. For more details, see
// http://en.wikipedia.org/wiki/CUDA
//
// I. Volobouev
// May 2008
//=========================================================================

#ifndef FFTJET_CUFFTFLOATENGINE_HH_
#define FFTJET_CUFFTFLOATENGINE_HH_

#include "fftjet/AbsFFTEngine.hh"
#include "cufft.h"

namespace fftjet {
    class CUFFTFloatEngine : public AbsFFTEngine<float, cufftComplex>
    {
    public:
        CUFFTFloatEngine(unsigned nEta, unsigned nPhi);
        virtual ~CUFFTFloatEngine();

        void transformForward(const float* from, cufftComplex* to) const;
        void transformBack(const cufftComplex* from, float* to) const;

        void multiplyTransforms(const cufftComplex* a, const cufftComplex* b,
                                cufftComplex* result) const;
        void addTransforms(const cufftComplex* a, const cufftComplex* b,
                           cufftComplex* result) const;
        void divideTransforms(const cufftComplex* a, const cufftComplex* b,
                              cufftComplex* result) const;
        void deconvolutionRatio(const cufftComplex* a, const cufftComplex* b,
                                cufftComplex* result) const;
        void subtractTransforms(const cufftComplex* a, const cufftComplex* b,
                                cufftComplex* result) const;
        void scaleTransform(const cufftComplex* a, float scale,
                            cufftComplex* result) const;
        void amplitudeSquared(const cufftComplex* a,
                              cufftComplex* result) const;
        double totalPower(const cufftComplex* a) const;
        void setToReal(cufftComplex* a, float value) const;

        void setTransformPoint(cufftComplex* points,
                               unsigned iEta, unsigned iPhi,
                               const std::complex<double>& value) const;
        std::complex<double> getTransformPoint(const cufftComplex* points,
                                               unsigned iEta,
                                               unsigned iPhi) const;

        cufftComplex* allocateComplex() const;
        void destroyComplex(cufftComplex*) const;        

    private:
        CUFFTFloatEngine(const CUFFTFloatEngine&);
        CUFFTFloatEngine& operator=(const CUFFTFloatEngine&);

        void checkFFTStatus(cufftResult status) const;

        float *in;
        cufftComplex *out;
        cufftHandle pf;
        cufftHandle pb;
    };
}

#endif // FFTJET_CUFFTFLOATENGINE_HH_
