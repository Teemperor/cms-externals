//=========================================================================
// FFTWEngine.hh
//
// FFT engine for jet reconstruction which uses the FFTW library
// from http://www.fftw.org. User code can not instantiate this
// class -- instead, use FFTWDoubleEngine and FFTWFloatEngine.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_FFTWENGINE_HH_
#define FFTJET_FFTWENGINE_HH_

#include "fftjet/AbsFFTEngine.hh"

namespace fftjet {
    template<typename Real, typename Complex, class Plan>
    class FFTWEngine : public AbsFFTEngine<Real, Complex>
    {
    public:
        virtual ~FFTWEngine() {}

        void transformForward(const Real* from, Complex* to) const;
        void transformBack(const Complex* from, Real* to) const;

        void multiplyTransforms(const Complex* a, const Complex* b,
                                Complex* result) const;
        void addTransforms(const Complex* a, const Complex* b,
                           Complex* result) const;
        void divideTransforms(const Complex* a, const Complex* b,
                              Complex* result) const;
        void deconvolutionRatio(const Complex* a, const Complex* b,
                                Complex* result) const;
        void subtractTransforms(const Complex* a, const Complex* b,
                                Complex* result) const;
        void scaleTransform(const Complex* a, Real scale,
                            Complex* result) const;
        void amplitudeSquared(const Complex* a,
                              Complex* result) const;
        double totalPower(const Complex* a) const;
        void setToReal(Complex* a, Real value) const;

        void setTransformPoint(Complex* points,
                               unsigned iEta, unsigned iPhi,
                               const std::complex<double>& value) const;
        std::complex<double> getTransformPoint(const Complex* points,
                                               unsigned iEta,
                                               unsigned iPhi) const;
        Complex* allocateComplex() const;
        void destroyComplex(Complex*) const;

    protected:
        FFTWEngine(unsigned nEta, unsigned nPhi);

        Real *in;
        Complex *out;
        Plan pf;
        Plan pb;

        void (*engineExecute)(const Plan);
        void* (*engineMalloc)(size_t);
        void (*engineFree)(void *);

    private:
        FFTWEngine(const FFTWEngine&);
        FFTWEngine& operator=(const FFTWEngine&);

        typedef AbsFFTEngine<Real, Complex> B;

        const unsigned nReal_;
        const unsigned nComplex_;
    };
}

#include "fftjet/FFTWEngine.icc"

#endif // FFTJET_FFTWENGINE_HH_
