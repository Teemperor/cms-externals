//=========================================================================
// AbsFFTEngine.hh
//
// Abstract FFT engine for jet reconstruction. User code should normally
// avoid dependence on concrete implementations of the engine and use
// this header instead.
//
// I. Volobouev
// January 2008
//=========================================================================

#ifndef FFTJET_ABSFFTENGINE_HH_
#define FFTJET_ABSFFTENGINE_HH_

#include <complex>
#include <cassert>

namespace fftjet {
    // Normally, each derived engine should have a constructor like this:
    //
    // ClassName(unsigned nEta, unsigned nPhi);
    //
    // so that it can call the base constructor with appropriate arguments.
    // If the number of eta points specified is 1, the engine should perform
    // 1-d transforms with "nPhi" points.
    //
    // Normally, "Real" type should be one of "float", "double", or
    // "long double". The "Complex" type should represent complex numbers
    // at the same level of precision.
    //
    template<typename Real, typename Complex>
    class AbsFFTEngine
    {
    public:
        inline AbsFFTEngine(const unsigned nEta, const unsigned nPhi)
            : nEta_(nEta), nPhi_(nPhi)
        {
            assert(nEta_);
            assert(nPhi_);
        }
        virtual ~AbsFFTEngine() {}

        // It is assumed that all complex buffers on which member
        // functions operate are created using "allocateComplex()"
        // method.
        virtual void transformForward(const Real* from, Complex* to) const=0;
        virtual void transformBack(const Complex* from, Real* to) const=0;

        // Arithmetic operations with complex 2d arrays. Note that
        // implementations of all of these operations should allow for
        // "result" to point into the same location as either "a" or "b".
        virtual void multiplyTransforms(const Complex* a, const Complex* b,
                                        Complex* result) const=0;
        virtual void addTransforms(const Complex* a, const Complex* b,
                                   Complex* result) const=0;
        // The following function calculates a/b
        virtual void divideTransforms(const Complex* a, const Complex* b,
                                      Complex* result) const=0;
        // The following function calculates a/b, but it does not
        // complain in case both a and b are complex zeros. In this
        // case the result is set to complex zero.
        virtual void deconvolutionRatio(const Complex* a, const Complex* b,
                                        Complex* result) const=0;
        // The following function calculates a - b
        virtual void subtractTransforms(const Complex* a, const Complex* b,
                                        Complex* result) const=0;
        // The following function multiplies a transform by a real number.
        // Note that "a" and "result" should be allowed to point at the same
        // location.
        virtual void scaleTransform(const Complex* a, Real scale,
                                    Complex* result) const=0;
        // The following function should allow "a" and "result"
        // to point at the same location
        virtual void amplitudeSquared(const Complex* a,
                                      Complex* result) const=0;
        // The following function computes the sum of absoule
        // values squared
        virtual double totalPower(const Complex* a) const=0;

        // The following function sets all array elements to Re=value, Im=0
        virtual void setToReal(Complex* a, Real value) const=0;

        // The following function sets all array elements to complex zero.
        // Override if a faster implementation is available.
        inline virtual void zeroOut(Complex* a) const
        {
            setToReal(a, static_cast<Real>(0));
        }

        // The following functions can be used to get/set transform points
        // using the standard C++ complex class which, in general, will
        // not be the class used by the engine to represent complex numbers
        virtual void setTransformPoint(Complex* points,
                                       unsigned iEta, unsigned iPhi,
                                       const std::complex<double>& p) const=0;
        virtual std::complex<double> getTransformPoint(const Complex* points,
                                                       unsigned iEta,
                                                       unsigned iPhi) const=0;

        // Allocation functions should generate 2d arrays of proper
        // size and alignment for the particular engine used. User code
        // should manage creation and destruction of such arrays by
        // calling these engine functions. It should be safe to call
        // "destroyComplex" on a NULL pointer.
        virtual Complex* allocateComplex() const=0;
        virtual void destroyComplex(Complex*) const=0;

        // Dimensionality of the real arrays. Note that the dimensionality
        // of the complex arrays is likely to be different, and it will
        // depend on particular engine implementation.
        inline unsigned nEta() const {return nEta_;}
        inline unsigned nPhi() const {return nPhi_;}

    protected:
        const unsigned nEta_;
        const unsigned nPhi_;

    private:
        // Disable the default constructor so that the derived classes
        // are forced to set "nEta_" and "nPhi_" correctly
        AbsFFTEngine();
    };
}

#endif // FFTJET_ABSFFTENGINE_HH_
