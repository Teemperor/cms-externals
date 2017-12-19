//=========================================================================
// Kernels1d.hh
//
// Simple 1d kernel functions. More complicated kernels are implemented
// in their own separate files.
//
// I. Volobouev
// November 2009
//=========================================================================

#ifndef FFTJET_KERNELS1D_HH_
#define FFTJET_KERNELS1D_HH_

#include <vector>
#include "fftjet/AbsScalableKernel1d.hh"

namespace fftjet {
    // Gaussian kernel
    class Gauss1d : public AbsScalableKernel1d
    {
    public:
        inline Gauss1d(double sx, int scalePow)
            : AbsScalableKernel1d(sx, scalePow) {}
        inline Gauss1d(double sx, int scalePow, const std::vector<double>&)
            : AbsScalableKernel1d(sx, scalePow) {}

        inline bool isDensity() const {return true;}

        inline static int nParameters() {return 0;}

    private:
        double calculate(double x) const;
        KernelSupportInterval unscaledInterval() const;
        double unscaledRandom(double r1) const;
    };

    // Symmetric beta kernel
    class SymmetricBeta1d : public AbsScalableKernel1d
    {
    public:
        SymmetricBeta1d(double sx, int scalePow, double power);
        SymmetricBeta1d(double sx, int scalePow, const std::vector<double>& p);
        virtual ~SymmetricBeta1d() {}

        inline bool isDensity() const {return true;}
        inline double power() const {return n_;}

        inline static int nParameters() {return 1;}

    private:
        double calculate(double x) const;
        KernelSupportInterval unscaledInterval() const;
        double unscaledRandom(double r1) const;

        double calculateNorm() const;

        double n_;
        double norm_;
    };

    class DeltaFunction1d : public AbsKernel1d
    {
    public:
        inline explicit DeltaFunction1d(double area) : value_(area) {}
        DeltaFunction1d(double sx, int scalePow,
                        const std::vector<double>& params);
        virtual ~DeltaFunction1d() {}

        inline double area() const {return value_;}

        double operator()(double x, double scale) const;
        KernelSupportInterval supportInterval(double scale) const;
        double intervalAverage(double x, double sc, double dx) const;
        inline bool isDensity() const {return true;}
        double random(double rnd, double scale) const;

        inline static int nParameters() {return 1;}

    private:
        double value_;
    };
}

#endif // FFTJET_KERNELS1D_HH_
