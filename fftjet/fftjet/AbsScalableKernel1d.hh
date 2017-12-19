//=========================================================================
// AbsScalableKernel1d.hh
//
// Interface class for scalable 1d kernel functions. Derive your kernel
// from this class if only the width of your kernel depends on scale but
// not the shape.
//
// I. Volobouev
// November 2009
//=========================================================================

#ifndef FFTJET_ABSSCALABLEKERNEL1D_HH_
#define FFTJET_ABSSCALABLEKERNEL1D_HH_

#include <cassert>
#include <cmath>

#include "fftjet/AbsKernel1d.hh"

namespace fftjet {
    class AbsScalableKernel1d : public AbsKernel1d
    {
    public:
        // xScaleFactor converts the scale into bandwidth value:
        // x_bandwidth = xScaleFactor*pow(scale, scalePower),
        inline AbsScalableKernel1d(const double xScaleFactor,
                                   const int scalePower)
            : xs_(xScaleFactor), scalePower_(scalePower)
        {
            assert(xs_);
        }
        virtual ~AbsScalableKernel1d() {}

        // Trivial accessors
        inline double xScaleFactor() const {return xs_;}
        inline int scalePower() const {return scalePower_;}

        // The function value
        inline double operator()(const double x, const double scale) const
        {
            double pscale;
            switch (scalePower_)
            {
            case -1:
                pscale = 1.0/scale;
                break;
            case 0:
                pscale = 1.0;
                break;
            case 1:
                pscale = scale;
                break;
            default:
                pscale = pow(scale, scalePower_);
                break;
            }
            const double bwx(xs_*pscale);
            return calculate(x/bwx)/bwx;
        }

        // The support interval
        inline KernelSupportInterval supportInterval(const double scale) const
        {
            double pscale;
            switch (scalePower_)
            {
            case -1:
                pscale = 1.0/scale;
                break;
            case 0:
                pscale = 1.0;
                break;
            case 1:
                pscale = scale;
                break;
            default:
                pscale = pow(scale, scalePower_);
                break;
            }
            return unscaledInterval() * (xs_*pscale);
        }

        // The randomizer
        inline double random(const double r1, const double scale) const
        {
            double pscale;
            switch (scalePower_)
            {
            case -1:
                pscale = 1.0/scale;
                break;
            case 0:
                pscale = 1.0;
                break;
            case 1:
                pscale = scale;
                break;
            default:
                pscale = pow(scale, scalePower_);
                break;
            }
            return xs_*pscale*unscaledRandom(r1);
        }

    private:
        AbsScalableKernel1d();

        virtual double calculate(double x) const = 0;
        virtual KernelSupportInterval unscaledInterval() const = 0;
        virtual double unscaledRandom(double r1) const = 0;

        const double xs_;
        const int scalePower_;
    };
}

#endif // FFTJET_ABSSCALABLEKERNEL1D_HH_
