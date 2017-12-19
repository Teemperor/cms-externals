//=========================================================================
// AbsScalableKernel.hh
//
// Interface class for scalable 2d kernel functions. Derive your kernel
// from this class if only the width of your kernel depends on scale but
// not the shape.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_ABSSCALABLEKERNEL_HH_
#define FFTJET_ABSSCALABLEKERNEL_HH_

#include <cassert>
#include <cmath>

#include "fftjet/AbsKernel2d.hh"

namespace fftjet {
    class AbsScalableKernel : public AbsKernel2d
    {
    public:
        // xScaleFactor and yScaleFactor convert the scale
        // into x and y bandwidth values, respectively:
        // x_bandwidth = xScaleFactor*pow(scale, scalePower),
        // y_bandwidth = yScaleFactor*pow(scale, scalePower)
        inline AbsScalableKernel(const double xScaleFactor,
                                 const double yScaleFactor,
                                 const int scalePower)
            : xs_(xScaleFactor),
              ys_(yScaleFactor),
              geomMean_(sqrt(xScaleFactor*yScaleFactor)),
              scalePower_(scalePower)
        {
            assert(xs_ > 0.0);
            assert(ys_ > 0.0);
        }
        virtual ~AbsScalableKernel() {}

        // The following member sets eta to phi (or x to y) scale ratio
        inline void setScaleRatio(const double r)
        {
            assert(r > 0.0);
            const double sqr = sqrt(r);
            xs_ = geomMean_*sqr;
            ys_ = geomMean_/sqr;
        }

        // Trivial accessors
        inline double xScaleFactor() const {return xs_;}
        inline double yScaleFactor() const {return ys_;}
        inline int scalePower() const {return scalePower_;}

        // The function value
        inline double operator()(const double x, const double y,
                                 const double scale) const
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
            const double bwy(ys_*pscale);
            return calculate(x/bwx, y/bwy)/bwx/bwy;
        }

        // The support rectangle
        inline void supportRectangle(const double scale,
                                     KernelSupportRectangle *r) const
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
            const double bwy(ys_*pscale);
            unscaledRectangle(r);
            r->xmin *= bwx;
            r->xmax *= bwx;
            r->ymin *= bwy;
            r->ymax *= bwy;
        }

        // The randomizer
        inline void random(const double r1, const double r2,
                           const double scale,
                           double* px, double* py) const
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
            unscaledRandom(r1, r2, px, py);
            *px *= (xs_*pscale);
            *py *= (ys_*pscale);
        }

    private:
        AbsScalableKernel();

        virtual double calculate(double x, double y) const = 0;
        virtual void unscaledRectangle(KernelSupportRectangle *r) const = 0;
        virtual void unscaledRandom(double r1, double r2,
                                    double* px, double* py) const = 0;
        double xs_;
        double ys_;
        const double geomMean_;
        const int scalePower_;
    };
}

#endif // FFTJET_ABSSCALABLEKERNEL_HH_
