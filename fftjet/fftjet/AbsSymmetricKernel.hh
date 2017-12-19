//=========================================================================
// AbsSymmetricKernel.hh
//
// Interface class for scalable 2d kernel functions with rotational
// symmetry
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_ABSSYMMETRICKERNEL_HH_
#define FFTJET_ABSSYMMETRICKERNEL_HH_

#include <cmath>
#include "fftjet/AbsScalableKernel.hh"

namespace fftjet {
    class AbsSymmetricKernel : public AbsScalableKernel
    {
    public:
        inline AbsSymmetricKernel(const double xScaleFactor,
                                  const double yScaleFactor,
                                  const int scalePower)
            : AbsScalableKernel(xScaleFactor, yScaleFactor, scalePower),
              halfW_(-1.0), reso_(-1.0) {}
        virtual ~AbsSymmetricKernel() {}

        // The following function returns the radius of the circle
        // enclosing one half of the unscaled kernel mass.
        // This quantity is useful for selecting compatible scales
        // for different kernels. This function is just a memoized
        // version of "calculateHalfWeightRadius()".
        inline double halfWeightRadius() const
        {
            if (halfW_ < 0.0)
                halfW_ = calculateHalfWeightRadius();
            return halfW_;
        }

        // The following function returns the unscaled kernel mass
        // within the circle of radius r. The "nsteps" parameter can be
        // used to specify the number of numerical integration steps.
        // Override if a better than default implementation is available.
        virtual double axialWeight(double r, unsigned nsteps=1000) const;

        // The following function returns the minimum value of r0
        // for which convolution of the unscaled kernel with the function
        // delta(r - r0)/(2 Pi r) results in a local minimum at the center
        // of coordinates (0, 0). Here, "delta" is the Dirac delta function.
        // This function is just a memoized version of "calculateResolution()".
        inline double circleResolution() const
        {
            if (reso_ < 0.0)
                reso_ = calculateResolution();
            return reso_;
        }

    protected:
        // The following function is used to implement the default
        // behavior of the "calculateResolution()" function. The scan
        // of the radial distance is performed from the lowest value
        // "rmin" to some reasonably selected large value by multiplying
        // the distance with "scaleFactor" on every cycle. "scaleFactor"
        // argument must be slightly larger than 1.0. If "rmin" argument
        // is 0.0 or negative, a reasonable starting radius will be
        // selected automatically. "nkernels" kernels are placed
        // around the circle at this increasing radial distance.
        // The cycling stops (it is switched to root finding) when
        // a minimum is found at the circle center.
        double scanForResolution(unsigned nkernels, double rmin=0.0,
                                 double scaleFactor=1.01,
                                 double eps=1.e-6) const;
    private:
        mutable double halfW_;
        mutable double reso_;

        inline double calculate(const double x, const double y) const
        {
            return eval(x*x + y*y);
        }
        inline void unscaledRectangle(KernelSupportRectangle *r) const
        {
            const double d(supportDistance());
            r->xmin = -d;
            r->xmax = d;
            r->ymin = -d;
            r->ymax = d;
        }
        inline void unscaledRandom(const double r1, const double r2,
                                   double* px, double* py) const
        {
            const double r(randomRadius(r1));
            const double phi(2.0*M_PI*r2);
            *px = r*cos(phi);
            *py = r*sin(phi);
        }

        // The following functions must be implemented
        // by the derived classes
        virtual double eval(double rsquared) const = 0;
        virtual double supportDistance() const = 0;
        virtual double randomRadius(double rnd) const = 0;

        // Override the following functions if better than default
        // implementations are available
        inline virtual double calculateHalfWeightRadius() const
        {
            return randomRadius(0.5);
        }
        inline virtual double calculateResolution() const
        {
            return scanForResolution(6);
        }

        // Helper function for calculating the sum value of n kernels
        // placed on a circle with radius "a"
        double circularResponse(double x, double y,
                                unsigned n, double a) const;
    };
}

#endif // FFTJET_ABSSYMMETRICKERNEL_HH_
