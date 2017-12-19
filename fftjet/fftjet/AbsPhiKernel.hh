//=========================================================================
// AbsPhiKernel.hh
//
// Interface class for kernel functions which are products of delta
// function in eta and some reasonable distribution in phi
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_ABSPHIKERNEL_HH_
#define FFTJET_ABSPHIKERNEL_HH_

#include <cfloat>

#include "fftjet/AbsKernel2d.hh"

namespace fftjet {
    class AbsPhiKernel : public AbsKernel2d
    {
    public:
        virtual ~AbsPhiKernel() {}

        // Attempt to set the scale ratio is ignored
        inline void setScaleRatio(double) {}

        // The function value
        inline double operator()(const double x, double, double) const
        {
            if (x == 0.0)
                assert(!"Attempt to evaluate delta function at 0");
            return 0.0;
        }

        // The support rectangle
        inline void supportRectangle(const double scale,
                                     KernelSupportRectangle *r) const
        {
            phiSupport(scale, &r->ymin, &r->ymax);
            const double eps = 1.0/sqrt(DBL_MAX);
            r->xmin = -eps;
            r->xmax = eps;
        }

        virtual inline double rectangleAverage(
            const double x, const double y, const double scale,
            const double dx, const double dy) const
        {
            const double halfx(fabs(dx/2.0));
            if (x - halfx <= 0.0 && x + halfx > 0.0)
            {
                // 0.28867513459481288225 is 1/(2 sqrt(3))
                const double yshift = 0.28867513459481288225*dy;
                return (phiFcn(y - yshift, scale) + 
                        phiFcn(y + yshift, scale))/2.0/dx;
            }
            else
                return 0.0;
        }

        // The randomizer
        inline void random(double, const double r2,
                           const double scale,
                           double* px, double* py) const
        {
            *px = 0.0;
            *py = phiRandom(r2, scale);
        }

    private:
        virtual double phiFcn(double phi, double scale) const = 0;
        virtual void phiSupport(double scale,
                                double *phimin, double *phimax) const = 0;
        virtual double phiRandom(double rnd, double scale) const = 0;
    };
}

#endif // FFTJET_ABSPHIKERNEL_HH_
