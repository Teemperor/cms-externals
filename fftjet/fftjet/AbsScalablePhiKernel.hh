//=========================================================================
// AbsScalablePhiKernel.hh
//
// Interface class for kernel functions which are products of delta
// function in eta and some reasonable distribution in phi. The shape
// of the distribution in phi is assumed to be independent of the scale.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_ABSSCALABLEPHIKERNEL_HH_
#define FFTJET_ABSSCALABLEPHIKERNEL_HH_

#include <cassert>

#include "fftjet/AbsPhiKernel.hh"

namespace fftjet {
    class AbsScalablePhiKernel : public AbsPhiKernel
    {
    public:
        inline AbsScalablePhiKernel(const double scaleFactor,
                                    const int scalePower)
            : phis_(scaleFactor), scalePower_(scalePower)
        {
            assert(phis_);
        }
        virtual ~AbsScalablePhiKernel() {}

    private:
        AbsScalablePhiKernel();

        inline double phiFcn(const double y, const double scale) const
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
            const double bwy(phis_*pscale);
            return unscaledPhiFcn(y/bwy)/bwy;
        }

        inline void phiSupport(const double scale,
                               double *phimin, double *phimax) const
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
            const double bwy(phis_*pscale);
            unscaledPhiSupport(phimin, phimax);
            *phimin *= bwy;
            *phimax *= bwy;
        }

        inline double phiRandom(const double rnd, const double scale) const
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
            const double bwy(phis_*pscale);
            return bwy*unscaledPhiRandom(rnd);
        }

        virtual double unscaledPhiFcn(double phi) const = 0;
        virtual void unscaledPhiSupport(double *phimin,
                                        double *phimax) const = 0;
        virtual double unscaledPhiRandom(double rnd) const = 0;

        const double phis_;
        const int scalePower_;
    };
}

#endif // FFTJET_ABSSCALABLEPHIKERNEL_HH_
