#include <cassert>
#include <cmath>

#include "fftjet/AbsSymmetricKernel.hh"

namespace fftjet {
    double AbsSymmetricKernel::axialWeight(const double rmax,
                                           const unsigned nsteps) const
    {
        assert(nsteps);
        assert(rmax >= 0.0);
        if (rmax == 0.0)
            return 0.0;
        const double step = rmax/nsteps;
        const double quarter = step/4.0;
        long double sum = 0.0L;
        double r = 0.0;
        // Bode's integration rule is used in which we do not
        // need to add the very first point (it is 0).
        for (unsigned i=0; i<nsteps; ++i)
        {
            r = i*step;
            r += quarter;
            sum += 64*r*eval(r*r);
            r += quarter;
            sum += 24*r*eval(r*r);
            r += quarter;
            sum += 64*r*eval(r*r);
            r += quarter;
            sum += 28*r*eval(r*r);
        }
        sum -= 14*r*eval(r*r);
        return M_PI/90.0*step*sum;
    }

    double AbsSymmetricKernel::circularResponse(
        const double x, const double y,
        const unsigned n, const double a) const
    {
        long double sum = 0.0L;
        for (unsigned i=0; i<n; ++i)
        {
            const double phi = i*2.0*M_PI/n;
            sum += calculate(x - a*cos(phi), y - a*sin(phi));
        }
        return sum;
    }

    double AbsSymmetricKernel::scanForResolution(
        const unsigned nkernels, const double r0,
        const double scaleFactor, const double releps) const
    {
        const double sFMinusOne = scaleFactor - 1.0;
        assert(sFMinusOne > 0.0);
        assert(nkernels > 1);
        assert(releps > 0.0);

        const double peakDetectionEps = 1.e-10;
        const double hrad = halfWeightRadius();
        const double scanStart = r0 > 0.0 ? r0 : hrad*sFMinusOne;
        const double scanEnd = hrad/sFMinusOne;
        const double eps = hrad*sFMinusOne*0.001;

        // Check if there is a minimum at the center of the kernel
        {
            const double v0 = eval(0.0);
            assert(v0 >= 0.0);
            if (v0 == 0.0)
                return 0.0;
            int dir = 0;
            for (double epsfac=1.0; epsfac*eps<hrad; epsfac*=2.0)
            {
                const double v1 = eval(epsfac*eps*epsfac*eps);
                if (fabs((v1 - v0)/v0) > peakDetectionEps)
                {
                    if (v1 > v0)
                        dir = 1;
                    else
                        dir = -1;
                    break;
                }
            }
            if (dir > 0)
                return 0.0;
        }

        // Some useful variables
        const double phi = M_PI/nkernels;
        const double cphi = cos(phi);
        const double sphi = sin(phi);

        // Scan the scale
        double rmin = 0.0, rmax = 0.0;
        for (double a=scanStart; a<scanEnd; a*=scaleFactor)
        {
            const double v0 = circularResponse(0.0, 0.0, nkernels, a);
            assert(v0 >= 0.0);
            if (v0 == 0.0)
            {
                rmax = a;
                break;
            }
            int dir = 0;
            for (double epsfac=1.0; epsfac*eps<a; epsfac*=2.0)
            {
                const double dr = epsfac*eps;
                const double v1 = circularResponse(dr, 0.0, nkernels, a);
                const double v2 = circularResponse(dr*cphi, dr*sphi, nkernels, a);
                const double diff1 = fabs((v1 - v0)/v0);
                const double diff2 = fabs((v2 - v0)/v0);

                // Do we have a maximum?
                if (diff1 > peakDetectionEps && v1 < v0 &&
                    diff2 > peakDetectionEps && v2 < v0)
                {
                    dir = -1;
                    break;
                }

                // Do we have a minimum in at least one direction?
                if ((diff1 > peakDetectionEps && v1 > v0) ||
                    (diff2 > peakDetectionEps && v2 > v0))
                {
                    dir = 1;
                    break;
                }
            }
            if (dir > 0)
            {
                // Not a maximum
                rmax = a;
                break;
            }
            rmin = a;
        }

        if (rmin < rmax)
        {
            // We have managed to bracket the resolution distance
            while (2.0*(rmax - rmin)/(rmax + rmin) > releps)
            {
                const double a = (rmax + rmin)/2.0;
                const double v0 = circularResponse(0.0, 0.0, nkernels, a);
                assert(v0 >= 0.0);
                if (v0 == 0.0)
                {
                    rmax = a;
                    continue;
                }
                int dir = 0;
                for (double epsfac=1.0; epsfac*eps<a; epsfac*=2.0)
                {
                    const double dr = epsfac*eps;
                    const double v1 = circularResponse(dr, 0.0, nkernels, a);
                    const double v2 = circularResponse(dr*cphi, dr*sphi, nkernels, a);
                    const double diff1 = fabs((v1 - v0)/v0);
                    const double diff2 = fabs((v2 - v0)/v0);

                    // Do we have a maximum?
                    if (diff1 > peakDetectionEps && v1 < v0 &&
                        diff2 > peakDetectionEps && v2 < v0)
                    {
                        dir = -1;
                        break;
                    }

                    // Do we have a minimum in at least one direction?
                    if ((diff1 > peakDetectionEps && v1 > v0) ||
                        (diff2 > peakDetectionEps && v2 > v0))
                    {
                        dir = 1;
                        break;
                    }
                }
                if (dir > 0)
                {
                    // Not a maximum
                    rmax = a;
                    continue;
                }
                rmin = a;
            }
            return (rmax + rmin)/2.0;
        }

        assert(!"Default resolution calculation function failed");
        return -1.0;
    }
}
