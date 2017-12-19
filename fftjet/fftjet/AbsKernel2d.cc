#include <cmath>
#include <cassert>

#include "fftjet/AbsKernel2d.hh"

inline static double deltaPhi(const double phi1, const double phi2)
{
    double delta = phi1 - phi2;
    while (delta > M_PI)
        delta -= 2.0*M_PI;
    while (delta < -M_PI)
        delta += 2.0*M_PI;
    return delta;
}

namespace fftjet {
    double AbsKernel2d::rectangleAverage(const double x, const double y,
                                         const double scale,
                                         const double dx, const double dy) const
    {
        // Abramowitz, Stegun 4-point formula from section 25.4.62
        const double ave_delta = 0.28867513459481288225; // this is 1/(2 sqrt(3))
        const double xshift = ave_delta*dx;
        const double yshift = ave_delta*dy;
        return (operator()(x - xshift, y - yshift, scale) + 
                operator()(x - xshift, y + yshift, scale) +
                operator()(x + xshift, y - yshift, scale) +
                operator()(x + xshift, y + yshift, scale))/4.0;
    }

    void AbsKernel2d::optimalLineScan(
        const unsigned nx, const double xmin, const double xmax,
        const double rmin, const double rmax,
        unsigned* nxmin, unsigned* nxmax) const
    {
        const double step = (xmax - xmin)/nx;
        if (rmin > xmin)
            *nxmin = static_cast<unsigned>((rmin - xmin)/step);
        else
            *nxmin = 0;
        if (rmax < xmax)
        {
            if (rmax >= xmin)
            {
                *nxmax = static_cast<unsigned>((rmax - xmin)/step);
                if (*nxmax < nx)
                    ++*nxmax;
            }
            else
                *nxmax = 0;
        }
        else
            *nxmax = nx;
    }

    bool AbsKernel2d::optimalPlaneScan(
        const unsigned nx, const unsigned ny,
        const double x0, const double y0, const double scale,
        const double xminscan, const double xmaxscan,
        const double yminscan, const double ymaxscan,
        unsigned* nxmin, unsigned* nxmax,
        unsigned* nymin, unsigned* nymax) const
    {
        KernelSupportRectangle r;
        supportRectangle(scale, &r);
        const double xmin = x0 + r.xmin;
        const double xmax = x0 + r.xmax;
        const double ymin = y0 + r.ymin;
        const double ymax = y0 + r.ymax;

        optimalLineScan(nx, xminscan, xmaxscan, xmin, xmax, nxmin, nxmax);
        optimalLineScan(ny, yminscan, ymaxscan, ymin, ymax, nymin, nymax);

        return *nxmin > 0 || *nxmax < nx || 
               *nymin > 0 || *nymax < ny;
    }

    bool AbsKernel2d::optimalCylinderScan(
        const unsigned nEta, const unsigned nPhi,
        const double eta0, const double /* phi0 */, const double scale,
        const double etaMin, const double etaMax,
        const double /* phiBin0Edge */,
        unsigned* netaMin, unsigned* netaMax,
        unsigned* nphiMin, unsigned* nphiMax) const
    {
        KernelSupportRectangle r;
        supportRectangle(scale, &r);
        const double xmin = eta0 + r.xmin;
        const double xmax = eta0 + r.xmax;

        optimalLineScan(nEta, etaMin, etaMax, xmin, xmax, netaMin, netaMax);

        // Figuring out the optimal scan in phi is not easy.
        // The following things have to be taken into account:
        //  1) "deltaPhi" is a complicated function
        //  2) The support rectangle can be very small -- yet,
        //     the integral can be finite
        //  3) The support rectangle can be very large, so that
        //     various ratios can overflow integer numbers
        // Because of these complications, we will for now require
        // a full scan in phi.

        *nphiMin = 0;
        *nphiMax = nPhi;
        return *netaMin > 0 || *netaMax < nEta;
    }
}
