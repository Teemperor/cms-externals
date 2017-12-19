#include <cmath>
#include <cassert>

#include "fftjet/PhiKernels.hh"
#include "fftjet/interpolate.hh"
#include "fftjet/cdfCellNumber.hh"
#include "fftjet/SpecialFunctions.hh"

#define ONEOVERSQR2PI 0.398942280401432678

// The following number is 1/(2 sqrt(3))
#define AVE_DELTA 0.28867513459481288225

namespace fftjet {
    double PhiGauss::unscaledPhiFcn(const double phi) const
    {
        return ONEOVERSQR2PI*exp(-phi*phi/2.0);
    }

    void PhiGauss::unscaledPhiSupport(double *phimin, double *phimax) const
    {
        *phimax = inverseGaussCdf(1.0);
        *phimin = -*phimax;
    }

    double PhiGauss::unscaledPhiRandom(const double r) const
    {
        return inverseGaussCdf(r);
    }

    void PhiProfileKernel::unscaledPhiSupport(double *phimin,
                                              double *phimax) const
    {
        *phimin = -M_PI/2.0;
        *phimax = M_PI/2.0;
    }

    double PhiProfileKernel::unscaledPhiFcn(const double yin) const
    {
        const double r = fabs(yin);
        if (r < M_PI/2.0)
        {
            double value;
            const unsigned n = params_.size();
            switch (n)
            {
            case 0:
                assert(0);
                value = 0.0;
                break;

            case 1:
                // Linear function
                value = lin_interpolate_1d(0.0, M_PI/2.0, params_[0], 0.0, r);
                break;

            case 2:
                // Quadratic function
                value = interpolate_quadratic(0.0, M_PI/2.0, params_[0],
                                              params_[1], 0.0, r);
                break;

            default:
                const double di = r/(M_PI/2.0)*n;
                unsigned ibelow = static_cast<unsigned>(di);
                if (ibelow >= n)
                    ibelow = n - 1;
                const double delta = di - ibelow;
                if (ibelow == 0)
                    value = interpolate_quadratic(
                        0, 2, params_[0], params_[1], params_[2], delta);
                else if (ibelow == n - 1)
                    value = interpolate_quadratic(
                        -1, 1, params_[n-2], params_[ibelow], 0.0, delta);
                else if (ibelow == n - 2)
                {
                    const double v0 = interpolate_quadratic(-1, 1,
                        params_[n-3], params_[ibelow], params_[n-1], delta);
                    const double v1 = interpolate_quadratic(
                        0, 2, params_[ibelow], params_[n-1], 0.0, delta);
                    value = v0*(1.0 - delta) + v1*delta;
                }
                else
                {
                    // Cubic interpolation scheme
                    const double v1 = interpolate_quadratic(
                        -1, 1, params_[ibelow-1],
                        params_[ibelow], params_[ibelow+1], delta);
                    const double v2 = interpolate_quadratic(
                        0, 2, params_[ibelow], params_[ibelow+1],
                        params_[ibelow+2], delta);
                    const double w1 = (2.0 - delta)/3.0;
                    const double w2 = (delta + 1.0)/3.0;
                    value = v1*w1 + v2*w2;
                }
                break;
            }
            if (value > 0.0)
                return value/normfactor_;
            else
                return 0.0;
        }
        else
            return 0.0;
    }

    PhiProfileKernel::PhiProfileKernel(double, const double sy,
                                       const int scalePower,
                                       const std::vector<double>& params)
        : AbsScalablePhiKernel(sy, scalePower),
          params_(params),
          cdf_(0),
          normfactor_(1.0)
    {
        const unsigned nsteps = params_.size();
        assert(nsteps);
        cdf_ = new double[nsteps];
        const double step = M_PI/2.0/nsteps;
        const double delta = AVE_DELTA*step;
        long double sum = 0.0L;
        for (unsigned i=0; i<nsteps; ++i)
        {
            assert(params[i] >= 0.0);
            const double center = (i + 0.5)*step;
            sum += unscaledPhiFcn(center - delta);
            sum += unscaledPhiFcn(center + delta);
            cdf_[i] = static_cast<double>(sum);
        }
        const unsigned nm1 = nsteps-1;
        const double normfactor = cdf_[nm1];
        assert(normfactor > 0.0);
        for (unsigned i=0; i<nm1; ++i)
            cdf_[i] /= normfactor;
        cdf_[nm1] = 1.0;
        normfactor_ = normfactor*step;
    }

    double PhiProfileKernel::unscaledPhiRandom(const double rnd) const
    {
        if (rnd >= 0.5)
            return randomDistance(2.0*rnd-1.0);
        else
            return -randomDistance(1.0-2.0*rnd);
    }

    unsigned PhiProfileKernel::rCellNumber(const double r, double *delta) const
    {
        assert(r >= 0.0);
        const unsigned ncells = params_.size();
        unsigned cell = static_cast<unsigned>(r*ncells/(M_PI/2.0));
        if (r < M_PI/2.0)
            if (cell >= ncells)
                cell = ncells - 1;
        if (delta)
        {
            *delta = r - cell*(M_PI/2.0/ncells);
            if (*delta < 0.0)
                *delta = 0.0;
        }
        return cell;
    }

    double PhiProfileKernel::randomDistance(const double rnd) const
    {
        const unsigned ncells = params_.size();
        const double step = M_PI/2.0/ncells;
        const unsigned icell = cdfCellNumber(rnd, cdf_, ncells);
        const double base = icell*step;
        const double ydelta = rnd - (icell ? cdf_[icell-1] : 0.0);
        assert(ydelta >= 0.0);
        const double eps = 8.0*DBL_EPSILON;
        double xlow = base;
        double xhi = (icell+1)*step;
        double half = (xlow + xhi)/2.0;
        while ((xhi - xlow)/(half + step) > eps)
        {
            const double delta = half - base;
            const double yshift = AVE_DELTA*delta;
            const double y = base+delta/2.0;
            const double remainder = delta*
                (unscaledPhiFcn(y - yshift) + 
                 unscaledPhiFcn(y + yshift));
            if (remainder > ydelta)
                xhi = half;
            else
                xlow = half;
            half = (xlow + xhi)/2.0;
        }
        return half;
    }
}
