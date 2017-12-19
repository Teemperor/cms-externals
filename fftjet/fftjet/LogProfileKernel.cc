#include <cfloat>
#include <cmath>
#include <cassert>

#include "fftjet/LogProfileKernel.hh"
#include "fftjet/interpolate.hh"
#include "fftjet/cdfCellNumber.hh"

namespace fftjet {
    double LogProfileKernel::annulusIntegral(const double middleR,
                                             const double width) const
    {
        // It makes sense to use a relatively high order formula here
        // because exponent has an infinite number of derivatives
        const double x[] = {
            -0.45308992296933199640,
            -0.26923465505284154552,
            0.0,
            0.26923465505284154552,
            0.45308992296933199640
        };
        const double w[] = {
            0.11846344252809454376,
            0.23931433524968323402,
            0.28444444444444444444,
            0.23931433524968323402,
            0.11846344252809454376,
        };
        double sum = 0.0;
        for (unsigned i=0; i<sizeof(x)/sizeof(x[0]); ++i)
        {
            const double r = middleR + width*x[i];
            assert(r >= 0.0 && r <= maxRadius_);
            sum += w[i]*r*eval(r*r);
        }
        return width*2.0*M_PI*sum;
    }

    void LogProfileKernel::normalize()
    {
        normfactor_ = 1.0;
        const unsigned n = params_.size();
        const double step = maxRadius_/n;
        const double halfstep = step/2.0;
        const double rmax = maxRadius_ - halfstep;

        // Calculate the integral from 0 to rmax
        long double sum = annulusIntegral(halfstep/2.0, halfstep);
        for (unsigned i=1; i<n; ++i)
            sum += annulusIntegral(step*i, step);
        weightBeforeMR_ = static_cast<double>(sum);

        // Calculate the integral from rmax to infinity
        const double slope = (params_[n-2] - params_[n-1])/step;
        const double fmax = eval(rmax*rmax);
        sum += 2.0*M_PI*fmax*(1.0 + slope*rmax)/slope/slope;

        normfactor_ = 1.0/static_cast<double>(sum);
        weightBeforeMR_ /= static_cast<double>(sum);
    }

    double LogProfileKernel::axialWeight(const double r, unsigned) const
    {
        assert(r >= 0.0);
        const unsigned n = params_.size();
        const double step = maxRadius_/n;
        const double halfstep = step/2.0;
        const double rlim = maxRadius_ - halfstep;
        if (r >= rlim)
        {
            const double slope = (params_[n-2] - params_[n-1])/step;
            const double fmax = eval(r*r);
            return 1.0 - 2.0*M_PI*fmax*(1.0 + slope*r)/slope/slope;
        }
        else if (r < halfstep)
        {
            return annulusIntegral(r/2.0, r);
        }
        else
        {
            if (cdf_ == 0)
                const_cast<LogProfileKernel*>(this)->buildCdf();
            const unsigned icdf = static_cast<unsigned>((r - halfstep)/step);
            assert(icdf <= n - 1);
            if (icdf == n - 1)
                return weightBeforeMR_;
            const double delta = r - halfstep - icdf*step;
            return cdf_[icdf] + annulusIntegral(
                halfstep + icdf*step + delta/2.0, delta);
        }
    }

    void LogProfileKernel::buildCdf()
    {
        if (cdf_ == 0)
        {
            const unsigned n = params_.size();
            assert(n);
            assert(normfactor_ > 0.0);
            cdf_ = new double[n];
            const double step = maxRadius_/n;
            const double halfstep = step/2.0;

            long double sum = annulusIntegral(halfstep/2.0, halfstep);
            cdf_[0] = sum;
            for (unsigned i=1; i<n; ++i)
            {
                sum += annulusIntegral(step*i, step);
                cdf_[i] = sum;
            }
            assert(fabs(weightBeforeMR_ - cdf_[n-1]) < 1.e-10);
        }
    }

    LogProfileKernel::LogProfileKernel(const double sx, const double sy,
                                       const int scalePower,
                                       const std::vector<double>& params,
                                       const double rmax)
        : AbsSymmetricKernel(sx, sy, scalePower),
          params_(params),
          maxRadius_(rmax),
          normfactor_(0.0),
          weightBeforeMR_(0.0),
          support_(0.0),
          cdf_(0)
    {
        assert(maxRadius_ > 0.0);

        // Must have at least two values in the profile
        // in order to determine the decay slope
        const unsigned n = params_.size();
        assert(n > 1);

        // Must decay exponentially beyond rmax
        assert(params_[n-1] < params_[n-2]);

        normalize();
        support_ = estimateSupport();
        assert(support_ > 0.0);
    }

    double LogProfileKernel::eval(const double rsquared) const
    {
        assert(rsquared >= 0.0);
        assert(normfactor_ > 0.0);
        const double r = sqrt(rsquared);
        const unsigned n = params_.size();
        const double step = maxRadius_/n;
        const double halfstep = step/2.0;
        const double rmax = maxRadius_ - halfstep;
        if (r >= rmax)
        {
            const double logval = lin_interpolate_1d(
                rmax-step, rmax, params_[n-2], params_[n-1], r);
            return normfactor_*exp(logval);
        }
        else
        {
            // n > 1 from the assert in the constructor
            double logval;
            switch (n)
            {
            case 2:
                logval = lin_interpolate_1d(
                    halfstep, rmax, params_[0], params_[1], r);
                break;

            case 3:
                logval = interpolate_quadratic(
                    halfstep, rmax, params_[0], params_[1], params_[2], r);
                break;

            default:
                if (r <= step + halfstep)
                {
                    logval = interpolate_quadratic(
                        halfstep, halfstep + 2.0*step,
                        params_[0], params_[1], params_[2], r);
                }
                else if (r >= rmax - step)
                {
                    logval = interpolate_quadratic(
                        rmax - 2.0*step, rmax,
                        params_[n-3], params_[n-2], params_[n-1], r);
                }
                else
                {
                    const unsigned bin = static_cast<unsigned>(r/step);
                    unsigned basebin;
                    if (bin <= 1)
                        basebin = 1;
                    else if (bin >= n-2)
                        basebin = n-3;
                    else if (r - bin*step >= halfstep)
                        basebin = bin;
                    else
                        basebin = bin - 1;
                    const double rbase = halfstep + basebin*step;
                    const double v1 = interpolate_quadratic(
                        rbase-step, rbase+step, params_[basebin-1],
                        params_[basebin], params_[basebin+1], r);
                    const double v2 = interpolate_quadratic(
                        rbase, rbase+2*step, params_[basebin],
                        params_[basebin+1], params_[basebin+2], r);
                    const double del = (r - rbase)/step;
                    const double w1 = (2.0 - del)/3.0;
                    const double w2 = (del + 1.0)/3.0;
                    logval = v1*w1 + v2*w2;
                }
                break;
            }
            return normfactor_*exp(logval);
        }
    }

    double LogProfileKernel::estimateSupport() const
    {
        const unsigned n = params_.size();
        const double step = maxRadius_/n;
        const double halfstep = step/2.0;
        const double rmax = maxRadius_ - halfstep;
        return lin_interpolate_1d(params_[n-2], params_[n-1],
                                  rmax - step, rmax,
                                  -800.0 - log(normfactor_));
    }

    double LogProfileKernel::randomRadius(const double rnd) const
    {
        if (rnd <= 0.0)
            return 0.0;
        if (rnd >= 1.0)
            return support_;

        const unsigned n = params_.size();
        const double step = maxRadius_/n;
        const double halfstep = step/2.0;
        const double rlim = maxRadius_ - halfstep;
        const double targetEps = 8.0*DBL_EPSILON;

        if (rnd >= weightBeforeMR_)
        {
            const double tail = 1.0 - rnd;
            const double slope = (params_[n-2] - params_[n-1])/step;
            double rmin = rlim;
            double rmax = rmin;
            double integ = 2.0*M_PI*eval(rmax*rmax)*
                (1.0 + slope*rmax)/slope/slope;
            while (integ > tail)
            {
                rmin = rmax;
                rmax *= 2.0;
                integ = 2.0*M_PI*eval(rmax*rmax)*
                    (1.0 + slope*rmax)/slope/slope;
            }
            while (2.0*(rmax-rmin)/(rmax+rmin) > targetEps)
            {
                double r = (rmax + rmin)/2.0;
                integ = 2.0*M_PI*eval(r*r)*(1.0 + slope*r)/slope/slope;
                if (integ > tail)
                    rmin = r;
                else
                    rmax = r;
            }
            return (rmax + rmin)/2.0;
        }
        else
        {
            if (cdf_ == 0)
                const_cast<LogProfileKernel*>(this)->buildCdf();
            if (rnd < cdf_[0])
            {
                double deltaMin = 0.0;
                double deltaMax = halfstep;
                while ((deltaMax - deltaMin)/halfstep > targetEps)
                {
                    const double delta = (deltaMin + deltaMax)/2.0;
                    if (annulusIntegral(delta/2.0, delta) > rnd)
                        deltaMax = delta;
                    else
                        deltaMin = delta;
                }
                return (deltaMin + deltaMax)/2.0;
            }
            else if (rnd >= cdf_[n-1])
            {
                return rlim;
            }
            else
            {
                const unsigned icdf = cdfCellNumber(rnd, cdf_, n) - 1;
                const double rbase = halfstep + icdf*step;
                const double deltaCdf = rnd - cdf_[icdf];
                double deltaMin = 0.0;
                double deltaMax = step;
                while ((deltaMax - deltaMin)/rbase > targetEps)
                {
                    const double delta = (deltaMin + deltaMax)/2.0;
                    if (annulusIntegral(rbase + delta/2.0, delta) > deltaCdf)
                        deltaMax = delta;
                    else
                        deltaMin = delta;
                }
                return rbase + (deltaMin + deltaMax)/2.0;
            }
        }
    }
}
