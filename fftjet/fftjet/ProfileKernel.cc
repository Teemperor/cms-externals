#include <cfloat>
#include <cmath>
#include <cassert>

#include "fftjet/ProfileKernel.hh"
#include "fftjet/interpolate.hh"
#include "fftjet/cdfCellNumber.hh"

// The following number is sqrt(0.6)
#define SQR3OF5 0.77459666924148337704

namespace fftjet {
    ProfileKernel::ProfileKernel(const double sx, const double sy,
                                 const int scalePower,
                                 const std::vector<double>& params)
        : AbsSymmetricKernel(sx, sy, scalePower),
          params_(params),
          normfactor_(1.0),
          cdf_(0)
    {
        normalizeProfile();
    }

    void ProfileKernel::normalizeProfile()
    {
        const unsigned n = params_.size();
        assert(n);

        switch (n)
        {
        case 1:
            assert(params_[0] > 0.0);
            break;
        case 2:
            // Assumes parabolic modeling
            assert(params_[0] >= 0.0);
            assert(params_[1] > 0.0);
            normfactor_ = 3.0/(2.0*M_PI*params_[1]);
            break;
        default:
            const double integ = gaussIntegral();
            normfactor_ /= integ;
            break;
        }
    }

    double ProfileKernel::gaussIntegral() const
    {
        const unsigned n = params_.size();
        assert(n);
        long double sum = 0.0L;
        const double step = 1.0/n;
        const double delta = SQR3OF5*step/2.0;

        // 3-point Gauss-Legendre quadrature on each interval
        for (unsigned i=0; i<n; ++i)
        {
            assert(params_[i] >= 0.0);
            const double r0 = (i+0.5)*step;
            sum += 8.0*r0*rawInterpolatedValue(r0);
            double rd = r0 - delta;
            sum += 5.0*rd*rawInterpolatedValue(rd);
            rd = r0 + delta;
            sum += 5.0*rd*rawInterpolatedValue(rd);
        }
        assert(sum > 0.0L);
        return static_cast<double>(sum)*M_PI/9.0*step;
    }

    double ProfileKernel::rawInterpolatedValue(const double r) const
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
            // This kernel turns into the linear kernel
            value = 3.0*(1.0 - r)/M_PI;
            break;

        case 2:
            // This profile is parabolic
            value = interpolate_quadratic(
                0.0, 1.0, params_[0], params_[1], 0.0, r);
            break;

        default:
            const double di = r*n;
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
                const double v0 = interpolate_quadratic(
                    -1, 1, params_[n-3], params_[ibelow], params_[n-1], delta);
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
            return value;
        else
            return 0.0;
    }

    double ProfileKernel::eval(const double rsquared) const
    {
        assert(rsquared >= 0.0);
        if (rsquared >= 1.0)
            return 0.0;
        else
            return normfactor_*rawInterpolatedValue(sqrt(rsquared));
    }

    void ProfileKernel::buildCdf()
    {
        if (cdf_ == 0)
        {
            const unsigned n = params_.size();
            assert(n);
            assert(normfactor_ > 0.0);
            cdf_ = new double[n];
            long double sum = 0.0L;
            const double step = 1.0/n;
            const double delta = SQR3OF5*step/2.0;
            const double k = M_PI/9.0*step*normfactor_;

            for (unsigned i=0; i<n-1; ++i)
            {
                const double r0 = (i+0.5)*step;
                sum += 8.0*r0*rawInterpolatedValue(r0);
                double rd = r0 - delta;
                sum += 5.0*rd*rawInterpolatedValue(rd);
                rd = r0 + delta;
                sum += 5.0*rd*rawInterpolatedValue(rd);
                cdf_[i] = static_cast<double>(sum)*k;
            }
            cdf_[n-1] = 1.0;
        }
    }

    unsigned ProfileKernel::rCellNumber(const double r, double *delta) const
    {
        assert(r >= 0.0);
        const unsigned ncells = params_.size();
        unsigned cell = static_cast<unsigned>(r*ncells);
        if (delta)
        {
            *delta = r - cell*1.0/ncells;
            if (*delta < 0.0)
                *delta = 0.0;
        }
        return cell;
    }

    double ProfileKernel::axialWeight(const double r, unsigned) const
    {
        if (cdf_ == 0)
            const_cast<ProfileKernel*>(this)->buildCdf();
        assert(r >= 0.0);
        if (r >= 1.0)
            return 1.0;
        double delta;
        const unsigned cell = rCellNumber(r, &delta);
        const double step = 1.0/params_.size();
        const double r0 = cell*step + delta/2.0;
        const double d = SQR3OF5*delta/2.0;
        double sum = 8.0*r0*rawInterpolatedValue(r0);
        double rd = r0 - d;
        sum += 5.0*rd*rawInterpolatedValue(rd);
        rd = r0 + d;
        sum += 5.0*rd*rawInterpolatedValue(rd);
        sum *= (M_PI/9.0*delta*normfactor_);
        if (cell)
            return cdf_[cell-1] + sum;
        else
            return sum;
    }

    double ProfileKernel::randomRadius(const double rnd) const
    {
        if (cdf_ == 0)
            const_cast<ProfileKernel*>(this)->buildCdf();

        const unsigned n = params_.size();
        const unsigned icell = cdfCellNumber(rnd, cdf_, n);
        const double step = 1.0/n;
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
            const double d = SQR3OF5*delta/2.0;
            const double r0 = base + delta/2.0;;
            double sum = 8.0*r0*rawInterpolatedValue(r0);
            double rd = r0 - d;
            sum += 5.0*rd*rawInterpolatedValue(rd);
            rd = r0 + d;
            sum += 5.0*rd*rawInterpolatedValue(rd);
            sum *= (M_PI/9.0*delta*normfactor_);
            if (sum > ydelta)
                xhi = half;
            else
                xlow = half;
            half = (xlow + xhi)/2.0;
        }
        return half;
    }
}
