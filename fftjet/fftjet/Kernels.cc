#include <cassert>
#include <cfloat>
#include <cmath>

#include "fftjet/quartic_lib.hh"
#include "fftjet/Kernels.hh"
#include "fftjet/cdfCellNumber.hh"
#include "fftjet/SpecialFunctions.hh"

#define INVP2D_CDF_INTERVALS 64

namespace fftjet {
    SymmetricBeta::SymmetricBeta(const double sx, const double sy,
                                 const int scalePow,
                                 const std::vector<double>& params)
        : AbsSymmetricKernel(sx, sy, scalePow)
    {
        assert(params.size() == static_cast<unsigned>(nParameters()));
        n_ = params[0];
        assert(n_ >= 0.0);
    }

    double SymmetricBeta::eval(const double rsquare) const
    {
        if (rsquare >= 1.0)
            return 0.0;
        else if (n_ == 1.0)
            return 0.636619772367581343*(1.0 - rsquare);
        else if (n_ == 0.0)
            return 1.0/M_PI;
        else
            return (n_ + 1.0)/M_PI*pow(1.0 - rsquare, n_);
    }

    double SymmetricBeta::calculateResolution() const
    {
        if (n_ <= 1.0)
            return 1.0;
        else
        {
            const double rmin = sqrt(1.0/(2.0*n_ - 1.0));
            return scanForResolution(6, rmin);
        }
    }

    double SymmetricBeta::randomRadius(const double rnd) const
    {
        return sqrt(1.0 - pow(rnd, 1.0/(1.0 + n_)));
    }

    double Gauss2d::supportDistance() const
    {
        return inverseGaussCdf(1.0);
    }

    double Gauss2d::eval(const double rsquare) const
    {
        return exp(-rsquare/2.0)/(2.0*M_PI);
    }

    double Gauss2d::axialWeight(double r, unsigned) const
    {
        assert(r >= 0.0);
        return 1.0 - exp(-r*r/2.0);
    }

    double Gauss2d::randomRadius(const double rnd) const
    {
        const double maxsigma = inverseGaussCdf(1.0);
        if (rnd)
        {
            const double r = sqrt(-2.0*log(rnd));
            if (r < maxsigma)
                return r;
        }
        return maxsigma;
    }

    double Linear2d::eval(const double rsquare) const
    {
        if (rsquare >= 1.0)
            return 0.0;
        else
            return 3.0*(1.0 - sqrt(rsquare))/M_PI;
    }

    double Linear2d::randomRadius(const double rnd) const
    {
        double r = 0.0, v3[3];
        const int nroots = cubic(-1.5, 0.0, rnd/2.0, v3);
        switch (nroots)
        {
        case 1:
            r = v3[0];
            break;

        case 3:
            {
                double d, distance = DBL_MAX;
                for (unsigned i=0; i<3; ++i)
                {
                    d = fabs(v3[i] - 0.5);
                    if (d < distance)
                    {
                        r = v3[i];
                        distance = d;
                    }
                }
            }
            break;

        default:
            assert(0);
        }

        if (r < 0.0)
            r = 0.0;
        if (r > 1.0)
            r = 1.0;
        return r;
    }

    SubGauss::SubGauss(const double sx, const double sy,
                       const int scalePow, const double alpha)
        : AbsSymmetricKernel(sx, sy, scalePow),
          alpha_(alpha)
    {
        assert(alpha_ > 0.0);
        normfactor_ = calculateNormfactor(alpha_);
    }

    SubGauss::SubGauss(const double sx, const double sy,
                       const int scalePow,
                       const std::vector<double>& params)
        : AbsSymmetricKernel(sx, sy, scalePow)
    {
        assert(params.size() == static_cast<unsigned>(nParameters()));
        alpha_ = params[0];
        assert(alpha_ > 0.0);
        normfactor_ = calculateNormfactor(alpha_);
    }

    double SubGauss::eval(const double rsquared) const
    {
        return normfactor_*exp(-pow(rsquared/2.0, alpha_/2.0));
    }

    double SubGauss::calculateNormfactor(const double alpha)
    {
        return 1.0/(2.0*M_PI*Gamma((alpha + 2.0)/alpha));
    }

    double SubGauss::supportDistance() const
    {
        const double a = alpha_ > 0.01 ? alpha_ : 0.01;
        // The returned result is sqrt(2.0)*exp(log(800.0)/a);
        return 1.41421356237*exp(6.68461172767/a);
    }

    double SubGauss::randomRadius(const double rnd) const
    {
        assert(rnd >= 0.0 && rnd <= 1.0);
        if (rnd == 0.0)
            return 0.0;
        else if (rnd == 1.0)
            return supportDistance();
        else
        {
            const double sqrof2 = 1.4142135623730950488;
            return sqrof2*pow(inverseIncompleteGamma(2.0/alpha_, rnd),
                              1.0/alpha_);
        }
    }

    Huber2d::Huber2d(const double sx, const double sy,
                     const int scalePow,
                     const std::vector<double>& params)
        : AbsSymmetricKernel(sx, sy, scalePow)
    {
        assert(params.size() == static_cast<unsigned>(nParameters()));
        tWeight_ = params[0];
        setup();
    }

    Huber2d::Huber2d(const double sx, const double sy,
                     const int scalePow,
                     const double iTailWeight)
        : AbsSymmetricKernel(sx, sy, scalePow), tWeight_(iTailWeight)
    {
        setup();
    }

    double Huber2d::supportDistance() const
    {
        if (a_ == DBL_MAX)
            return 40.0;
        else
            return 800.0/a_ + a_/2.0;
    }

    double Huber2d::calculateHalfWeightRadius() const
    {
        if (tWeight_ == 0.0)
            return 1.177410022515474691;
        else if (tWeight_ < 0.5)
        {
            const double l = log(1.0 - 0.5/(1.0-tWeight_)*
                                 (1.0-exp(-a_*a_/2.0)));
            return sqrt(-2.0*l);
        }
        else if (tWeight_ == 0.5)
            return a_;
        else
        {
            const double eps = 8.0*DBL_EPSILON;
            const double wratio = 0.5/tWeight_;
            double xmin = a_;
            double xmax = 2.0*a_;
            while ((1 + a_*xmax)/(1 + a_*a_)*exp(a_*(a_ - xmax)) >= wratio)
            {
                xmin = xmax;
                xmax *= 2.0;
            }
            while (2.0*(xmax - xmin)/(xmax + xmin) > eps)
            {
                const double x = (xmax + xmin)/2.0;
                const double f = (1 + a_*x)/(1 + a_*a_)*exp(a_*(a_ - x));
                if (f < wratio)
                    xmax = x;
                else
                    xmin = x;
            }
            return (xmax + xmin)/2.0;
        }
    }

    double Huber2d::weight_(const double asquared) const
    {
        return (1.0 + asquared)/(1.0 + asquared*exp(asquared/2.0));
    }

    void Huber2d::setup()
    {
        assert(tWeight_ >= 0.0 && tWeight_ < 1.0);
        if (tWeight_ == 0.0)
        {
            // Pure Gaussian
            a_ = DBL_MAX;
            normfactor_ = 1.0/(2.0*M_PI);
        }
        else
        {
            // Solve the equation for "a" by bisection
            const double eps = 8.0*DBL_EPSILON;
            double b = -2.0*log(tWeight_);
            assert(b > 0.0);
            assert(weight_(b) >= tWeight_);
            double c = 2.0*b;
            while (weight_(c) >= tWeight_)
                c *= 2.0;
            while ((c - b)/b > eps)
            {
                const double half = (c + b)/2.0;
                if (weight_(half) >= tWeight_)
                    b = half;
                else
                    c = half;
            }
            a_ = sqrt(b);
            normfactor_ = 1.0/(2.0*M_PI*(1.0 + exp(-b/2.0)/b));
        }
    }

    double Huber2d::eval(const double rsquared) const
    {
        const double r = sqrt(rsquared);
        if (r <= a_)
            return normfactor_*exp(-rsquared/2.0);
        else
            return normfactor_*exp(a_*(a_/2.0-r));
    }

    double Huber2d::randomRadius(const double rnd) const
    {
        if (tWeight_ == 0.0)
        {
            // Pure Gaussian
            const double maxsigma = inverseGaussCdf(1.0);
            const double del = 1.0 - rnd;
            if (del)
            {
                const double r = sqrt(-2.0*log(del));
                if (r < maxsigma)
                    return r;
            }
            return maxsigma;
        }
        else if (rnd <= 1.0 - tWeight_)
        {
            // Gaussian part of Huber
            const double del = 1.0 - rnd/(2.0*M_PI*normfactor_);
            assert(del > 0.0);
            const double r = sqrt(-2.0*log(del));
            if (r < a_)
                return r;
            else
                return a_;
        }
        else
        {
            // Exponential part of Huber
            const double del = rnd - (1.0 - tWeight_);
            const double b = a_*a_;
            const double expb = exp(b);
            const double rhs = (1.0 + b - del*(1.0 + b*exp(b/2.0)))/expb;
            if (rhs <= 0.0)
                return 800.0/a_ + a_/2.0;
            // Bracket the root
            double left = -log(rhs);
            if (left < b)
                left = b;
            assert(exp(-left)*(1.0 + left) >= rhs);
            double right = 2*left;
            while (exp(-right)*(1.0 + right) >= rhs)
            {
                left = right;
                right *= 2.0;
            }
            // Solve by bisection
            const double eps = 8.0*DBL_EPSILON;
            double half = (right + left)/2.0;
            while ((right - left)/half > eps)
            {
                if (exp(-half)*(1.0 + half) >= rhs)
                    left = half;
                else
                    right = half;
                half = (right + left)/2.0;
            }
            return half/a_;
        }
    }

    DeltaFunctionKernel::DeltaFunctionKernel(double, double, int,
                                             const std::vector<double>& params)
    {
        assert(params.size() == static_cast<unsigned>(nParameters()));
        value_ = params[0];
    }

    double DeltaFunctionKernel::operator()(
        const double x, const double y, double) const
    {
        if (x == 0.0 && y == 0.0)
            assert(!"Attempt to evaluate 2d delta function at (0, 0)");
        return 0.0;
    }

    void DeltaFunctionKernel::supportRectangle(
        double, KernelSupportRectangle *r) const
    {
        const double eps = 1.0/sqrt(DBL_MAX);
        r->xmin = -eps;
        r->xmax = eps;
        r->ymin = -eps;
        r->ymax = eps;
    }

    double DeltaFunctionKernel::rectangleAverage(
        const double x, const double y, double,
        const double dx, const double dy) const
    {
        const double halfx(fabs(dx/2.0));
        if (x - halfx <= 0.0 && x + halfx > 0.0)
        {
            const double halfy(fabs(dy/2.0));
            if (y - halfy <= 0.0 && y + halfy > 0.0)
                return value_/dx/dy;
        }
        return 0.0;
    }

    void DeltaFunctionKernel::random(
        double, double, double,
        double* px, double* py) const
    {
        *px = 0.0;
        *py = 0.0;
    }

    InvPower2d::InvPower2d(const double sx, const double sy,
                           const int scalePow,
                           const std::vector<double>& params)
        : AbsSymmetricKernel(sx, sy, scalePow), cdf_(0)
    {
        assert(params.size() == static_cast<unsigned>(nParameters()));
        n_ = params[0];
        setup();
    }

    InvPower2d::InvPower2d(const double sx, const double sy,
                           const int scalePow,
                           const double power)
        : AbsSymmetricKernel(sx, sy, scalePow), n_(power), cdf_(0)
    {
        setup();
    }

    void InvPower2d::setup()
    {
        assert(n_ >= 1.05);
        normfactor_ = sin(M_PI/n_)*n_/(M_PI*M_PI);

        lowbound_ = pow(5.0, -1.0/(2.0*n_));
        lowweight_ = headIntegral(lowbound_);
        upbound_ = 1.0/lowbound_;
        upweight_ = tailIntegral(upbound_);
    }

    // The following function returns the integral
    // from x to infinity of N (2 Pi y)/(1 + y^(2n)) dy
    // for sufficiently large values of x (so that
    // the series actually converge). N here is the kernel
    // normalization factor.
    double InvPower2d::tailIntegral(const double x) const
    {
        assert(x > 1.0);
        long double sum = 0.0L;
        for (unsigned k=0; k<100000; ++k)
        {
            const long double oldsum = sum;
            sum += (k % 2 ? -1.0 : 1.0)*pow(x, -2.0*n_*(k+1))/(2*n_*(k+1)-2);
            if (sum == oldsum)
                break;
        }
        return 2.0*M_PI*x*x*sum*normfactor_;
    }

    // The following function returns the integral
    // from 0 to x of N (2 Pi y)/(1 + y^(2n)) dy
    // for sufficiently small values of x (so that
    // the series actually converge). N here is the kernel
    // normalization factor.
    double InvPower2d::headIntegral(const double x) const
    {
        assert(x < 1.0);
        long double sum = 0.0L;
        for (unsigned k=0; k<100000; ++k)
        {
            const long double oldsum = sum;
            sum += (k % 2 ? -1.0 : 1.0)*pow(x, 2*k*n_)/(2*n_*k+2);
            if (sum == oldsum)
                break;
        }
        return 2.0*M_PI*x*x*sum*normfactor_;
    }

    double InvPower2d::eval(const double rsquared) const
    {
        return normfactor_/(1.0 + pow(rsquared, n_));
    }

    double InvPower2d::supportDistance() const
    {
        return exp(400.0/n_);
    }

    double InvPower2d::randomRadius(const double rnd) const
    {
        if (rnd == 0.0)
            return 0.0;
        if (rnd == 1.0)
            return exp(400.0/n_);

        if (n_ == 2.0)
            return sqrt(tan(rnd*M_PI/2.0));

        const double targetEps = 8.0*DBL_EPSILON;
        if (rnd <= lowweight_)
        {
            double xmin = 0.0;
            double xmax = lowbound_;
            while ((xmax - xmin)/lowbound_ > targetEps)
            {
                const double xtry = (xmin + xmax)/2.0;
                if (headIntegral(xtry) > rnd)
                    xmax = xtry;
                else
                    xmin = xtry;
            }
            return (xmin + xmax)/2.0;
        }
        else if (rnd >= 1.0 - upweight_)
        {
            const double tail = 1.0 - rnd;
            double xmin = upbound_;
            double xmax = xmin;
            double integ = upweight_;
            while (integ > tail)
            {
                xmin = xmax;
                xmax *= 2.0;
                integ = tailIntegral(xmax);
            }
            while (2.0*(xmax-xmin)/(xmax + xmin) > targetEps)
            {
                const double x = (xmax + xmin)/2.0;
                if (tailIntegral(x) > tail)
                    xmin = x;
                else
                    xmax = x;
            }
            return (xmax + xmin)/2.0;
        }
        else
        {
            if (cdf_ == 0)
                const_cast<InvPower2d*>(this)->buildCdf();
            const double cdffactor = 1.0/(1.0 - (upweight_ + lowweight_));
            const double newr = (rnd - lowweight_)*cdffactor;
            const unsigned icell = cdfCellNumber(
                newr, cdf_, INVP2D_CDF_INTERVALS);
            const double step = (upbound_ - lowbound_)/INVP2D_CDF_INTERVALS;
            const double xcdf = lowbound_ + icell*step;
            const double remainder = icell ? newr - cdf_[icell-1] : newr;
            double deltaMin = 0.0;
            double deltaMax = step;
            while ((deltaMax - deltaMin)/step > targetEps)
            {
                const double delta = (deltaMin + deltaMax)/2.0;
                const double integ = cdffactor*annulusIntegral(
                    xcdf + delta/2.0, delta);
                if (integ > remainder)
                    deltaMax = delta;
                else
                    deltaMin = delta;
            }
            return xcdf + (deltaMin + deltaMax)/2.0;
        }
    }

    void InvPower2d::buildCdf()
    {
        if (cdf_ == 0)
        {
            cdf_ = new double[INVP2D_CDF_INTERVALS];
            const double step = (upbound_ - lowbound_)/INVP2D_CDF_INTERVALS;
            long double sum = 0.0L;
            for (unsigned i=0; i<INVP2D_CDF_INTERVALS; ++i)
            {
                sum += annulusIntegral(lowbound_ + (i+0.5)*step, step);
                cdf_[i] = sum;
            }
            const double dsum = static_cast<double>(sum);
            for (unsigned i=0; i<INVP2D_CDF_INTERVALS-1; ++i)
                cdf_[i] /= dsum;
            cdf_[INVP2D_CDF_INTERVALS - 1] = 1.0;
        }
    }

    double InvPower2d::annulusIntegral(const double r0,
                                       const double width) const
    {
        double sum = 0.0;
        // 0.77459666924148337704 is actually sqrt(0.6)
        const double delta = 0.77459666924148337704*width/2.0;
        sum += 8.0*r0*eval(r0*r0);
        double rd = r0 - delta;
        sum += 5.0*rd*eval(rd*rd);
        rd = r0 + delta;
        sum += 5.0*rd*eval(rd*rd);
        return sum*M_PI/9.0*width;
    }

    double InvPower2d::calculateHalfWeightRadius() const
    {
        if (n_ == 2.0)
            return 1.0;
        else
            return randomRadius(0.5);
    }

    double InvPower2d::axialWeight(const double r, unsigned) const
    {
        if (r <= lowbound_)
            return headIntegral(r);
        else if (r >= upbound_)
            return 1.0 - tailIntegral(r);
        else
        {
            if (cdf_ == 0)
                const_cast<InvPower2d*>(this)->buildCdf();
            const double step = (upbound_ - lowbound_)/INVP2D_CDF_INTERVALS;
            const unsigned icell = static_cast<unsigned>((r - lowbound_)/step);
            const double delta = r - lowbound_ - icell*step;
            const double add = icell ? 
                cdf_[icell-1]*(1.0 -(upweight_+lowweight_)) : 0.0;
            return lowweight_ + add + annulusIntegral(
                lowbound_ + icell*step + delta/2.0, delta);
        }
    }
}
