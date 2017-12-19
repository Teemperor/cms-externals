#include <cmath>
#include <cassert>
#include <cfloat>

#include "fftjet/Kernels1d.hh"
#include "fftjet/SpecialFunctions.hh"

namespace fftjet {
    double Gauss1d::calculate(const double x) const
    {
        static const double sqrtwopi = 2.5066282746310005;
        return exp(-x*x/2.0)/sqrtwopi;
    }

    KernelSupportInterval Gauss1d::unscaledInterval() const
    {
        return KernelSupportInterval(inverseGaussCdf(0.0),
                                     inverseGaussCdf(1.0));
    }

    double Gauss1d::unscaledRandom(const double r1) const
    {
        return inverseGaussCdf(r1);
    }

    SymmetricBeta1d::SymmetricBeta1d(double sx, int scalePow, double power)
        : AbsScalableKernel1d(sx, scalePow),
          n_(power)
    {
        norm_ = calculateNorm();
    }

    SymmetricBeta1d::SymmetricBeta1d(double sx, int scalePow,
                                     const std::vector<double>& p)
        : AbsScalableKernel1d(sx, scalePow)
    {
        assert(static_cast<int>(p.size()) == nParameters());
        n_ = p[0];
        norm_ = calculateNorm();
    }

    double SymmetricBeta1d::calculateNorm() const
    {
        static const double normcoeffs[11] = {
            0.5, 0.75, 0.9375, 1.09375, 1.23046875, 1.353515625,
            1.46630859375, 1.571044921875, 1.6692352294921875,
            1.76197052001953125, 1.85006904602050781};

        assert(n_ > -1.0);
        
        const int intpow = static_cast<int>(floor(n_));
        if (static_cast<double>(intpow) == n_ &&
            intpow >= 0 && intpow <= 10)
            return normcoeffs[intpow];
        else
            return Gamma(1.5 + n_)/sqrt(M_PI)/Gamma(1.0 + n_);
    }

    double SymmetricBeta1d::calculate(const double x) const
    {
        const double oneminusrsq = 1.0 - x*x;
        if (oneminusrsq <= 0.0)
            return 0.0;
        else
            return norm_*pow(oneminusrsq, n_);
    }

    KernelSupportInterval SymmetricBeta1d::unscaledInterval() const
    {
        return KernelSupportInterval(-1.0, 1.0);
    }

    double SymmetricBeta1d::unscaledRandom(const double r1) const
    {
        double r;
        if (n_ == 0.0)
            r = r1*2.0 - 1.0;
        else
            r = 2.0*inverseIncompleteBeta(n_+1.0, n_+1.0, r1) - 1.0;
        if (r < -1.0)
            r = -1.0;
        else if (r > 1.0)
            r = 1.0;
        return r;
    }

    DeltaFunction1d::DeltaFunction1d(double, int,
                                     const std::vector<double>& params)
    {
        assert(params.size() == static_cast<unsigned>(nParameters()));
        value_ = params[0];
    }

    double DeltaFunction1d::operator()(const double x, double) const
    {
        if (x == 0.0)
            assert(!"Attempt to evaluate 1d delta function at 0");
        return 0.0;
    }

    KernelSupportInterval DeltaFunction1d::supportInterval(double) const
    {
        const double eps = 1.0/sqrt(DBL_MAX);
        return KernelSupportInterval(-eps, eps);
    }

    double DeltaFunction1d::intervalAverage(
        const double x, const double /* scale */, const double dx) const
    {
        const double halfx(fabs(dx/2.0));
        if (x - halfx <= 0.0 && x + halfx > 0.0)
            return value_/dx;
        else
            return 0.0;
    }

    double DeltaFunction1d::random(double, double) const
    {
        return 0.0;
    }
}
