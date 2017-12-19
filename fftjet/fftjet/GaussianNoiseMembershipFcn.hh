//=========================================================================
// GaussianNoiseMembershipFcn.hh
//
// Noise membership functor for a simple model in which the noise
// in each digitization cell is Gaussian, pedestal is 0, and energy
// spectrum of the signal does not change substantially near 0.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_GAUSSIANNOISEMEMBERSHIPFCN_HH_
#define FFTJET_GAUSSIANNOISEMEMBERSHIPFCN_HH_

#include <cmath>
#include <cassert>

#include "fftjet/SimpleFunctors.hh"

namespace fftjet {
    class GaussianNoiseMembershipFcn : public Functor2<double,double,double>
    {
    public:
        // The "prior" paramerer should be inversely proportional
        // to an apriori estimate of the signal occupancy.
        inline GaussianNoiseMembershipFcn(const double minWeight,
                                          const double prior)
            : minWeight_(minWeight), prior_(prior), sqrt2pi_(sqrt(2.0*M_PI))
        {
            assert(minWeight_ >= 0.0);
            assert(prior_ >= 0.0);
        }

        inline double operator()(const double& et, const double& sigma) const
        {
            if (prior_)
            {
                assert(et >= 0.0 && sigma > 0.0);
                const double x = et/sigma;
                return minWeight_  + prior_/sqrt2pi_/sigma*exp(-x*x/2.0);
            }
            else
                return minWeight_;
        }

    private:
        GaussianNoiseMembershipFcn();

        const double minWeight_;
        const double prior_;
        const double sqrt2pi_;
    };
}

#endif // FFTJET_GAUSSIANNOISEMEMBERSHIPFCN_HH_
