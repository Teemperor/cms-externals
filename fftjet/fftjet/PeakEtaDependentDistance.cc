#include <cmath>
#include <cassert>

#include "fftjet/PeakEtaDependentDistance.hh"

namespace fftjet {
    PeakEtaDependentDistance::PeakEtaDependentDistance(
        const LinearInterpolator1d& tableOfEtaToPhiBandwidthRatios)
        : table_(tableOfEtaToPhiBandwidthRatios)
    {
        // Make sure the bandwidth ratio stays positive for any argument
        const unsigned n = table_.nx();
        const double* data = table_.getData();
        for (unsigned i=0; i<n; ++i)
            assert(data[i] > 0.0);
        assert(table_.fLow() > 0.0);
        assert(table_.fHigh() > 0.0);
    }

    double PeakEtaDependentDistance::distanceBetweenClusters(
        const double smallerScale, const Peak& smallerJet,
        const double biggerScale, const Peak& biggerJet) const
    {
        // Which location in the interpolation table should we pick?
        // We will weight eta values by the product of scale and magnitude
        // (this quantity, for an appropriate scale, is a scale-invariant
        // measure of the peak volume).
        const double smallerEta = smallerJet.eta();
        const double biggerEta = biggerJet.eta();
        double smallerMag = smallerJet.magnitude();
        if (smallerMag < 0.0)
            smallerMag = 0.0;
        double biggerMag = biggerJet.magnitude();
        if (biggerMag < 0.0)
            biggerMag = 0.0;
        const double denom = smallerMag*smallerScale + biggerMag*biggerScale;
        const double lookupEta = denom > 0.0 ?
            (smallerMag*smallerScale*smallerEta + 
             biggerMag*biggerScale*biggerEta)/denom : 
            (smallerEta + biggerEta)/2.0;
        const double bwEta = sqrt(table_(lookupEta));
        const double bwPhi = 1.0/bwEta;

        // Now, proceed with the normal calculation
        const double deta = (smallerEta - biggerEta)/bwEta;
        double dphi = smallerJet.phi() - biggerJet.phi();
        if (dphi > M_PI)
            dphi -= 2.0*M_PI;
        else if (dphi < -M_PI)
            dphi += 2.0*M_PI;
        dphi /= bwPhi;
        return sqrt(deta*deta + dphi*dphi);
    }
}
