#include <cmath>
#include "fftjet/PeakEtaPhiDistance.hh"

namespace fftjet {
    double PeakEtaPhiDistance::distanceBetweenClusters(
        double, const Peak& smallerJet,
        double, const Peak& biggerJet) const
    {
        const double deta = (smallerJet.eta() - biggerJet.eta())/bwEta;
        double dphi = smallerJet.phi() - biggerJet.phi();
        if (dphi > M_PI)
            dphi -= 2.0*M_PI;
        else if (dphi < -M_PI)
            dphi += 2.0*M_PI;
        dphi /= bwPhi;
        return sqrt(deta*deta + dphi*dphi);
    }
}
