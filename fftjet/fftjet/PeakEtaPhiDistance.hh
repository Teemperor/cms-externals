//=========================================================================
// PeakEtaPhiDistance.hh
//
// A simple implementation of the AbsDistanceCalculator class for use with
// ProximityClusteringTree<Peak, ...>.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_PEAKETAPHIDISTANCE_HH_
#define FFTJET_PEAKETAPHIDISTANCE_HH_

#include <cassert>

#include "fftjet/AbsDistanceCalculator.hh"
#include "fftjet/Peak.hh"

namespace fftjet {
    class PeakEtaPhiDistance : public AbsDistanceCalculator<Peak>
    {
    public:
        inline explicit PeakEtaPhiDistance(const double etaToPhiBandwidthRatio=1.0)
            : bwEta(sqrt(etaToPhiBandwidthRatio)), bwPhi(1.0/bwEta)
        {
            assert(etaToPhiBandwidthRatio > 0.0);
        }
        virtual ~PeakEtaPhiDistance() {}

    private:
        double distanceBetweenClusters(
            double, const Peak& smallerJet,
            double, const Peak& biggerJet) const;

        const double bwEta;
        const double bwPhi;
    };
}

#endif // FFTJET_PEAKETAPHIDISTANCE_HH_
