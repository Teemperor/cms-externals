//=========================================================================
// PeakEtaDependentDistance.hh
//
// An implementation of the AbsDistanceCalculator class similar to
// PeakEtaPhiDistance but in which the ratio between eta and phi
// bandwidth values is given by an eta-dependent interpolation table.
// Note that it is possible to (mis)specify the interpolation table
// in such a way that this function is no longer a pseudo-metric.
// Therefore, use at your own risk.
//
// I. Volobouev
// June 2010
//=========================================================================

#ifndef FFTJET_PEAKETADEPENDENTDISTANCE_HH_
#define FFTJET_PEAKETADEPENDENTDISTANCE_HH_

#include "fftjet/AbsDistanceCalculator.hh"
#include "fftjet/LinearInterpolator1d.hh"
#include "fftjet/Peak.hh"

namespace fftjet {
    class PeakEtaDependentDistance : public AbsDistanceCalculator<Peak>
    {
    public:
        explicit PeakEtaDependentDistance(
            const LinearInterpolator1d& tableOfEtaToPhiBandwidthRatios);

        virtual ~PeakEtaDependentDistance() {}

    private:
        PeakEtaDependentDistance();

        double distanceBetweenClusters(
            double smallerScale, const Peak& smallerJet,
            double biggerScale, const Peak& biggerJet) const;

        LinearInterpolator1d table_;
    };
}

#endif // FFTJET_PEAKETADEPENDENTDISTANCE_HH_
