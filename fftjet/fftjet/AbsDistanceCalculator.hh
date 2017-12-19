//=========================================================================
// AbsDistanceCalculator.hh
//
// Interface class for the distance calculator in the proximity-matched
// clustering tree
//
// I. Volobouev
// April 2008
//=========================================================================

#ifndef FFTJET_ABSDISTANCECALCULATOR_HH_
#define FFTJET_ABSDISTANCECALCULATOR_HH_

namespace fftjet {
    template<typename Cluster>
    class AbsDistanceCalculator
    {
    public:
        virtual ~AbsDistanceCalculator() {}

        inline double operator()(const double scale1,
                                 const Cluster& jet1,
                                 const double scale2,
                                 const Cluster& jet2) const
        {
            if (scale1 < scale2)
                return distanceBetweenClusters(scale1, jet1, scale2, jet2);
            else
                return distanceBetweenClusters(scale2, jet2, scale1, jet1);
        }

    private:
        virtual double distanceBetweenClusters(
            double smallerScale, const Cluster& smallerScaleCluster,
            double biggerScale, const Cluster& biggerScaleCluster) const = 0;
    };
}

#endif // FFTJET_ABSDISTANCECALCULATOR_HH_
