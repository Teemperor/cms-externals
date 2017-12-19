//=========================================================================
// ExtrapolatedClusterIntegrator.hh
//
// This class implements a specific way to calculate the cluster
// energy from the information stored in the clustering tree.
// The energy integral is approximated by the value of
// scale*scale*clusterMagnitude
// extrapolated towards infinite scale (in the code, to zero of 1/scale).
// This method essentially assumes that cluster magnitudes are peak
// values of convoluted Et-like quantity.
//
// In order to use such integrals in a meaningful way, the user must
// multiply them by the correct conversion factor which, in general,
// depends on the kernel used to build the clustering tree and the degree
// of fitted polynomial.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef EXTRAPOLATEDCLUSTERINTEGRATOR_HH_
#define EXTRAPOLATEDCLUSTERINTEGRATOR_HH_

#include <vector>
#include "AbsClusterIntegrator.hh"

template<class Cluster, typename LevelInfo, class BooleanPredicate>
class ExtrapolatedClusterIntegrator : 
    public AbsClusterIntegrator<Cluster,LevelInfo>
{
public:
    // This class will not own the object to which "pred" points
    ExtrapolatedClusterIntegrator(BooleanPredicate* pred,
                         double maxRadius, double scaleInversionPower,
                         unsigned maxInterpolationDegree);
    virtual ~ExtrapolatedClusterIntegrator() {}

    // The ClusterIntegral status field is set as follows:
    //
    //   -1 : the original cluster does not satisfy the predicate
    //
    //    0 : only the original cluster satisfies the predicate/distance
    //        conditions and not any parents/daughters. In this case
    //        the extrapolation to 0 scale is not performed. The
    //        ClusterIntegral "value" field is set to
    //        scale*scale*clusterMagnitude for the original cluster.
    //
    //   >0 : this is the degree of polynomial fitted
    //
    virtual ClusterIntegral operator()(
        const fftjet::AbsClusteringTree<Cluster,LevelInfo>& tree,
        const fftjet::TreeNodeId& id);

private:
    ExtrapolatedClusterIntegrator();

    BooleanPredicate* pred_;
    const double maxDistance_;
    const double sInv_;
    const unsigned maxDeg_;

    std::vector<Cluster> clusters;
    std::vector<double> scales;
    std::vector<double> magnitudes;
    std::vector<double> coeffs;
};

#include "ExtrapolatedClusterIntegrator.icc"

#endif // EXTRAPOLATEDCLUSTERINTEGRATOR_HH_
