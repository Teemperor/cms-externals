//=========================================================================
// ClusteringTreeSparsifier.hh
//
// This class converts a dense clustering tree (with all levels) into
// a sparse clustering tree by selecting peaks with slow coordinate
// drift in the scale space (thus avoiding bifurcation points).
//
// I. Volobouev
// June 2010
//=========================================================================

#ifndef FFTJET_CLUSTERINGTREESPARSIFIER_HH_
#define FFTJET_CLUSTERINGTREESPARSIFIER_HH_

#include "fftjet/AbsClusteringTree.hh"
#include "fftjet/SparseClusteringTree.hh"

namespace fftjet {
    template<class Cluster, typename LevelInfo>
    class ClusteringTreeSparsifier
    {
    public:
        // Constructor arguments are as follows:
        //
        // maxLevelNumber -- Maximum level to write into the sparsified
        //                   tree. If this number is 0 or negative, this
        //                   number will be added to the maximum level
        //                   number in the input tree in order to obtain
        //                   the sparsified maximum level number. Note that
        //                   the trees we process will not necessarily
        //                   have the same number of levels (e.g., adaptive
        //                   trees). Also, the level with the highest
        //                   number may actually contain the complete set
        //                   of energy flow objects (however, this typically
        //                   does not affect the drift speed calculated for
        //                   the next-to-last level). This is why the default
        //                   value is -1 and not 0.
        //
        // filterMask     -- Some bits of this mask can be set to 0.
        //                   Corresponding reasons to write out the cluster
        //                   will be ignored (per enum at the beginning of the
        //                   public section of the "SparseClusteringTree"
        //                   class). Do not mask out either MIDRANGE_NODE or
        //                   SLOWEST_DRIFT, or you will be killing the very
        //                   cluster you want to look at!
        //
        // userScales     -- Scales for which all peaks should be written
        //                   out. The code will pick the nearest level.
        //                   Note, again, that these are not level numbers
        //                   because the number of levels in the input
        //                   trees is not necessarily fixed.
        //
        // nUserScales    -- Number of scales in the "userScales" array.
        //
        ClusteringTreeSparsifier(int maxLevelNumber = -1,
                                 unsigned filterMask = 0xffffffffU,
                                 const double* userScales = 0,
                                 unsigned nUserScales = 0);

        inline virtual ~ClusteringTreeSparsifier() {}

        virtual void sparsify(
            const AbsClusteringTree<Cluster,LevelInfo>& in,
            SparseClusteringTree<Cluster,LevelInfo>* output) const;

    protected:
        void transferScaleInfo(
            const AbsClusteringTree<Cluster,LevelInfo>& in,
            SparseClusteringTree<Cluster,LevelInfo>* output) const;

        virtual void descendTreeBranch(
            const AbsClusteringTree<Cluster,LevelInfo>& input,
            unsigned maxLevel, const TreeNodeId& topId,
            unsigned parentIndex,
            SparseClusteringTree<Cluster,LevelInfo>* output) const;

        std::vector<double> userScales_;
        mutable std::vector<unsigned> userLevels_;
        int maxLevelNumber_;
        unsigned filterMask_;
    };
}

#include "fftjet/ClusteringTreeSparsifier.icc"

#endif // FFTJET_CLUSTERINGTREESPARSIFIER_HH_
