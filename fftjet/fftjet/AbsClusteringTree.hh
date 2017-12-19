//=========================================================================
// AbsClusteringTree.hh
//
// Interface class for the clustering tree
//
// I. Volobouev
// May 2008
//=========================================================================

#ifndef FFTJET_ABSCLUSTERINGTREE_HH_
#define FFTJET_ABSCLUSTERINGTREE_HH_

#include <vector>
#include <utility>
#include <climits>

namespace fftjet {
    // Forward declarations
    class StatAccumulator;

    // The tree nodes are identified by their scale space level
    // and their sequential number within the level
    typedef std::pair<unsigned,unsigned> TreeNodeId;

    template<class Cluster, typename LevelInfo>
    class AbsClusteringTree
    {
    public:
        typedef Cluster cluster_type;
        typedef LevelInfo info_type;
        typedef TreeNodeId NodeId;

        AbsClusteringTree();
        AbsClusteringTree(const AbsClusteringTree&);
        virtual ~AbsClusteringTree() {}

        // Assignment operator
        AbsClusteringTree& operator=(const AbsClusteringTree&);

        // This function inserts a new level and returns its number.
        // The concrete implementations of this class should normally
        // arrange levels in such a way that the levels with the
        // lagest scales have smallest numbers.
        virtual unsigned insert(double scale,
                                const std::vector<Cluster>& clusters,
                                const LevelInfo& levelInfo) = 0;

        // Change the level info
        virtual void setLevelInfo(unsigned level, const LevelInfo& info) = 0;

        // Remove the given level. An attempt to remove level 0
        // should be ignored. An attempt to remove level out of
        // range should result in a run-time error.
        virtual void remove(unsigned level) = 0;

        // Reset the whole tree leaving only the root node
        virtual void clear();

        // Number of scale levels in the tree. Level 0 corresponds
        // to the root level (it will always be there).
        virtual unsigned nLevels() const = 0;

        // Level corresponding to the given scale
        virtual unsigned getLevel(double scale) const = 0;

        // Check whether a certain scale value was used in building the tree
        virtual bool isScaleUsed(double scale) const;

        // Return the data which was provided during the level
        // creation. Calling this function on level 0 should
        // result in a run-time error.
        virtual void getLevelData(unsigned level, double* scale,
                                  std::vector<Cluster>* clustersToFill,
                                  LevelInfo* levelInfo) const = 0;

        // Return the list of nodes whose associated clusters satisfy
        // a certain predicate. The predicate must define the operator
        // "bool operator()(const Cluster&) const" which should return
        // "true" in case the node is to be included in the output.
        template<class BooleanPredicate>
        void getPassingNodes(unsigned level, const BooleanPredicate& pred,
                             std::vector<NodeId>* nodesToFill) const;

        // Number of clusters on the given level
        virtual unsigned nClusters(unsigned level) const = 0;

        // Total number of clusters in the tree. This count should
        // include the top-level cluster on the root level.
        virtual unsigned nClusters() const;

        // Number of clusters on the given level satisfying
        // a certain predicate. The predicate must define the operator
        // "bool operator()(const Cluster&) const" which should return
        // "true" in case the cluster is to be included in the count.
        template<class BooleanPredicate>
        unsigned clusterCount(unsigned level, const BooleanPredicate&) const;

        // Same thing for all levels. The answers are placed in a vector.
        // By convention, occupancy of level 0 will be set to 1.
        template<class BooleanPredicate>
        void occupancyInScaleSpace(const BooleanPredicate& pred,
                                   std::vector<unsigned>* occupancy) const;

        // Statistics for clusters on the given level satisfying
        // a certain predicate. The predicate must define the operator
        // "bool operator()(const Cluster&) const".
        // The "PropertySelector" functor must define the operator
        // "double operator()(const Cluster&) const".
        template<typename PropertySelector, class BooleanPredicate>
        void clusterStats(unsigned level, const PropertySelector&,
                          const BooleanPredicate&, StatAccumulator*) const;

        // The following function finds minimum and maximum tree levels
        // for which the number of clusters satisfying a certain predicate
        // equals the desired count. If no levels have the requested number
        // of clusters, zeros are returned for both min and max levels.
        // This function returns "true" if all the levels between
        // minimum and maximum have the same number of clusters
        // (equal to requested), otherwise the function returns "false".
        template<class BooleanPredicate>
        bool clusterCountLevels(unsigned desiredCount, const BooleanPredicate&,
                                unsigned* minLevel, unsigned* maxLevel) const;

        // The following function returns the number of clusters N for which
        // stability of the configuration is the highest. Here, "stability"
        // is defined as pow(N,alpha)*log(scaleMax/scaleMin), where "scaleMax"
        // is the maximum scale with exactly N clusters satisfying the
        // predicate, "scaleMin" is the minimum scale with N such clusters,
        // and all tree levels with scales between scaleMax and scaleMin
        // have the same number of clusters, N.
        //
        // The input variables "minStartingLevel" and "maxStartingLevel"
        // can be used to limit the searched range of scales. The scale
        // which corresponds to "maxStartingLevel" is used in the search,
        // and by default the full range of tree scales is used. The levels
        // which correspond to the bounds of the most stable configuration
        // are returned in *minLevel and *maxLevel.
        //
        // If there are no stable configurations (that is, the number of
        // clusters is different on each level) the function returns 0
        // and sets *minLevel and *maxLevel to 0.
        //
        // Note that this function is going to be useful only when the
        // clustering kernel (as well as the input predicate) are
        // "reasonable". In particular, it is highly desirable that
        // the number of clusters decreases monotonously as the scale
        // increases.
        //
        template<class BooleanPredicate>
        unsigned stableClusterCount(const BooleanPredicate& pred,
                                    unsigned* minLevel, unsigned* maxLevel,
                                    double alpha = 0.0,
                                    unsigned minStartingLevel=0,
                                    unsigned maxStartingLevel=UINT_MAX) const;

        // Scale corresponding to the given level
        virtual double getScale(unsigned level) const = 0;

        // Minimum and maximum scales. The default implementations
        // of these methods assume that level 1 has the largest scale.
        // 0 is returned in case the tree is not populated.
        virtual double minScale() const;
        virtual double maxScale() const;

        // Scales for all levels except level 0
        virtual void getAllScales(std::vector<double>* scales,
                                  bool increasingOrder=false) const;

        // Level-wide information (e.g., unclustered energy) provided
        // during the level creation. Note that the reference can be
        // invalidated by insert/remove/clear operations.
        virtual const LevelInfo& getLevelInfo(unsigned level) const = 0;

        // Cluster data for the given tree node. The "unchecked"
        // function should work faster by assuming that the cluster
        // with the given id exists.
        virtual const Cluster& getCluster(const NodeId& id) const=0;
        virtual const Cluster& uncheckedCluster(const NodeId& id) const=0;

        // "Distance" between the two clusters with given ids.
        // The exact meaning of this distance is up to derived classes.
        // It is very important, however, that it is symmetric and that
        // it satisfies the triangle inequality.
        virtual double clusterDistance(const NodeId& id1,
                                       const NodeId& id2) const = 0;

        // Same thing, but should work faster by assuming that
        // the clusters with the given ids exist
        virtual double uncheckedClusterDistance(const NodeId& id1,
                                                const NodeId& id2) const=0;

        // The extent of the node descendants (think balltree). The size
        // of the whole tree is undefined because the toplevel cluster
        // does not have a meaningful location.
        virtual double nodeRadius(const NodeId& id) const = 0;
        virtual double uncheckedNodeRadius(const NodeId& id) const = 0;

        // "Distance" to the parent for the given cluster.
        // The smaller the distance, the better the match between
        // the daughter and its parent. The derived classes are
        // encouraged to override the default implementation because
        // distances to parents are likely to be calculated during
        // the tree construction anyway.
        virtual double distanceToParent(const NodeId& id) const;

        // The following function returns the id of the closest cluster
        // at the same tree level which satisfies the predicate, together
        // with the distance to that cluster. Bad node id is returned
        // in case there are no neighbor clusters satisfying the predicate.
        template<class BooleanPredicate>
        NodeId closestNeighbor(const NodeId& id,
                               const BooleanPredicate& pred,
                               double* distance) const;

        // The following function returns the id of the closest cluster
        // at the same or lower tree level which satisfies the predicate
        // and which is not one of the descendants of this cluster.
        // Bad node id is returned in case there are no such clusters.
        //
        // The use of this function should be limited to predicates
        // which can be satisfied by a resonably large fraction of
        // clusters, otherwise the code may end up searching almost the
        // entire tree at and below the level of the node with given id.
        template<class BooleanPredicate>
        NodeId closestNonDescendant(const NodeId& id,
                                    const BooleanPredicate& pred,
                                    double* distance) const;

        // The following function returns the node radius using
        // only the clusters satisfying the predicate
        template<class BooleanPredicate>
        double conditionalNodeRadius(const NodeId& id,
                                     const BooleanPredicate& pred) const;

        // Tree navigation. Note that insert and remove operations
        // can invalidate all node ids and numbers of daughters
        // obtained earlier.
        //
        // The daughters must be arranged in the order of increasing
        // distance to parent. badId should be returned as a parent
        // of the root node.
        //
        virtual NodeId parentId(const NodeId& id) const = 0;
        virtual unsigned nDaughters(const NodeId& id) const = 0;
        virtual NodeId daughterId(const NodeId& id,
                                  unsigned idau) const = 0;

        // The following function should be implemented if the tree
        // itself knows how to grow optimally. It should return
        // the scale which would improve the level matching in the best
        // possible place (plug the biggest hole). The function should
        // return 0.0 when there is nothing more to gain from adding
        // intermediate levels.
        //
        // Note that the tree should be initialized with at least two
        // data levels (with the largest and the smallest scales) for
        // the default implementation of this function to work.
        //
        // The "minRatioLog" parameter prevents the code from trying
        // to request extremely close scale values in pathological
        // configurations. If log(parentLevelScale/daughterLevelScale)
        // is less than "minRatioLog" then the function should not
        // request a new scale between these existing scales.
        virtual double nextBestScale(double minRatioLog) const;

        // The default implementation of the nextBestScale() function
        // uses the following function to figure out how well a level
        // is matched to its parent level. The function should return 0
        // for a good match. If you want to improve on the default tree
        // growing algorithm but do not want to create a fancy multi-step
        // optimization strategy then it should be sufficient to override
        // just this function and rely on the default nextBestScale().
        virtual double levelMatchDistance(unsigned daughterLevel) const;

        // Various useful utilities
        virtual NodeId closestDaughter(const NodeId& id) const;
        bool operator==(const AbsClusteringTree& r) const;
        bool operator!=(const AbsClusteringTree& r) const;

        // The following function says whether node with id
        // id1 is an ancestor of the node with id id2
        virtual bool isAncestor(const NodeId& id1,
                                const NodeId& id2) const;

        // The following function fills a continuous sequence of
        // clusters (in the order of increasing level number) which are
        // all related to each other through the parent/closest daughter
        // relationship, satisfy the given predicate, and are all
        // within the given distance from the given initial cluster.
        // The function returns the position of the initial cluster
        // inside the sequence.
        //
        // -1 is returned (and the sequence is cleared) if the initial
        // cluster does not satisfy the predicate
        //
        template<class BooleanPredicate>
        int heritageLine(const NodeId& initialId,
                         const BooleanPredicate& pred,
                         double maxDistance,
                         std::vector<Cluster>* clusterSequence) const;

        // Number of clusters on the given level with the
        // number of daughters equal to or above the given
        // minimum number
        virtual unsigned nParentsWithDaus(unsigned level,
                                          unsigned nMinDau) const;

        // Number of clusters connected directly to the root node.
        // All clusters are like that at the level with the
        // largest scale.
        virtual unsigned nRootDaughters(unsigned level) const;

        // The following function calculates and remembers
        // the drift speed for all clusters. Call this function
        // as needed after completing the tree construction.
        virtual void calculateDriftSpeeds();

        // The following function calculates and remembers
        // the magnitude change speed for all clusters. Call this
        // function as needed after completing the tree construction.
        virtual void calculateMagnitudeSpeeds();

        // The following function calculates and remembers
        // the lifetime for all clusters. Call this function
        // as needed after completing the tree construction.
        virtual void calculateLifetimes();

        // The following function calculates and remembers the cluster
        // distance to the nearest neighbor cluster at the same scale.
        // Call this function as needed after completing the tree
        // construction. The "noNeighborDistance" argument specifies
        // what number to use when the cluster has no neighbors.
        virtual void calculateNearestNeighbors(double noNeighborDistance=-2.0);

        // The following functions transfer the information from the tree
        // nodes to the clusters. The "updateClusterRadiusInfo" method uses
        // "nodeRadius" method to evaluate the radii, and the
        // "updateClusterSeparationInfo" method uses "closestNonDescendant"
        // with an "all pass" predicate. The "failDistance" argument tells
        // us how to set the result if the quantity can not be calculated
        // (for example, there is only one cluster in the whole tree).
        // The results are not produced for the top level and for the level
        // with the smallest scale (for these levels they are meaningless).
        virtual void updateClusterRadiusInfo();
        virtual void updateClusterSeparationInfo(double failDistance=-2.0);

        // A convenience function which invokes the default sequence
        // of the tree post-processing steps. If you do not like
        // this sequence, just call the functions you need directly.
        virtual inline void postProcess()
        {
            calculateDriftSpeeds();
            calculateMagnitudeSpeeds();
            calculateLifetimes();
            calculateNearestNeighbors();
        }

        // Bad node id. Should be returned when the user requests
        // the parent id of the root node or the closest daughter
        // of a daughterless node.
        const NodeId badId;

    private:
        virtual unsigned traceLifetime(const TreeNodeId& id, unsigned topLev);

        template<class BooleanPredicate>
        void closestDauOrSelf(const TreeNodeId& target,
                              const TreeNodeId& id,
                              const BooleanPredicate& pred,
                              double* bestDistanceSoFar,
                              TreeNodeId* bestIdSoFar) const;

        template<class BooleanPredicate>
        void farthestDauOrSelf(const TreeNodeId& target,
                               const TreeNodeId& id,
                               const BooleanPredicate& pred,
                               double* bestDistanceSoFar,
                               TreeNodeId* bestIdSoFar) const;

        mutable std::vector<std::pair<double, unsigned> > neigbors_;
    };
}

#include "fftjet/AbsClusteringTree.icc"

#endif // FFTJET_ABSCLUSTERINGTREE_HH_
