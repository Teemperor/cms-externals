//=========================================================================
// SparseClusteringTree.hh
//
// Hierarchical clustering tree without complete level-by-level structure.
//
// Do not create this object every event -- this will lead to many
// unnecessary memory allocations and deallocations. Instead, create it
// at the start of your program and reuse.
//
// I. Volobouev
// June 2010
//=========================================================================

#ifndef FFTJET_SPARSECLUSTERINGTREE_HH_
#define FFTJET_SPARSECLUSTERINGTREE_HH_

#include <vector>
#include <utility>
#include <climits>

#include "fftjet/SmallVector.hh"

namespace fftjet {
    // Same template parameters as in AbsClusteringTree. The second
    // parameter is currently unused, but included here anyway for
    // uniformity with AbsClusteringTree (this is useful for generic
    // code in the I/O classes).
    template<class Cluster, typename LevelInfo=long>
    class SparseClusteringTree
    {
    public:
        typedef Cluster cluster_type;
        typedef LevelInfo info_type;
        typedef unsigned NodeId;

        // Node mask bits
        enum {
            MIN_LEVEL     = 1U,    // level with the largest scale
            MAX_LEVEL     = 2U,    // level with the smallest scale
            USER_LEVEL    = 4U,    // user-specified scale to preserve
            SPLIT_NODE    = 8U,    // more than one daughter
            TERMINAL_NODE = 16U,   // no daughters
            MIDRANGE_NODE = 32U,   // node in the middle of the branch
            SLOWEST_DRIFT = 64U,   // node with the slowest drift speed
            RESERVED_LOCK = 128U,  // reserved for internal use
            RESERVED_CNT  = 256U,  // reserved for internal use
            RESERVED_TAG  = 512U,  // reserved for internal use
            RESERVED_PASS = 1024U  // reserved for internal use
        };

        class Node
        {
        public:
            Node(const Cluster& jet, unsigned level, unsigned mask);

            inline const Cluster& getCluster() const {return jet_;}
            inline unsigned originalLevel() const {return originalLevel_;}
            inline unsigned mask() const {return nodeMask_;}
            inline unsigned parent() const {return parent_;}
            inline const SmallVector<unsigned,3U>& daus() const {return daus_;}
            inline unsigned nDaus() const {return daus_.size();}

            // Non-const access to the payload
            inline Cluster& mutableCluster() {return jet_;}

            bool operator==(const Node& r) const;
            inline bool operator!=(const Node& r) const {return !(*this == r);}

        private:
            Node();

            // Node payload
            Cluster jet_;

            // Original level number helps to speed-up various algorithms
            unsigned originalLevel_;

            // The reason why this node was included in this tree.
            // Should be constructed by bitwise "OR" using the enums
            // declared at the beginning of the public section.
            // Also used by several SparseClusteringTree algorithms.
            mutable unsigned nodeMask_;

            // Tree navigation. It is assumed that the node with number 0
            // is the top node. It will have a dummy payload.
            unsigned parent_;
            SmallVector<unsigned,3U> daus_;

            friend class SparseClusteringTree;
        };

        // Default constructor builds an empty tree with a top node
        SparseClusteringTree();
        SparseClusteringTree(const SparseClusteringTree&);
        virtual ~SparseClusteringTree() {}

        SparseClusteringTree& operator=(const SparseClusteringTree&);

        // The following function returns the node index which
        // can be used for quick lookup of node contents. Note
        // that the parent index must be correct. Direct daughters
        // of top node must have parent index 0.
        //
        // It is assumed that this tree is produced from another,
        // larger tree via some kind of a pruning process. We do not
        // anticipate the need to shrink this tree even further,
        // and this is why there is no "removeNode" method. If this
        // becomes really necessary, this tree should be cleared
        // completely and recreated (or use another, more appropriate
        // data structure).
        //
        NodeId addNode(const Node& node, NodeId parent);

        // Remove all nodes except the top one and forget all scales
        void clear();

        // If you know in advance how many nodes/scales you are going
        // to make, use the following functions
        inline void reserveNodes(const unsigned n) {nodes_.reserve(n);}
        inline void reserveScales(const unsigned n) {scales_.reserve(n);}

        // Add a scale. For correct operation, all scales should be added,
        // starting with the level 0.
        void addScale(double scale);

        // Sort the contents of the tree by the product of cluster
        // magnitude and scale squared. If necessary, this sorting
        // should be performed after all nodes have been added.
        void sortNodes();
        inline bool isSorted() const {return sorted_;}

        // Functions for retrieving the maximum level number stored
        inline unsigned maxStoredLevel() const {return maxLevel_;}

        // The following functions have the same meaning as in
        // AbsClusteringTree. It is assumed that the scale info
        // was correctly transferred from AbsClusteringTree into
        // this object.
        inline unsigned nLevels() const {return scales_.size();}
        inline double getScale(const unsigned level) const
            {return scales_.at(level);}
        unsigned getLevel(double scale) const;

        // Minimum and maximum scales.
        // 0 is returned in case the scales are not populated.
        double minScale() const;
        double maxScale() const;

        // Scales for all levels except level 0
        void getAllScales(std::vector<double>* scales,
                          bool increasingOrder=false) const;

        // How many nodes are there (including the top node)
        inline unsigned size() const {return nodes_.size();}

        // Look up the node. Index 0 corresponds to the top node.
        inline const Node& getNode(const NodeId index) const
            {return nodes_.at(index);}
        inline const Node& uncheckedNode(const NodeId index) const
            {return nodes_[index];}

        // Non-const access to the nodes
        inline Node& mutableNode(const NodeId index)
            {return nodes_.at(index);}
        inline Node& uncheckedMutableNode(const NodeId index)
            {return nodes_[index];}

        // Some useful utilities
        bool operator==(const SparseClusteringTree& r) const;
        inline bool operator!=(const SparseClusteringTree& r) const
            {return !(*this == r);}

        // Number of clusters on the given level. Note that
        // this operation is a lot slower than the similar method
        // of the "AbsClusteringTree" class: the code here actually
        // has to drill down the branch structure and count branches
        // which intersect the given level.
        unsigned nClusters(unsigned level) const;

        // Number of (slowest drift) clusters on the branches which
        // intersect the given level and which satisfy a certain
        // predicate. The predicate must define the operator
        // "bool operator()(const Cluster&) const" which should return
        // "true" in case the cluster is to be included in the count.
        template<class BooleanPredicate>
        unsigned clusterCount(unsigned level, const BooleanPredicate&) const;

        // Fill a vector of (slowest drift) nodes on the branches
        // which intersect the given level
        void getLevelNodes(unsigned level, std::vector<NodeId>* nodes) const;

        // Fill a vector of nodes whose associated clusters satisfy
        // a certain predicate. The predicate must define the operator
        // "bool operator()(const Cluster&) const" which should return
        // "true" in case the node is to be included in the output.
        template<class BooleanPredicate>
        void getPassingNodes(unsigned level, const BooleanPredicate& pred,
                             std::vector<NodeId>* nodes) const;

        // The following method returns the number of tree branches
        // which satisfy a certain predicate as a function of
        // level number. The slowest drift speed node on each branch
        // will be tested by the predicate to determine if the branch
        // satisfies the predicate. The occupancy vector will have
        // maxStoredLevel()+1 elements. By convention, occupancy of
        // level 0 will be set to 1.
        template<class BooleanPredicate>
        void occupancyInScaleSpace(const BooleanPredicate& pred,
                                   std::vector<unsigned>* occupancy) const;

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

        // The following function has the same meaning as the function
        // with the same name in the "AbsClusteringTree" class
        template<class BooleanPredicate>
        unsigned stableClusterCount(const BooleanPredicate& pred,
                                    unsigned* minLevel, unsigned* maxLevel,
                                    double alpha = 0.0,
                                    unsigned minStartingLevel=0,
                                    unsigned maxStartingLevel=UINT_MAX) const;

        // The following function determines the largest number of
        // (slowest drift) clusters which can be constructed in such
        // a way that the value of cluster magnitude times the scale
        // squared, MagS2, is above a certain threshold, and the
        // clusters satisfy the given predicate. The threshold values
        // corresponding to the given cluster count are returned in
        // the "thresholds" vector. For example, thresholds[N] == t
        // means that we can not construct more than N clusters with
        // MagS2 >= t, no matter what the branch selection is. The
        // parameter "maxClusterN" tells the function to stop when the
        // result for thresholds[maxClusterN] is obtained. The size
        // of the output vector will be maxClusterN + 1 or less.
        // 
        // The tree must be sorted before this method is called.
        //
        template<class BooleanPredicate>
        void nClustersAboveMagS2Threshold(
            const BooleanPredicate& pred,
            std::vector<double>* thresholds,
            unsigned maxClusterN = UINT_MAX-1U) const;

        // The following function returns the set of node numbers
        // optimized in such a manner that, when clusters are arranged
        // in the decreasing order using the MagS2 variable, the cluster
        // with number "nOpt" has the largest possible value of MagS2.
        // If the tree has nodes satisfying the predicate beyond nOpt,
        // they will be taken from the requested level. Note that the
        // levels (and, therefore, the resolution scales) for the first
        // nOpt nodes will be chosen automatically.
        //
        // The function can also fill the "thresholds" vector, in a manner
        // identical to "nClustersAboveMagS2Threshold" method. The size of
        // the filled vector will normally be nOpt + 2.
        //
        // The tree must be sorted before this method is called.
        //
        template<class BooleanPredicate>
        void getMagS2OptimalNodes(const BooleanPredicate& pred,
                                  unsigned nOpt, unsigned level,
                                  std::vector<NodeId>* nodes,
                                  std::vector<double>* thresholds = 0) const;

        // The following function says whether node with id
        // id1 is an ancestor of the node with id id2
        bool isAncestor(NodeId id1, NodeId id2) const;

        // Bad id will be returned for the parent index of the top node
        const NodeId badId;

    protected:
        std::vector<Node> nodes_;
        std::vector<double> scales_;
        unsigned maxLevel_;

    private:
        bool sorted_;

        // Branches are identified by the index of the parent node
        // and by daughter number in the parent. The parent node itself
        // does not belong to the branch.
        void getLevelNodesOnBranch(unsigned level,
                                   std::vector<unsigned>* nodes,
                                   unsigned parentIndex,
                                   unsigned dauNumber,
                                   unsigned vetoBits = 0U) const;
        unsigned countLevelNodesOnBranch(unsigned level,
                                         unsigned parentIndex,
                                         unsigned dauNumber) const;
        void scaledMagnitudeOnBranch(unsigned parentIndex,
                                     unsigned dauNumber,
                                     double currentValue);
        void lockDaus(unsigned parentIndex, unsigned dauNumber,
                      bool lock) const;
        unsigned relockDausAboveId(unsigned parentIndex, unsigned dauNumber,
                                   NodeId id) const;

        // Helper methods for magS2-related calculations
        template<class BooleanPredicate>
        void magS2Init(const BooleanPredicate& pred,
                       std::vector<double>* thresholds,
                       unsigned maxClusterN) const;

        // In the function below, nodeNumber is the in/out parameter.
        // The function returns the actual number of fills made (or
        // would be made if thresholds is NULL). The first "dummy"
        // fill is not counted.
        unsigned magS2Cycle(std::vector<double>* thresholds,
                            unsigned* nodeNumber,
                            unsigned maxFills) const;

        // Check whether the daughter with number "dau" of the given parent
        // and/or at least one of the daughter descendants satisfy the
        // following conditions:
        //  1) "slowest drift speed"
        //  2) not an ancestor of the node with id "id"
        //  3) id number is smaller than "id"
        bool slowDausAboveId(unsigned parent, unsigned dau, NodeId id) const;

        template<class BooleanPredicate>
        void branchOccupancy(const BooleanPredicate& pred,
                             std::vector<unsigned>* occupancy,
                             unsigned parentIndex, unsigned dauNumber) const;

        // These are essentially memory buffers for various methods
        std::vector<Node> swapped_;
        mutable std::vector<unsigned> inverseMap_;

        // The following buffer should not be modified anywhere outside
        // the "sortNodes" (and "scaledMagnitudeOnBranch" called by the
        // "sortNodes") method. This is because its contents are used
        // by other functions when the tree is in the sorted state.
        std::vector<std::pair<double,unsigned> > moveMap_;
    };
}

#include "fftjet/SparseClusteringTree.icc"

#endif // FFTJET_SPARSECLUSTERINGTREE_HH_
