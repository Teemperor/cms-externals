//=========================================================================
// ProximityClusteringTree.hh
//
// Clustering tree based on simple proximity matching in the eta-phi space
// from the clusters at lower scale to clusters at higher scale
//
// I. Volobouev
// April 2008
//=========================================================================

#ifndef FFTJET_PROXIMITYCLUSTERINGTREE_HH_
#define FFTJET_PROXIMITYCLUSTERINGTREE_HH_

#include "fftjet/SmallVector.hh"
#include "fftjet/AbsClusteringTree.hh"
#include "fftjet/AbsDistanceCalculator.hh"

namespace fftjet {
    template<class Cluster, typename LevelInfo>
    class ProximityClusteringTree : public AbsClusteringTree<Cluster,LevelInfo>
    {
    public:
        explicit ProximityClusteringTree(
            const AbsDistanceCalculator<Cluster>* calc,
            bool treeOwnsCalc=false);
        virtual ~ProximityClusteringTree();

        // Copy constructor
        ProximityClusteringTree(const ProximityClusteringTree&);

        // Assignment operator
        ProximityClusteringTree& operator=(const ProximityClusteringTree&);

        // Inserts a new level and returns its number
        unsigned insert(double scale,
                        const std::vector<Cluster>& clusters,
                        const LevelInfo& levelInfo);

        // Change the level info
        void setLevelInfo(unsigned level, const LevelInfo& info);

        // Remove the given level. An attempt to remove level 0
        // should be ignored. An attempt to remove level out of
        // range should result in a run-time error.
        void remove(unsigned level);

        // Reset the whole tree leaving only the root node
        void clear();

        // Number of scale levels in the tree. Level 0 corresponds
        // to the root level (it will always be there).
        unsigned nLevels() const;

        // Level corresponding to the given scale
        unsigned getLevel(double scale) const;

        // Return the data which was provided during the level
        // creation. Calling this function on level 0 should
        // result in a run-time error.
        void getLevelData(unsigned level, double* scale,
                          std::vector<Cluster>* clustersToFill,
                          LevelInfo* levelInfo) const;

        // Number of clusters on the given level
        unsigned nClusters(unsigned level) const;

        // Scale corresponding to the given level
        double getScale(unsigned level) const;

        // Level-wide information (e.g., unclustered energy) provided
        // during the level creation
        const LevelInfo& getLevelInfo(unsigned level) const;

        // Cluster data for the given tree node. The "unchecked"
        // function should work faster by assuming that the cluster
        // with the given id exists.
        const Cluster& getCluster(const TreeNodeId&) const;
        const Cluster& uncheckedCluster(const TreeNodeId& id) const;

        // "Distance" between the two clusters with given ids
        double clusterDistance(const TreeNodeId& id1,
                               const TreeNodeId& id2) const;

        // Same thing, but should work faster by assuming that
        // the clusters with the given ids exist
        double uncheckedClusterDistance(const TreeNodeId& id1,
                                        const TreeNodeId& id2) const;

        // The extent of the node descendants (think balltree)
        double nodeRadius(const TreeNodeId& id) const;
        double uncheckedNodeRadius(const TreeNodeId& id) const;

        // "Distance" to parent for the given cluster.
        // The smaller the distance, the better the match
        // between the daughter and its parent.
        double distanceToParent(const TreeNodeId& id) const;

        // Tree navigation. Note that insert and remove operations
        // can invalidate all node ids obtained earlier.
        TreeNodeId parentId(const TreeNodeId& id) const;
        unsigned nDaughters(const TreeNodeId& id) const;
        TreeNodeId daughterId(const TreeNodeId& id, unsigned idau) const;

    private:
        ProximityClusteringTree();

        struct TreeNode;
        typedef TreeNode* TreeNodePtr;

        struct TreeNode
        {
            TreeNode(const TreeNodeId& id,
                     const Cluster& jet);

            // Balltree construction
            void updateRadius(const ProximityClusteringTree& tree);

            // Add daughter in the right order
            void addDaughter(TreeNodePtr dau);

            // Sorting will be according to the cluster magnitude
            inline bool operator<(const TreeNode& r) const
                {return jet < r.jet;}
            inline bool operator>(const TreeNode& r) const
                {return jet > r.jet;}

            // node payload
            Cluster jet;

            // tree navigation
            TreeNodeId id;
            TreeNodePtr parent;
            double distanceToParent;
            double radius;
            SmallVector<TreeNodePtr,3U> daus;

        private:
            TreeNode();

            // Largest distance from this node and from all
            // of its subnodes to some other node
            double largestToNode(const ProximityClusteringTree& tree,
                                 const TreeNode* other) const;
        };

        struct TreeLevel
        {
            TreeLevel(const TreeLevel&);
            TreeLevel(double scale, const LevelInfo& unclus, unsigned njets);

            const double scale;
            LevelInfo unclustered;
            std::vector<TreeNode> nodes;

        private:
            TreeLevel();
            TreeLevel& operator=(const TreeLevel&);
        };

        void renumberLevels(unsigned startingLevel);
        void connectToParent(TreeLevel* dau, TreeLevel* parent);
        void disconnectFromParent(TreeLevel* dau);
        TreeNode* getNodeById(const TreeNodeId& id) const;
        void reestablishParenthood();
        void calculateRadii();

        std::vector<TreeLevel*> levels;
        const AbsDistanceCalculator<Cluster>* distance;
        bool treeOwnsCalc;
        mutable bool recalculateRadii;
    };
}

#include "fftjet/ProximityClusteringTree.icc"

#endif // FFTJET_PROXIMITYCLUSTERINGTREE_HH_
