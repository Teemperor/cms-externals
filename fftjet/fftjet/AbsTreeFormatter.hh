//=========================================================================
// AbsTreeFormatter.hh
//
// This class provides an interface for writing out AbsClusteringTree
// and SparseClusteringTree information for subsequent visualization.
//
// The derived classes must implement the "write" function. After this,
// formatters can be written out to an ostream, and it is not necessary
// to implement separate operator<< for them.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_ABSTREEFORMATTER_HH_
#define FFTJET_ABSTREEFORMATTER_HH_

#include <cassert>
#include <iostream>

namespace fftjet {
    template
    <
        typename Cluster,
        typename LevelInfo,
        template <typename, typename> class AbsTree
    >
    class AbsTreeFormatter
    {
    public:
        inline AbsTreeFormatter() : tree_(0), run_(0), n_(0) {}

        // The "runNum" and "evNum" arguments are some kind of
        // an identifier for the tree, typically run and event numbers
        inline AbsTreeFormatter(
            const AbsTree<Cluster,LevelInfo>& tree,
            const unsigned long runNum, const unsigned long evNum)
            : tree_(&tree), run_(runNum), n_(evNum) {}

        inline virtual ~AbsTreeFormatter() {}

        inline void setTree(const AbsTree<Cluster,LevelInfo>& tree,
                            unsigned long runNum, unsigned long evNum)
            {tree_ = &tree; run_ = runNum; n_ = evNum;}

        inline void writeTree(std::ostream& os) const
            {assert(tree_); this->write(*tree_, run_, n_, os);}

    private:
        const AbsTree<Cluster,LevelInfo>* tree_;
        unsigned long run_;
        unsigned long n_;

        virtual void write(const AbsTree<Cluster,LevelInfo>& tree,
                           unsigned long runNum,
                           unsigned long evNum,
                           std::ostream& os) const = 0;
    };
}

template
<
    typename Cluster,
    typename LevelInfo,
    template <typename, typename> class AbsTree
>
inline std::ostream& operator<<(
    std::ostream& strm,
    const fftjet::AbsTreeFormatter<Cluster,LevelInfo,AbsTree>& f)
{
    f.writeTree(strm);
    return strm;
}

#endif // FFTJET_ABSTREEFORMATTER_HH_
