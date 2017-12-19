//=========================================================================
// SimplePredicates.hh
//
// Some trivial predicates for use with AbsClusteringTree
// and other classes. Usually, "SimplePeakSelector" class will
// be more useful for the purpose of cluster filtering.
//
// I. Volobouev
// February 2009
//=========================================================================

#ifndef FFTJET_SIMPLEPREDICATES_HH_
#define FFTJET_SIMPLEPREDICATES_HH_

#include "fftjet/SimpleFunctors.hh"

namespace fftjet {
    template <bool result, typename Argument>
    struct Always : public Functor1<bool, Argument>
    {
        inline bool operator()(const Argument&) const {return result;}
    };
}

#endif // FFTJET_SIMPLEPREDICATES_HH_
