//=========================================================================
// peakEtLifetime.hh
//
// Cluster lifetime weighted by the local cluster Et fraction.
// Intended for use with sparse clustering trees.
//
// I. Volobouev
// May 2013
//=========================================================================

#ifndef FFTJET_PEAKETLIFETIME_HH_
#define FFTJET_PEAKETLIFETIME_HH_

namespace fftjet {
    template<class SparseTree>
    double peakEtSplitTime(const SparseTree& tree,
                           typename SparseTree::NodeId id,
                           double minScale);

    template<class SparseTree>
    double peakEtMergeTime(const SparseTree& tree,
                           typename SparseTree::NodeId id,
                           double maxScale);

    template<class SparseTree>
    void updateSplitMergeTimes(SparseTree& tree, double minScale, double maxScale);
}

#include "fftjet/peakEtLifetime.icc"

#endif // FFTJET_PEAKETLIFETIME_HH_
