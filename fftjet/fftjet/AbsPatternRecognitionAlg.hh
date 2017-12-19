//=========================================================================
// AbsPatternRecognitionAlg.hh
//
// Interface class for multiscale pattern recognition algorithms.
//
// The utility of this class is not expected to be great -- it is likely
// that different pattern recognition strategies will need different
// concepts to describe their success or failure (the "PatternInfo" type
// on which this class is templated). In such a case an abstract base
// class no longer serves its main purpose very well. Anyway, it may be
// still useful to derive concrete pattern recognition strategies from
// this class. This standardizes data flow through the sequence of
// algorithm steps and makes the code somewhat easier to read and
// understand (and, of course, it makes sure that different algorithms
// using the same PatternInfo will be interchangeable).
//
// I. Volobouev
// May 2009
//=========================================================================

#ifndef FFTJET_ABSPATTERNRECOGNITIONALG_HH_
#define FFTJET_ABSPATTERNRECOGNITIONALG_HH_

#include "fftjet/AbsClusteringTree.hh"
#include "fftjet/Peak.hh"

namespace fftjet {
    template<typename PatternInfo>
    struct AbsPatternRecognitionAlg
    {
        virtual ~AbsPatternRecognitionAlg() {}

        // The "run" function should examine the given clustering
        // tree and fill out the collection of preclusters (peaks).
        // The "PatternInfo" class should provide the necessary
        // information about the "quality" (stability in the scale
        // space, etc) of the precluster collection found.
        //
        // The function should return a status word. 0 means everything
        // is fine. The meaning of other codes is up to concrete
        // implementations.
        //
        virtual int run(const AbsClusteringTree<Peak,long>& tree,
                        std::vector<Peak>* peaks, PatternInfo* info) = 0;
    };
}

#endif // FFTJET_ABSPATTERNRECOGNITIONALG_HH_
