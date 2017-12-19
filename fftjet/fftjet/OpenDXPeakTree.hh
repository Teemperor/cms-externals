//=========================================================================
// OpenDXPeakTree.hh
//
// This class provides a way to write out AbsClusteringTree<Peak, ...>
// and SparseClusteringTree<Peak, ...> information in a form suitable
// for visualization by OpenDX.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_OPENDXPEAKTREE_HH_
#define FFTJET_OPENDXPEAKTREE_HH_

#include <cmath>

#include "fftjet/AbsOpenDXTreeFormatter.hh"
#include "fftjet/Peak.hh"
#include "fftjet/SimpleFunctors.hh"

namespace fftjet {
    template
    <
        typename LevelInfo,
        template <typename, typename> class Tree
    >
    class OpenDXPeakTree : public AbsOpenDXTreeFormatter<Peak,LevelInfo,Tree>
    {
    public:
        inline OpenDXPeakTree(
            const Functor1<double,Peak>* glyphSizeFunctor,
            const Functor1<double,Peak>* glyphColorFunctor,
            const Tree<Peak,LevelInfo>& tree,
            const unsigned long runNum,
            const unsigned long evNum,
            const double etaRange = 0.0,
            const double scaleRange = 0.0)
            : AbsOpenDXTreeFormatter<Peak,LevelInfo,Tree>(
                tree, runNum, evNum, etaRange, scaleRange),
              glyphSizeFunctor_(glyphSizeFunctor),
              glyphColorFunctor_(glyphColorFunctor) {}

        inline OpenDXPeakTree(
            const Functor1<double,Peak>* glyphSizeFunctor,
            const Functor1<double,Peak>* glyphColorFunctor,
            const double etaRange = 0.0, const double scaleRange = 0.0)
            : AbsOpenDXTreeFormatter<Peak,LevelInfo,Tree>(etaRange,scaleRange),
              glyphSizeFunctor_(glyphSizeFunctor),
              glyphColorFunctor_(glyphColorFunctor) {}

        inline virtual ~OpenDXPeakTree() {}

    private:
        OpenDXPeakTree();

        const Functor1<double,Peak>* const glyphSizeFunctor_;
        const Functor1<double,Peak>* const glyphColorFunctor_;

        inline virtual DXGlyphInfo buildGlyph(const Peak& peak,
                                              const int parent) const
        {
            return DXGlyphInfo(peak.eta(), peak.phi(), log(peak.scale()),
                               (*glyphSizeFunctor_)(peak),
                               (*glyphColorFunctor_)(peak),
                               parent);
        }
    };
}

#endif // FFTJET_OPENDXPEAKTREE_HH_
