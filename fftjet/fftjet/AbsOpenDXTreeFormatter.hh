//=========================================================================
// AbsOpenDXTreeFormatter.hh
//
// This class provides a way to write out AbsClusteringTree and
// SparseClusteringTree information in a form suitable for visualization
// by OpenDX ( www.opendx.org ).
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_ABSOPENDXTREEFORMATTER_HH_
#define FFTJET_ABSOPENDXTREEFORMATTER_HH_

#include "fftjet/AbsTreeFormatter.hh"
#include "fftjet/DXGlyphInfo.hh"

namespace fftjet {
    // Forward declaration
    template
    <
        typename Cluster,
        typename LevelInfo,
        template <typename, typename> class AbsTree
    >
    class AbsOpenDXTreeFormatter;

    // We need to specialize AbsOpenDXTreeFormatter class slightly
    // so that it works a little bit differently for AbsClusteringTree
    // and for SparseClusteringTree. Unfortunately, C++ does not allow
    // us to partially specialize member (or non-member) functions.
    // This is why the following work-around class is needed.
    namespace Private {
        template
        <
            typename Cluster,
            typename LevelInfo,
            template <typename, typename> class AbsTree
        >
        struct AbsOpenDXTreeFormatterWriter
        {
            void operator()(
                const AbsOpenDXTreeFormatter<Cluster,LevelInfo,AbsTree>&,
                const AbsTree<Cluster,LevelInfo>& tree,
                unsigned long runNum, unsigned long evNum, std::ostream& os);
        };
    }

    // The following class performs all the work necessary
    // to construct the OpenDX "field" object except building
    // the glyphs themselves. The task of building glyphs
    // is delegated to derived classes.
    template
    <
        typename Cluster,
        typename LevelInfo,
        template <typename, typename> class AbsTree
    >
    class AbsOpenDXTreeFormatter :
        public AbsTreeFormatter<Cluster,LevelInfo,AbsTree>
    {
    public:
        //
        // The "etaRange" parameter can be useful to specify whether
        // you would like the trees to be displayed in the eta
        // range at least from -etaRange to etaRange. If it is
        // specified as 0.0 (default) then the range of the tree
        // itself will be used as the range of the plot. This may be
        // inconvenient if you plan on looking at many plots.
        //
        // The "scaleRange" parameter can be used to change the range
        // of scales produced by the code. If this parameter
        // is left at its default value, 0.0, then the code will
        // simply keep the calculated scale values unchanged.
        // If this parameter is not 0.0 then these values will be
        // remapped linearly to the range from 0.0 to "scaleRange".
        // This may be useful to get a better picture in OpenDX
        // (however, the plotted scales will be somewhat arbitrary).
        // If you want to do this remapping, a good value of
        // "scaleRange" is around 6 or so -- then the scale range
        // in the visualization will be similar to the phi range.
        //
        inline AbsOpenDXTreeFormatter(
            const AbsTree<Cluster,LevelInfo>& tree,
            const unsigned long runNum,
            const unsigned long evNum,
            const double etaRange = 0.0,
            const double scaleRange = 0.0)
            : AbsTreeFormatter<Cluster,LevelInfo,AbsTree>(tree, runNum, evNum),
              etaRange_(fabs(etaRange)),
              scaleRange_(scaleRange)
        {
        }

        inline explicit AbsOpenDXTreeFormatter(
            const double etaRange = 0.0, const double scaleRange = 0.0)
            : AbsTreeFormatter<Cluster,LevelInfo,AbsTree>(),
              etaRange_(fabs(etaRange)), scaleRange_(scaleRange) {}

        inline virtual ~AbsOpenDXTreeFormatter() {}

        inline double getEtaRange() const {return etaRange_;}
        inline double getScaleRange() const {return scaleRange_;}

    private:
        friend struct Private::AbsOpenDXTreeFormatterWriter<
            Cluster,LevelInfo,AbsTree>;

        AbsOpenDXTreeFormatter();

        const double etaRange_;
        const double scaleRange_;

        void write(const AbsTree<Cluster,LevelInfo>& tree,
                   unsigned long runNum, unsigned long evNum,
                   std::ostream& os) const;

        virtual DXGlyphInfo buildGlyph(const Cluster& clus,int parent) const=0;
    };
}

#include "fftjet/AbsOpenDXTreeFormatter.icc"

#endif // FFTJET_ABSOPENDXTREEFORMATTER_HH_
