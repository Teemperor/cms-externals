#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>

#include "fftjet/DXGlyphInfo.hh"

namespace fftjet {
    static void fixConnection(std::vector<DXGlyphInfo>* glyphInfo,
                              const unsigned index,
                              const double minColor,
                              const double minSize)
    {
        DXGlyphInfo& glyph1((*glyphInfo)[index]);
        const unsigned ipar(glyph1.parent());
        const DXGlyphInfo& glyph2((*glyphInfo)[ipar]);

        // This function assumes that the phi angles are non-negative
        assert(glyph1.phi() >= 0.0);
        assert(glyph2.phi() >= 0.0);

        // Check which glyph is closer to 0 angle boundary
        const bool dauCloser = glyph1.phi() < glyph2.phi();
        DXGlyphInfo close0, far0;
        if (dauCloser)
        {
            close0 = glyph1;
            far0 = glyph2;
        }
        else
        {
            close0 = glyph2;
            far0 = glyph1;
        }

        // Figure out the intersection point with the boundary
        const double w = close0.phi()/(close0.phi() - far0.phi() + 2.0*M_PI);

        // Build two dummy glyphs at the boundaries
        const double eta = close0.eta() + w*(far0.eta() - close0.eta());
        const double s = close0.scale() + w*(far0.scale() - close0.scale());
        DXGlyphInfo dummy0(eta, 0.0, s, minSize, minColor, -1);
        DXGlyphInfo dummy2pi(eta, 2.0*M_PI, s, minSize, minColor, -1);

        // Make correct connections
        const unsigned nGlyphs = glyphInfo->size();
        if (dauCloser)
        {
            glyph1.setParent(nGlyphs);     // dummy0 is the new parent
            dummy2pi.setParent(ipar);
        }
        else
        {
            glyph1.setParent(nGlyphs + 1); // dummy2pi is the new parent
            dummy0.setParent(ipar);
        }
        glyphInfo->push_back(dummy0);
        glyphInfo->push_back(dummy2pi);
    }


    void fixDXCylindricalGeometry(std::vector<DXGlyphInfo>* glyphInfo,
                                  const double etaRange)
    {
        std::vector<DXGlyphInfo>& glyphs(*glyphInfo);
        const unsigned nGlyphs = glyphs.size();

        // Figure out bad connections, min values
        // for glyph size and color
        double minColor = DBL_MAX, minSize = DBL_MAX;
        double minScale = DBL_MAX, maxScale = -DBL_MAX;
        std::vector<unsigned> badConnections;
        for (unsigned i=0; i<nGlyphs; ++i)
        {
            const DXGlyphInfo& glyph(glyphs[i]);
            const int ipar(glyph.parent());
            if (ipar >= 0)
                if (fabs(glyph.phi() - glyphs[ipar].phi()) > M_PI)
                    badConnections.push_back(i);
            if (glyph.color() < minColor)
                minColor = glyph.color();
            if (glyph.size() < minSize)
                minSize = glyph.size();
            const double s = glyph.scale();
            if (s < minScale)
                minScale = s;
            if (s > maxScale)
                maxScale = s;
        }

        const unsigned nBad = badConnections.size();
        for (unsigned ibad=0; ibad<nBad; ++ibad)
            fixConnection(glyphInfo, badConnections[ibad],
                          minColor, minSize);

        // Add small glyphs to keep the plot size constant
        // between different trees
        if (etaRange)
        {
            DXGlyphInfo dummy0(-etaRange, 0.0, minScale,
                               minSize, minColor, -1);
            DXGlyphInfo dummy2pi(etaRange, 2.0*M_PI, maxScale,
                                 minSize, minColor, -1);
            glyphInfo->push_back(dummy0);
            glyphInfo->push_back(dummy2pi);
        }
    }


    void writeDXGlyphVector(const std::vector<DXGlyphInfo>& glyphInfo,
                            const unsigned long runNum,
                            const unsigned long evNum,
                            std::ostream& os)
    {
        // First, dump glyph coordinates.
        const unsigned nGlyphs = glyphInfo.size();
        {
            os << "object 1 class array type float rank 1 shape 3 items "
               << nGlyphs << " data follows\n";
            for (unsigned i=0; i<nGlyphs; ++i)
            {
                const DXGlyphInfo& glyph(glyphInfo[i]);
                os << glyph.eta() << ' ' << glyph.phi() << ' '
                   << glyph.scale() << '\n';
            }
            os << "attribute \"dep\" string \"positions\"\n\n";
        }

        // Dump glyph sizes
        {
            os << "object 2 class array type float rank 0 items "
               << nGlyphs << " data follows\n";
            for (unsigned i=0; i<nGlyphs; ++i)
                os << glyphInfo[i].size() << '\n';
            os << "attribute \"dep\" string \"positions\"\n\n";
        }

        // Dump glyph colors
        {
            os << "object 3 class array type float rank 0 items "
               << nGlyphs << " data follows\n";
            for (unsigned i=0; i<nGlyphs; ++i)
                os << glyphInfo[i].color() << '\n';
            os << "attribute \"dep\" string \"positions\"\n\n";
        }

        // Dump glyph connections
        {
            unsigned totalDaus = 0;
            for (unsigned i=0; i<nGlyphs; ++i)
                if (glyphInfo[i].parent() >= 0)
                    ++totalDaus;
            os << "object 4 class array type integer rank 1 shape 2 items "
               << totalDaus << " data follows\n";
            for (unsigned i=0; i<nGlyphs; ++i)
                if (glyphInfo[i].parent() >= 0)
                    os << glyphInfo[i].parent() << ' ' << i << '\n';
            os << "attribute \"dep\" string \"connections\"\n";
            os << "attribute \"element type\" string \"lines\"\n";
            os << "attribute \"ref\" string \"positions\"\n\n";
        }
        
        // Collect the data into an OpenDX field object
        {
            os << "object 5 class field\n";
            os << "component \"positions\" value 1\n";
            os << "component \"connections\" value 4\n";
            os << "component \"sizedata\" value 2\n";
            os << "component \"colordata\" value 3\n";
            os << "attribute \"runnumber\" number " << runNum << '\n';
            os << "attribute \"eventnumber\" number " << evNum << '\n';
            os << "\nend";
        }
    }


    void remapDXscales(std::vector<DXGlyphInfo>* glyphs, const double range)
    {
        if (range == 0.0)
            return;

        const unsigned nGlyphs = glyphs->size();
        double minScale = DBL_MAX, maxScale = -DBL_MAX;
        for (unsigned i=0; i<nGlyphs; ++i)
        {
            const double s = (*glyphs)[i].scale();
            if (s > maxScale)
                maxScale = s;
            if (s < minScale)
                minScale = s;
        }
        if (maxScale <= minScale)
            return;

        const double a = range/(maxScale - minScale);
        const double b = -a*minScale;

        for (unsigned i=0; i<nGlyphs; ++i)
        {
            DXGlyphInfo& glyph((*glyphs)[i]);
            glyph.setScale(a*glyph.scale() + b);
        }
    }
}
