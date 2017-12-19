//=========================================================================
// DXGlyphInfo.hh
//
// A little helper class which colects all the info necessary to build
// OpenDX glyphs. Also, this header declares several stand-alone
// glyph-related utility functions.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_DXGLYPHINFO_HH_
#define FFTJET_DXGLYPHINFO_HH_

#include <vector>
#include <iosfwd>

namespace fftjet {
    class DXGlyphInfo
    {
    public:
        inline DXGlyphInfo()
            : eta_(0.0), phi_(0.0), scale_(0.0),
              size_(0.0), color_(0.0), parent_(-1) {}
        inline DXGlyphInfo(double eta, double phi, double scale,
                           double size, double color, int parent)
            : eta_(eta), phi_(phi), scale_(scale),
              size_(size), color_(color), parent_(parent) {}

        // Accessors
        inline double eta() const {return eta_;}
        inline double phi() const {return phi_;}
        inline double scale() const {return scale_;}
        inline double size() const {return size_;}
        inline double color() const {return color_;}
        inline int parent() const {return parent_;}

        // Necessary modifiers
        inline void setParent(const int par) {parent_ = par;}
        inline void setScale(const double s) {scale_ = s;}

    private:
        double eta_;
        double phi_;
        double scale_;
        double size_;
        double color_;
        int parent_;
    };

    // The following function modifies the collection of DXGlyphInfo
    // objects in order to produce short connections through the
    // phi = 0 = 2*pi section.
    void fixDXCylindricalGeometry(std::vector<DXGlyphInfo>* glyphs,
                                  double etaRange);

    // The following function remaps the scales. It does nothing
    // if the "range" parameter is 0.0.
    void remapDXscales(std::vector<DXGlyphInfo>* glyphs, double range);

    // The following function writes out a vector of DXGlyphInfo objects
    void writeDXGlyphVector(const std::vector<DXGlyphInfo>& glyphs,
                            unsigned long runNumber, unsigned long eventNumber,
                            std::ostream& os);
}

#endif // FFTJET_DXGLYPHINFO_HH_
