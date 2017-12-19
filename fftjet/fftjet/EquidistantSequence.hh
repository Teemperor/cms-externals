//========================================================================
// EquidistantSequence.hh
//
// A sequence of points equidistant in linear or log space
//
// I. Volobouev
// March 2009
//========================================================================

#ifndef FFTJET_EQUIDISTANTSEQUENCE_HH_
#define FFTJET_EQUIDISTANTSEQUENCE_HH_

#include <vector>

namespace fftjet {
    class EquidistantInLinearSpace : public std::vector<double>
    {
    public:
        EquidistantInLinearSpace(double minScale, double maxScale,
                                 unsigned nScales);
        virtual ~EquidistantInLinearSpace() {}

    private:
        EquidistantInLinearSpace();
    };

    class EquidistantInLogSpace : public std::vector<double>
    {
    public:
        EquidistantInLogSpace(double minScale, double maxScale,
                              unsigned nScales);
        virtual ~EquidistantInLogSpace() {}

    private:
        EquidistantInLogSpace();
    };
}

#endif // FFTJET_EQUIDISTANTSEQUENCE_HH_
