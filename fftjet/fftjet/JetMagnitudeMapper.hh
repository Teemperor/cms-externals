//=========================================================================
// JetMagnitudeMapper.hh
//
// Simple jet corrector used to invert a jet response curve
// which is independent from either location or scale
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_JETMAGNITUDEMAPPER_HH_
#define FFTJET_JETMAGNITUDEMAPPER_HH_

#include "fftjet/AbsJetCorrector.hh"
#include "fftjet/LinearInterpolator1d.hh"

namespace fftjet {
    template <typename Jet>
    class JetMagnitudeMapper : public AbsJetCorrector<Jet>
    {
    public:
        // The constructor arguments are as follows:
        //
        //   f               is the jet response curve (ratio of
        //                   reconstructed/actual jet magnitude). Should
        //                   define "double operator()(double) const".
        //
        //   maxMagnitude    is the maximum value of _actual_
        //                   magnitude for which the jet response
        //                   curve is build. Beyond this magnitude
        //                   the jet response is assumed to be constant
        //                   and equal to its value at maxMagnitude
        //
        //   npoints         number of points to use for modeling
        //                   the inverse curve between 0 and the
        //                   reconstructed value which corresponds
        //                   to maxMagnitude
        //
        template <class Functor>
        JetMagnitudeMapper(const Functor& f, double maxMagnitude,
                           unsigned npoints);

        inline virtual ~JetMagnitudeMapper() {delete interp;}

        inline double operator()(const double& magnitude, const Jet&) const
            {return magnitude*(*interp)(magnitude);}

    private:
        JetMagnitudeMapper();
        JetMagnitudeMapper(const JetMagnitudeMapper&);
        JetMagnitudeMapper& operator=(const JetMagnitudeMapper&);

        LinearInterpolator1d* interp;
    };
}

#include "fftjet/JetMagnitudeMapper.icc"

#endif // FFTJET_JETMAGNITUDEMAPPER_HH_
