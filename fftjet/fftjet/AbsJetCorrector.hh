//=========================================================================
// AbsJetCorrector.hh
//
// Interface class for jet correction-type calculations. It is assumed
// that we do not want to adjust the jet location, only its magnitude.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_ABSJETCORRECTOR_HH_
#define FFTJET_ABSJETCORRECTOR_HH_

#include <cassert>

namespace fftjet {
    template <typename Jet>
    struct AbsJetCorrector
    {
        virtual ~AbsJetCorrector() {}

        // It is assumed that "magnitude" is some kind of a guess
        // for the jet energy-like variable (Et, Pt, etc) which we
        // want to correct. The Jet class should provide all
        // other info about the jet (eta, phi, scale, etc).
        // The function should return adjusted magnitude.
        virtual double operator()(const double& magnitude,
                                  const Jet& jet) const = 0;
    };

    // The following class does not perform any correction
    template <typename Jet>
    struct IdleJetCorrector : public AbsJetCorrector<Jet>
    {
        inline double operator()(const double& magnitude, const Jet&) const
            {return magnitude;}
    };

    // The following can sometimes be used for debugging purposes
    template <typename Jet>
    struct InvalidJetCorrector : public AbsJetCorrector<Jet>
    {
        inline double operator()(const double&, const Jet&) const
        {
            assert(!"InvalidJetCorrector::operator() triggered");
            return 0.0;
        }
    };
}

#endif // FFTJET_ABSJETCORRECTOR_HH_
