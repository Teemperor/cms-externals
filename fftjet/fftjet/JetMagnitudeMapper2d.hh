//=========================================================================
// JetMagnitudeMapper2d.hh
//
// This jet corrector can used to invert a jet response curve
// which depends on one of the jet properties (for example, scale or eta)
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_JETMAGNITUDEMAPPER2D_HH_
#define FFTJET_JETMAGNITUDEMAPPER2D_HH_

#include "fftjet/AbsJetCorrector.hh"
#include "fftjet/LinearInterpolator2d.hh"
#include "fftjet/SimpleFunctors.hh"

namespace fftjet {
    template <typename Jet>
    class JetMagnitudeMapper2d : public AbsJetCorrector<Jet>
    {
    public:
        // The constructor arguments are as follows:
        //
        //   f               is the jet response curve (ratio of
        //                   reconstructed/actual jet magnitude). Should
        //                   define "double operator()(double, double) const".
        //                   The first variable is assumed to be the
        //                   predictor on which the response depends
        //                   (such as scale or eta) and the second
        //                   variable is the jet magnitude.
        //
        //   predictorCalc   Pointer to the functor which calculates
        //                   the predictor variable using the jet.
        //
        //   ownPredictor    Set this to "true" if the predictor should be
        //                   deleted in the destructor.
        //
        //   minPredictor,   Together, these three variables describe
        //   maxPredictor,   how to create the grid in the response
        //   nPredPoints     predictor dimension.
        //
        //   maxMagnitude    is the maximum value of _actual_
        //                   magnitude for which the jet response
        //                   curve is build. Beyond this magnitude
        //                   the jet response is assumed to be constant
        //                   and equal to its value at maxMagnitude
        //                   (for each predictor value separately).
        //
        //   nMagPoints      number of reconstructed magnitude points
        //                   to use for modeling the inverse curve between 0
        //                   and the reconstructed value which corresponds
        //                   to maxMagnitude
        //
        template <class Functor2d>
        JetMagnitudeMapper2d(const Functor2d& f,
                             const Functor1<double,Jet>* predictorCalc,
                             bool ownPredictor,
                             double minPredictor,
                             double maxPredictor,
                             unsigned nPredPoints,
                             double maxMagnitude,
                             unsigned nMagPoints);

        virtual ~JetMagnitudeMapper2d();

        inline double operator()(const double& magnitude, const Jet& j) const
            {return magnitude*(*interp)((*fcn)(j), magnitude);}

    private:
        JetMagnitudeMapper2d();
        JetMagnitudeMapper2d(const JetMagnitudeMapper2d&);
        JetMagnitudeMapper2d& operator=(const JetMagnitudeMapper2d&);

        LinearInterpolator2d* interp;
        const Functor1<double,Jet>* fcn;
        const bool ownFcn;
    };
}

#include "fftjet/JetMagnitudeMapper2d.icc"

#endif // FFTJET_JETMAGNITUDEMAPPER2D_HH_
