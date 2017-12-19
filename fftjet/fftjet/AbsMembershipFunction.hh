//=========================================================================
// AbsMembershipFunction.hh
//
// Interface class for jet recombination membership functions
// which depend on both angular and energy variables
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_ABSMEMBERSHIPFUNCTION_HH_
#define FFTJET_ABSMEMBERSHIPFUNCTION_HH_

#include "fftjet/ScaleSpaceKernel.hh"

namespace fftjet {
    struct AbsMembershipFunction : public ScaleSpaceKernel
    {
        virtual ~AbsMembershipFunction() {}

        // The following member sets eta to phi (or x to y) scale ratio
        virtual void setScaleRatio(double r) = 0;

        // The function value
        virtual double operator()(double eta, double phi,
                                  double et, double scale) const = 0;

        // Maximum energy (Et, Pt, etc) which can be absorbed by
        // the jet at the given angle and scale. Override this function
        // if your jet recombination algorithm can meaningfully employ
        // the model in which the cell energy can be shared (actually
        // split, not just use probabilities that the cell belongs
        // to one jet or another) between several jets.
        virtual double absorbableEnergy(double /* eta */, double /* phi */,
                                        double /* scale */) const {return 0.0;}
    };
}

#endif // FFTJET_ABSMEMBERSHIPFUNCTION_HH_
