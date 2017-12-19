//=========================================================================
// PeakSelectors.hh
//
// Some simple functors for peak selection
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_PEAKSELECTORS_HH_
#define FFTJET_PEAKSELECTORS_HH_

#include <cfloat>
#include "fftjet/SimpleFunctors.hh"

namespace fftjet {
    // Simplistic, scale and data-independent peak selector
    class SimplePeakSelector : public Functor1<bool,Peak>
    {
    public:
        explicit SimplePeakSelector(double magCut,
                                    double driftSpeedCut=DBL_MAX,
                                    double magSpeedCut=-DBL_MAX,
                                    double lifeTimeCut=-DBL_MAX,
                                    double NNDCut=-DBL_MAX,
                                    double etaCut=DBL_MAX,
                                    double splitTimeCut=-DBL_MAX,
                                    double mergeTimeCut=-DBL_MAX);
        virtual ~SimplePeakSelector() {}

        virtual bool operator()(const Peak& peak) const;

    private:
        SimplePeakSelector();

        double magCut_;
        double driftSpeedCut_;
        double magSpeedCut_;
        double lifeTimeCut_;
        double NNDCut_;
        double etaCut_;
        double splitTimeCut_;
        double mergeTimeCut_;
    };

    // Selector based on scale-dependent magnitude. Intended for use
    // before clustering tree is constructed. The magnitude cut formula
    // looks like this:
    //
    // a*pow(scale,p) + b
    //
    // Note that, normally, "a" should be positive, "b" should be
    // non-negative, and "p" should be a negative number smaller
    // or equal to 2 in magnitude. If these conditions are not
    // satisfied, the peaks produced by stand-alone energy depositions
    // will "pop up" with increasing scale rather than disappear.
    // 
    class ScalePowerPeakSelector : public Functor1<bool,Peak>
    {
    public:
        ScalePowerPeakSelector(double a, double p, double b,
                               double etaCut=DBL_MAX);
        virtual ~ScalePowerPeakSelector() {}

        virtual bool operator()(const Peak& peak) const;

    private:
        ScalePowerPeakSelector();

        double a_;
        double p_;
        double b_;
        double etaCut_;
    };
}

#include "fftjet/PeakSelectors.icc"

#endif // FFTJET_PEAKSELECTORS_HH_
