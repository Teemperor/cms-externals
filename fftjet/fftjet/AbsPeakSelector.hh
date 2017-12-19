//=========================================================================
// AbsPeakSelector.hh
//
// Interface for peak selectors in ClusteringSequencer and
// ConstScaleReconstruction. Derive you peak selector from this class
// if your selector needs to do some calculations using the discretized
// event data before the selection is made. Note that predicates used
// by clustering trees do not get the event data -- they should normally
// be derived from Functor1<bool,Peak>.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_ABSPEAKSELECTOR_HH_
#define FFTJET_ABSPEAKSELECTOR_HH_

#include "fftjet/Peak.hh"
#include "fftjet/Grid2d.hh"
#include "fftjet/SimpleFunctors.hh"

namespace fftjet {
    template<typename Real>
    struct AbsPeakSelector : public Functor1<bool,Peak>
    {
        virtual ~AbsPeakSelector() {}

        // The following function will be called once per event,
        // before the peak selection must be made
        virtual void setEventData(const Grid2d<Real>& eventData) = 0;

        // The following operator should return "true" if the peak
        // is to be included in the subsequent processing
        virtual bool operator()(const Peak& peak) const = 0;
    };

    // A simple selector which selects everything
    struct AllPeaksPass : public Functor1<bool,Peak>
    {
        inline bool operator()(const Peak&) const {return true;}
    };
}

#endif // FFTJET_ABSPEAKSELECTOR_HH_
