//=========================================================================
// ConstScaleReconstruction.hh
//
// Driver class for running both pattern recognition and jet recombination
// at the same constant scale throughout the whole data grid
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_CONSTSCALERECONSTRUCTION_HH_
#define FFTJET_CONSTSCALERECONSTRUCTION_HH_

#include <vector>

#include "fftjet/AbsConvolverBase.hh"
#include "fftjet/SimpleFunctors.hh"
#include "fftjet/PeakFinder.hh"
#include "fftjet/AbsRecombinationAlg.hh"

namespace fftjet {
    template <typename Real, typename VectorLike, typename BgData>
    class ConstScaleReconstruction
    {
    public:
        // This object will not own "convolver", "selector", and
        // "recombiner" pointers. It is a responsibility of the user
        // of this class to make sure those objects are still
        // alive when the "run" method is called.
        //
        // The code will attempt to cast the selector dynamically
        // to type AbsPeakSelector<Real>. If the cast succeeds,
        // the "setEventData" method of the selector will be called.
        //
        // It is possible to use other powers of the energy-like
        // variable besides 1 for pattern recognition (argument "etPower").
        // See the comment to the "power" function inside the "Grid2d.hh"
        // header file.
        //
        ConstScaleReconstruction(
            AbsConvolverBase<Real>* convolver,
            Functor1<bool,Peak>* selector,
            const PeakFinder& peakFinder,
            AbsRecombinationAlg<Real,VectorLike,BgData>* recombiner,
            double etPower = 1.0);
        virtual ~ConstScaleReconstruction();

        // The return value from the following function is the status
        // produced by the AbsRecombinationAlg "run" function
        virtual int run(double scale,
                        const Grid2d<Real>& eventData,
                        const BgData* bgData, unsigned nBgEta, unsigned nBgPhi,
                        std::vector<RecombinedJet<VectorLike> >* outputJets,
                        VectorLike* unclustered, double* unused);
    private:
        ConstScaleReconstruction();
        ConstScaleReconstruction(const ConstScaleReconstruction&);
        ConstScaleReconstruction& operator=(const ConstScaleReconstruction&);

        AbsConvolverBase<Real>* const convolver_;
        Functor1<bool,Peak>* const selector_;
        PeakFinder peakFinder_;
        AbsRecombinationAlg<Real,VectorLike,BgData>* const recombiner_;

        const unsigned nEta_;
        const unsigned nPhi_;
        Real* const convolvedData_;
        const double etPower_;

        std::vector<Peak> peaks;
        Grid2d<Real>* copyGrid;
    };
}

#include "fftjet/ConstScaleReconstruction.icc"

#endif // FFTJET_CONSTSCALERECONSTRUCTION_HH_
