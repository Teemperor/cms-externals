//=========================================================================
// AbsRecombinationAlg.hh
//
// Interface class for fuzzy/crisp jet recombination algorithms
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_ABSRECOMBINATIONALG_HH_
#define FFTJET_ABSRECOMBINATIONALG_HH_

#include <vector>

#include "fftjet/Grid2d.hh"
#include "fftjet/RecombinedJet.hh"

namespace fftjet {
    template <typename Real, typename VectorLike, typename BgData>
    class AbsRecombinationAlg
    {
    public:
        virtual ~AbsRecombinationAlg() {}

        // Prototype for the main recombination algorithm function.
        // The derived classes must provide an implementation. Arguments
        // are as follows:
        //
        //  peaks -- Collection of preclusters from the pattern recognition
        //           stage (produced by the "PeakFinder" or similar class).
        //           Implementations of this algorithm should ignore peaks
        //           for which "membershipScaleFactor()" member function
        //           returns 0.0.
        //
        //  eventData -- Actual discretized data from the current event.
        //               Usually this is the same Grid2d object as the one
        //               used for pattern recognition.
        //
        //  bgData  -- Some representation of the average pile-up energy,
        //             noise, background, etc. Exact interpretation of 
        //             this argument is up to concrete implementations
        //             of this class.
        //
        //  nBgEta, -- Dimensionalities of the "bgData" array. Normally,
        //  nBgPhi     the derived classes should support the following
        //             interpretation of these variables:
        //
        //             nBgEta == 1, nBgPhi == 1. This means that "bgData" is
        //             a scalar, and the same background description will be
        //             used for each grid cell.
        //
        //             nBgEta == eventData.nEta(), nBgPhi == 1. This means
        //             that "bgData" is a 1-d array. The same background
        //             description will be used for each phi cell number
        //             in each given eta, but backgrounds are different
        //             for different eta.
        //
        //             nBgEta == eventData.nEta(), nBgPhi == eventData.nPhi()
        //             This means that "bgData" is a 2-d array, and there is
        //             a separate background description object for every
        //             grid cell.
        //
        //             All other nBgEta, nBgPhi combinations should be
        //             deemed illegal. The "bgDimensionality" method
        //             can be used to enforce this policy.
        //
        //  jets   -- Output jets.
        //
        //  unclustered -- Vector sum of grid cells not included
        //                 in any particular jet (produced on output).
        //
        //  unused -- Scalar sum of values in the event grid not included
        //            in any particular jet (produced on output).
        //
        // The function should return a status word. 0 means everything
        // is fine. The meaning of other codes is up to concrete
        // implementations.
        //
        virtual int run(const std::vector<Peak>& peaks,
                        const Grid2d<Real>& eventData,
                        const BgData* bgData, unsigned nBgEta, unsigned nBgPhi,
                        std::vector<RecombinedJet<VectorLike> >* jets,
                        VectorLike* unclustered, double* unused) = 0;

        // Get some info about the data last processed.
        // Useful for retrieving the energy correlation
        // matrix and the cell-to-jet assignment mask.
        virtual unsigned getLastNEta() const = 0;
        virtual unsigned getLastNPhi() const = 0;
        virtual unsigned getLastNJets() const = 0;

        // The following function should return the jet energy
        // covariance matrix if it is available. In case it is not
        // available, the function should return the null pointer.
        // The data should remain valid at least until the next
        // "run" call (or until this object is destroyed).
        // The returned array should be dimensioned [njets+1][njets+1],
        // where "njets" is the number of jets returned by the
        // previous "run" call. C storage convention should be used
        // (row-major). Index 0 in this matrix should be reserved
        // for the unclustered energy.
        virtual const double* getEnergyCovarianceMatrix() const {return 0;}

        // The following call should return the cell-to-jet mask
        // if it is available. In case it is not available,
        // the function should return the null pointer. The data
        // should remain valid at least until the next "run" call.
        // The returned array should be dimensioned [nEta][nPhi]
        // with sizes corresponding to the data grid size from
        // the previous "run" call.
        //
        // clusterMask[ieta][iphi] should be set to the jet number
        // plus one. Jet numbering should be performed according to
        // the jet sequence returned by "run". Number 0 is reserved
        // for the unclustered energy.
        virtual const unsigned* getClusterMask() const {return 0;}

        // The following function will checks the "nBgEta" and
        // "nBgPhi" arguments and will return the effective
        // dimensionality of the "bgData" array (0, 1, or 2).
        // A run-time error will be generated in case nBgEta, nBgPhi
        // arguments are inconsistent with the Grid2d dimensions.
        static unsigned bgDimensionality(const Grid2d<Real>& eventData,
                                         unsigned nBgEta, unsigned nBgPhi);
    };
}

#include "fftjet/AbsRecombinationAlg.icc"

#endif // FFTJET_ABSRECOMBINATIONALG_HH_
