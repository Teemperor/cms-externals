//=========================================================================
// AbsVectorRecombinationAlg.hh
//
// Interface class for fuzzy/crisp jet recombination algorithms in which
// event energy flow is described by a collection of 4-vectors
//
// I. Volobouev
// June 2009
//=========================================================================

#ifndef FFTJET_ABSVECTORRECOMBINATIONALG_HH_
#define FFTJET_ABSVECTORRECOMBINATIONALG_HH_

#include <vector>

#include "fftjet/RecombinedJet.hh"

namespace fftjet {
    template <typename VectorLike, typename BgData>
    class AbsVectorRecombinationAlg
    {
    public:
        virtual ~AbsVectorRecombinationAlg() {}

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
        //  eventData -- A vector of elements which describe the
        //               energy flow in the event. "VectorLike" is
        //               a type (typically, a 4-vector) from which
        //               it is possible to extract Et (or pt),
        //               eta, and phi information.
        //
        //  bgData  -- Some representation of the average pile-up energy,
        //             noise, background, etc. Exact interpretation of 
        //             this argument is up to concrete implementations
        //             of this class.
        //
        //  bgDataLength -- Dimensionality of the "bgData" array. Normally,
        //             the derived classes should support the following
        //             interpretation of this variable:
        //
        //             bgDataLength == 1. This means that "bgData" is
        //             a scalar, and the same background description
        //             will be used for every data point in the event.
        //
        //             bgDataLength == eventData.size(). This means
        //             that "bgData" is a 1-d array. Each data point
        //             has its own associated background.
        //
        //             All other bgDataLength values should be
        //             deemed illegal.
        //
        //  jets   -- Output jets.
        //
        //  unclustered -- Vector sum of data elements not included
        //                 in any particular jet (produced on output).
        //
        //  unused -- Scalar sum of Et (or pt) values in the data not
        //            included in any particular jet (produced on output).
        //
        // The function should return a status word. 0 means everything
        // is fine. The meaning of other codes is up to concrete
        // implementations.
        //
        virtual int run(const std::vector<Peak>& peaks,
                        const std::vector<VectorLike>& eventData,
                        const BgData* bgData, unsigned bgDataLength,
                        std::vector<RecombinedJet<VectorLike> >* jets,
                        VectorLike* unclustered, double* unused) = 0;

        // Get some info about the data last processed.
        // Useful for retrieving the energy correlation
        // matrix and the cell-to-jet assignment mask.
        virtual unsigned getLastNData() const = 0;
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

        // The following call should return the assignments of the
        // input vectors to jets. In case these assignments are not
        // available, the function should return the null pointer.
        // The result should remain valid at least until the next
        // "run" call. The returned array should have "NData"
        // elements consistent with the data length from the previous
        // "run" call.
        //
        // clusterMask[i] should be set to the jet number
        // plus one. Jet numbering should be performed according to
        // the jet sequence returned by "run". Number 0 is reserved
        // for the unclustered energy.
        virtual const unsigned* getClusterMask() const {return 0;}
    };
}

#endif // FFTJET_ABSVECTORRECOMBINATIONALG_HH_
