//=========================================================================
// EtSumRecombinationAlg.hh
//
// Class for fuzzy/crisp recombination algorithms in which
// the proximity to the peak is determined by a kernel function.
// This particular class assumes that the data grid is populated with
// transverse energy values.
//
// It is not intended for this class to be constructed and destroyed
// often -- it does too many allocations/deallocations of memory
// buffers to work efficiently in this mode. Instead, create one
// instance of this class at the beginning of your event processing
// loop and call the "run" function for each event.
//
// The "VBuilder" functor on which this class is templated must
// implement a method with the following prototype:
// VectorLike operator()(Real Et, Real eta, Real phi) const;
//
// The "BgData" class should contain all the info necessary for
// calculating the background weight.
//
// In this algorithm, the eta and phi positions of the jet are taken 
// from the eta and phi of the precluster.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_ETSUMRECOMBINATIONALG_HH_
#define FFTJET_ETSUMRECOMBINATIONALG_HH_

#include "fftjet/KernelRecombinationAlg.hh"

namespace fftjet {
    template
    <
        typename Real,
        typename VectorLike,
        typename BgData,
        typename VBuilder
    >
    class EtSumRecombinationAlg :
        public KernelRecombinationAlg<Real, VectorLike, BgData, VBuilder>
    {
    public:
        inline EtSumRecombinationAlg(ScaleSpaceKernel* kernel,
                                const Functor2<double,double,BgData>* bgWeight,
                                double unlikelyBgWeight,
                                double dataCutoff,
                                bool winnerTakesAll,
                                bool buildCorrelationMatrix,
                                bool buildClusterMask,
                                unsigned etaBinMin=0,
                                unsigned etaBinMax=UINT_MAX)
            : KernelRecombinationAlg<Real,VectorLike,BgData,VBuilder>(
                kernel, bgWeight, unlikelyBgWeight, dataCutoff,
                winnerTakesAll, buildCorrelationMatrix,
                buildClusterMask, etaBinMin, etaBinMax) {}
        inline virtual ~EtSumRecombinationAlg() {}

    protected:
        inline virtual bool recombine4Vectors() {return false;}
        inline virtual bool useEtCentroid() {return false;}
    };
}

#endif // FFTJET_ETSUMRECOMBINATIONALG_HH_
