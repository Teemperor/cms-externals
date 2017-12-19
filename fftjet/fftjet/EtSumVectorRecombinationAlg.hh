//=========================================================================
// EtSumVectorRecombinationAlg.hh
//
// Class for fuzzy/crisp recombination algorithms in which
// the proximity to the peak is determined by a kernel function.
//
// It is not intended for this class to be constructed and destroyed
// often -- it does too many allocations/deallocations of memory
// buffers to work efficiently in this mode. Instead, create one
// instance of this class at the beginning of your event processing
// loop and call the "run" function for each event.
//
// The "VBuilder" functor on which this class is templated must
// implement a method with the following prototype:
//
// VectorLike operator()(double energyLike, double eta, double phi) const;
//
// The "BgData" class should contain all the info necessary for
// calculating the background weight.
//
// In this algorithm, the eta and phi positions of the jet are taken 
// from the eta and phi of the precluster.
//
// I. Volobouev
// June 2009
//=========================================================================

#ifndef FFTJET_ETSUMVECTORRECOMBINATIONALG_HH_
#define FFTJET_ETSUMVECTORRECOMBINATIONALG_HH_

#include "fftjet/KernelVectorRecombinationAlg.hh"

namespace fftjet {
    template <typename VectorLike, typename BgData, typename VBuilder>
    class EtSumVectorRecombinationAlg :
        public KernelVectorRecombinationAlg<VectorLike, BgData, VBuilder>
    {
    public:
        typedef double (VectorLike::* VectorLikeMemberFunction)() const;

        inline EtSumVectorRecombinationAlg(
            ScaleSpaceKernel* kernel,
            VectorLikeMemberFunction etFcn,
            VectorLikeMemberFunction etaFcn,
            VectorLikeMemberFunction phiFcn,
            const Functor2<double,double,BgData>* bgWeight,
            double unlikelyBgWeight,
            bool winnerTakesAll,
            bool buildCorrelationMatrix,
            bool buildClusterMask)
            : KernelVectorRecombinationAlg<VectorLike, BgData, VBuilder>(
                kernel, etFcn, etaFcn, phiFcn, bgWeight, unlikelyBgWeight,
                winnerTakesAll, buildCorrelationMatrix, buildClusterMask) {}
        inline virtual ~EtSumVectorRecombinationAlg() {}

    protected:
        inline virtual bool recombine4Vectors() {return false;}
        inline virtual bool useEtCentroid() {return false;}
    };
}

#endif // FFTJET_ETSUMVECTORRECOMBINATIONALG_HH_
