//=========================================================================
// VectorRecombinationAlgFactory.hh
//
// Factories for algorithms used to recombine collections of 4-vectors
//
// Intended usage is like this: once you know which classes to use
// in place of VectorLike, BgData, and VBuilder, instantiate
// the "DefaultVectorRecombinationAlgFactory". You can then choose
// which algorithm to use dynamically, by name, as follows:
//
// std::string algoName;
// ....  Get the algoName somehow ....
// AbsVectorRecombinationAlg<VectorLike,BgData>* alg = 0;
// DefaultVectorRecombinationAlgFactory<...> aFactory;
// if (aFactory[algoName] == NULL)
// {
//     .... Process invalid algorithm name error .....
// }
// else
//     alg = aFactory[algoName]->create(kernel, ...);
// .... Perform some work with alg ....
// delete alg;
//
// The correspondence between valid algorithm names and classes is
//
// algoName:        Actual class name:
//
// "Kernel"         KernelVectorRecombinationAlg
// "EtCentroid"     EtCentroidVectorRecombinationAlg
// "EtSum"          EtSumVectorRecombinationAlg
//
// I. Volobouev
// June 2009
//=========================================================================

#ifndef FFTJET_VECTORRECOMBINATIONALGFACTORY_HH_
#define FFTJET_VECTORRECOMBINATIONALGFACTORY_HH_

#include <string>
#include <map>

#include "fftjet/ScaleSpaceKernel.hh"
#include "fftjet/AbsVectorRecombinationAlg.hh"
#include "fftjet/SimpleFunctors.hh"

namespace fftjet {
    template<typename VectorLike, typename BgData>
    class AbsVectorRecombinationAlgFactory
    {
    public:
        typedef double (VectorLike::* VectorLikeMemberFunction)() const;

        virtual ~AbsVectorRecombinationAlgFactory() {}

        virtual AbsVectorRecombinationAlg<VectorLike,BgData>* create(
            ScaleSpaceKernel* kernel,
            VectorLikeMemberFunction etFcn,
            VectorLikeMemberFunction etaFcn,
            VectorLikeMemberFunction phiFcn,
            const Functor2<double,double,BgData>* bgWeight,
            double unlikelyBgWeight,
            bool winnerTakesAll,
            bool buildCorrelationMatrix,
            bool buildClusterMask) const = 0;
    };

    template
    <
        template
        <
            typename VectorLike,
            typename BgData,
            typename VBuilder
        > class VectorRecombinationAlg,
        typename VectorLike,
        typename BgData,
        typename VBuilder
    >
    class VectorRecombinationAlgFactory :
        public AbsVectorRecombinationAlgFactory<VectorLike, BgData>
    {
    public:
        typedef double (VectorLike::* VectorLikeMemberFunction)() const;

        virtual ~VectorRecombinationAlgFactory() {}

        VectorRecombinationAlg<VectorLike,BgData,VBuilder>* create(
            ScaleSpaceKernel* kernel,
            VectorLikeMemberFunction etFcn,
            VectorLikeMemberFunction etaFcn,
            VectorLikeMemberFunction phiFcn,
            const Functor2<double,double,BgData>* bgWeight,
            double unlikelyBgWeight,
            bool winnerTakesAll,
            bool buildCorrelationMatrix,
            bool buildClusterMask) const
        {
            return new VectorRecombinationAlg<VectorLike,BgData,VBuilder>(
                kernel, etFcn, etaFcn, phiFcn, bgWeight, unlikelyBgWeight,
                winnerTakesAll, buildCorrelationMatrix, buildClusterMask);
        }
    };

    template
    <
        typename VectorLike,
        typename BgData,
        typename VBuilder
    >
    class DefaultVectorRecombinationAlgFactory :
        public std::map<std::string,
                        AbsVectorRecombinationAlgFactory<VectorLike,BgData>*>
    {
    public:
        DefaultVectorRecombinationAlgFactory();
        virtual ~DefaultVectorRecombinationAlgFactory();

    private:
        typedef AbsVectorRecombinationAlgFactory<VectorLike,BgData> BaseFactory;

        DefaultVectorRecombinationAlgFactory(
            const DefaultVectorRecombinationAlgFactory&);
        DefaultVectorRecombinationAlgFactory& operator=(
            const DefaultVectorRecombinationAlgFactory&);
    };
}

#include "fftjet/VectorRecombinationAlgFactory.icc"

#endif // FFTJET_VECTORRECOMBINATIONALGFACTORY_HH_
