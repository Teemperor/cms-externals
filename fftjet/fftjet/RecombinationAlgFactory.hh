//=========================================================================
// RecombinationAlgFactory.hh
//
// Factories for recombination algorithms.
//
// Intended usage is like this: once you know which classes to use
// in place of Real, VectorLike, BgData, and VBuilder, instantiate
// the "DefaultRecombinationAlgFactory". You can then choose which
// algorithm to use dynamically, by name, as follows:
//
// std::string algoName;
// ....  Get the algoName somehow ....
// AbsRecombinationAlg<Real,VectorLike,BgData>* alg = 0;
// DefaultRecombinationAlgFactory<...> aFactory;
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
// algoName:           Actual class name:
//
// "Kernel"            KernelRecombinationAlg
// "EtCentroid"        EtCentroidRecombinationAlg
// "EtSum"             EtSumRecombinationAlg
// "FasterKernel"      FasterKernelRecombinationAlg
// "FasterEtCentroid"  FasterEtCentroidRecombinationAlg
// "FaterEtSum"        FasterEtSumRecombinationAlg
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_RECOMBINATIONALGFACTORY_HH_
#define FFTJET_RECOMBINATIONALGFACTORY_HH_

#include <string>
#include <map>

#include "fftjet/ScaleSpaceKernel.hh"
#include "fftjet/AbsRecombinationAlg.hh"
#include "fftjet/SimpleFunctors.hh"

namespace fftjet {
    template<typename Real, typename VectorLike, typename BgData>
    class AbsRecombinationAlgFactory
    {
    public:
        virtual ~AbsRecombinationAlgFactory() {}

        virtual AbsRecombinationAlg<Real,VectorLike,BgData>* create(
            ScaleSpaceKernel* kernel,
            const Functor2<double,double,BgData>* bgWeight,
            double unlikelyBgWeight,
            double dataCutoff,
            bool winnerTakesAll,
            bool buildCorrelationMatrix,
            bool buildClusterMask,
            unsigned etaBinMin=0,
            unsigned etaBinMax=UINT_MAX) const = 0;
    };

    template
    <
        template
        <
            typename Real,
            typename VectorLike,
            typename BgData,
            typename VBuilder
        > class RecombinationAlg,
        typename Real,
        typename VectorLike,
        typename BgData,
        typename VBuilder
    >
    class RecombinationAlgFactory :
        public AbsRecombinationAlgFactory<Real, VectorLike, BgData>
    {
    public:
        virtual ~RecombinationAlgFactory() {}

        RecombinationAlg<Real,VectorLike,BgData,VBuilder>* create(
            ScaleSpaceKernel* kernel,
            const Functor2<double,double,BgData>* bgWeight,
            double unlikelyBgWeight,
            double dataCutoff,
            bool winnerTakesAll,
            bool buildCorrelationMatrix,
            bool buildClusterMask,
            unsigned etaBinMin=0,
            unsigned etaBinMax=UINT_MAX) const
        {
            return new RecombinationAlg<Real,VectorLike,BgData,VBuilder>(
                kernel, bgWeight, unlikelyBgWeight, dataCutoff,
                winnerTakesAll, buildCorrelationMatrix, buildClusterMask,
                etaBinMin, etaBinMax);
        }
    };

    template
    <
        typename Real,
        typename VectorLike,
        typename BgData,
        typename VBuilder
    >
    class DefaultRecombinationAlgFactory :
        public std::map<std::string,
                        AbsRecombinationAlgFactory<Real,VectorLike,BgData>*>
    {
    public:
        DefaultRecombinationAlgFactory();
        virtual ~DefaultRecombinationAlgFactory();

    private:
        typedef AbsRecombinationAlgFactory<Real,VectorLike,BgData> BaseFactory;

        DefaultRecombinationAlgFactory(
            const DefaultRecombinationAlgFactory&);
        DefaultRecombinationAlgFactory& operator=(
            const DefaultRecombinationAlgFactory&);
    };
}

#include "fftjet/RecombinationAlgFactory.icc"

#endif // FFTJET_RECOMBINATIONALGFACTORY_HH_
