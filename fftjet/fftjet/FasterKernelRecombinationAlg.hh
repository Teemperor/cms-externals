//=========================================================================
// FasterKernelRecombinationAlg.hh
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
// VectorLike operator()(Real energyLike, Real eta, Real phi) const;
//
// This method builds VectorLike objects (e.g., 4-vectors) out of grid
// points. These objects are later agglomerated by the recombination
// algorithm. There is no abstract base class for VBuilder because it is
// used inside a pretty tight loop where execution speed is important.
//
// The "BgData" class should contain all the info necessary for
// calculating the background weight.
//
// This code usually runs faster than "KernelRecombinationAlg" because
// it maintains a collection of kernel scans internally. Because of this,
// it does not have to reevaluate the membership function for each
// jet, instead it looks it up in a table of precomputed values. This
// internal cacheing results in a useful acceleration if the kernel
// function is expensive to evaluate (for simple kernel functions
// the timing of the algorithm is dominated by other code). The
// catch is that the precluster positions are assumed to fall in the
// centers of the grid cells. This results in a small additional
// uncertainty due to binning (which can often be ignored).
//
// Note that this recombination algorithm can not meaningfully use
// kernels whose eta to phi bandwidth ratio is supposed to change
// from one peak to another.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_FASTERKERNELRECOMBINATIONALG_HH_
#define FFTJET_FASTERKERNELRECOMBINATIONALG_HH_

#include "fftjet/KernelRecombinationAlg.hh"

namespace fftjet {
    template
    <
        typename Real,
        typename VectorLike,
        typename BgData,
        typename VBuilder
    >
    class FasterKernelRecombinationAlg :
        public KernelRecombinationAlg<Real, VectorLike, BgData, VBuilder>
    {
    public:
        // The meaning of the constructor arguments is as follows:
        //
        // kernel -- The kernel functor used to represent the cluster
        //           membership function centered at each peak position.
        //           The functor associated with each peak must not
        //           be set, and the kernel itself must be derived
        //           from "AbsKernel2d".
        //
        // bgWeight  -- Background/pileup membership functor. Must
        //              implement the operator
        //              "double operator()(const double& et,
        //                                 const BgData& bg) const"
        //              which calculates the background weight.
        //              "et" argument will be set to the cell value
        //              from the event data grid. "bg" argument will be
        //              set to the corresponding background description.
        //
        // unlikelyBgWeight -- Reserved for future use (currently ignored)
        //
        // dataCutoff     -- Only the data points with values above
        //                   this cutoff will contribute to the final
        //                   clusters. The cutoff is not intended
        //                   for background suppression. Instead, it
        //                   should be used to skip a large number of
        //                   zero-level points in the data grid which
        //                   are present in case the data is sparse.
        //                   In this case the cutoff should be set to 0.0.
        //                   If the data is not sparse, it is best to set
        //                   this cutoff to a negative number of large
        //                   magnitude.
        //
        // winnerTakesAll -- If "true", there will be one-to-one
        //                   mapping of grid cells to clustered jets or
        //                   to background. That is, each cell will be
        //                   assigned to one jet only ("crisp" clustering).
        //                   If "false", each grid cell will be assigned
        //                   to multiple jets/background using weights.
        //
        // buildCorrelationMatrix -- This parameter is reserved for
        //                           future use (currently ignored).
        //
        // buildClusterMask -- If "true", the code will remember how
        //                     grid cells are assigned to clustered jets
        //                     for the last processed event. The assignments
        //                     are calculated as if the "winnerTakesAll"
        //                     parameter is "true" no matter what its actual
        //                     value is. The mask can be later retrieved
        //                     using the "getClusterMask" function.
        //
        // etaBinMin -- The minimum grid bin number in eta. The code
        //              will not attempt to cluster the cells below
        //              that bin number.
        //
        // etaBinMax -- The maximum grid bin number in eta. The code
        //              will not attempt to cluster the cells for that
        //              bin number or larger. If this parameter is
        //              larger than the grid size then the grid size
        //              will be used internally.
        //
        // This class does not assume ownership of the kernel object.
        // It is a responsibility of the user of this class to make
        // sure that the "run" method is called only when the kernel
        // object is still alive.
        //
        FasterKernelRecombinationAlg(ScaleSpaceKernel* kernel,
                                const Functor2<double,double,BgData>* bgWeight,
                                double unlikelyBgWeight,
                                double dataCutoff,
                                bool winnerTakesAll,
                                bool buildCorrelationMatrix,
                                bool buildClusterMask,
                                unsigned etaBinMin=0,
                                unsigned etaBinMax=UINT_MAX);
        virtual ~FasterKernelRecombinationAlg();

        virtual int run(const std::vector<Peak>& peaks,
                        const Grid2d<Real>& eventData,
                        const BgData* bgData, unsigned nBgEta, unsigned nBgPhi,
                        std::vector<RecombinedJet<VectorLike> >* outJets,
                        VectorLike* unclustered, double* unused);
    private:
        typedef KernelRecombinationAlg<Real,VectorLike,BgData,VBuilder> B;

        // Nested class to manage kernel scan collections
        class KernelScanCollection
        {
        public:
            explicit KernelScanCollection(const AbsKernel2d* kernel);
            virtual ~KernelScanCollection();
            
            void mapCoords(const Grid2d<Real>& eventData, unsigned etaBinMin,
                           unsigned etaBinMax, const Peak** peakPositions,
                           unsigned njets, unsigned **etaBuf, unsigned **phiBuf);
            const Real* getKernelScan(const Grid2d<Real>& eventData, double scale);
            
        private:
            KernelScanCollection();
            KernelScanCollection(const KernelScanCollection&);
            KernelScanCollection& operator=(const KernelScanCollection&);
            
            const AbsKernel2d* const kernel;
            
            unsigned *mappedEta;
            unsigned mapBufLen;
            
            std::map<double, Grid2d<Real>*> kernelScans;
        };

        FasterKernelRecombinationAlg();
        FasterKernelRecombinationAlg(const FasterKernelRecombinationAlg&);
        FasterKernelRecombinationAlg& operator=(
            const FasterKernelRecombinationAlg&);

        void setupTableBuf(unsigned njets);

        KernelScanCollection kernelScans;

        const Real** scanTable;
        unsigned scanTableLen;
    };
}

#include "fftjet/FasterKernelRecombinationAlg.icc"

#endif // FFTJET_FASTERKERNELRECOMBINATIONALG_HH_
