//=========================================================================
// KernelRecombinationAlg.hh
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
// I. Volobouev
// March 2008
//=========================================================================

#ifndef FFTJET_KERNELRECOMBINATIONALG_HH_
#define FFTJET_KERNELRECOMBINATIONALG_HH_

#include <climits>

#include "fftjet/AbsRecombinationAlg.hh"
#include "fftjet/AbsKernel2d.hh"
#include "fftjet/AbsMembershipFunction.hh"
#include "fftjet/SimpleFunctors.hh"

namespace fftjet {
    template
    <
        typename Real,
        typename VectorLike,
        typename BgData,
        typename VBuilder
    >
    class KernelRecombinationAlg : 
        public AbsRecombinationAlg<Real,VectorLike,BgData>
    {
    public:
        // The meaning of the constructor arguments is as follows:
        //
        // kernel  -- The default kernel functor used to represent
        //            the jet membership function centered at each
        //            peak position. This default functor will be
        //            used whenever the functor associated with a Peak
        //            object is not set. Here, we can use an instance
        //            of any class derived either from AbsKernel2d or
        //            from AbsMembershipFunction. The default functor
        //            can be NULL in which case all individual Peak
        //            functors must be provided.
        //
        // bgWeight  -- Background/noise membership functor. Must
        //              implement the operator
        //
        //              "double operator()(const double& et,
        //                                 const BgData& bg) const"
        //
        //              which returns the function value. "et" argument
        //              will be set to the cell value from the event
        //              energy discretization grid. "bg" argument will be
        //              set to the corresponding background description.
        //
        // unlikelyBgWeight -- If the membership function values
        //                     for each jet are 0 and the background
        //                     weight is smaller than this parameter,
        //                     the cell is a poor fit to the probability
        //                     model used. For membership functions with
        //                     explicit cell Et argument, this condition
        //                     triggers an alternative handling of cell
        //                     energy: the algorithm attempts to split
        //                     the energy between several sources.
        //
        // dataCutoff     -- Only the data points with values above
        //                   this cutoff will contribute to the final
        //                   clusters. The cutoff is not intended
        //                   for noise suppression. Instead, it
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
        //                   to multiple jets and background using weights
        //                   ("fuzzy" clustering).
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
        // This class does not assume ownership of the jet and background
        // membership functors. It is a responsibility of the user of this
        // class to make sure that the "run" method is called only when
        // the membership functor objects are still alive.
        //
        KernelRecombinationAlg(ScaleSpaceKernel* kernel,
                               const Functor2<double,double,BgData>* bgWeight,
                               double unlikelyBgWeight,
                               double dataCutoff,
                               bool winnerTakesAll,
                               bool buildCorrelationMatrix,
                               bool buildClusterMask,
                               unsigned etaBinMin=0,
                               unsigned etaBinMax=UINT_MAX);
        virtual ~KernelRecombinationAlg();

        // In this particular algorithm, the "etaWidth" and "phiWidth" of
        // the output jets will be estimated with respect to the Et centroid
        // direction (minimal width), not the jet direction. If you need it
        // w.r.t. the jet direction, add in quadrature the angular distance
        // from the jet direction to the centroid.
        virtual int run(const std::vector<Peak>& peaks,
                        const Grid2d<Real>& eventData,
                        const BgData* bgData, unsigned nBgEta, unsigned nBgPhi,
                        std::vector<RecombinedJet<VectorLike> >* outputJets,
                        VectorLike* unclustered, double* unused);

        virtual inline unsigned getLastNEta() const {return lastNEta;}
        virtual inline unsigned getLastNPhi() const {return lastNPhi;}
        virtual inline unsigned getLastNJets() const {return passedJets;}

        virtual const unsigned* getClusterMask() const;

    protected:
        // Override the following function to tell the
        // algorithm whether it should build jets using 4-vectors
        inline virtual bool recombine4Vectors() {return true;}

        // Override the following function to tell the
        // algorithm whether it should build 4-vectors
        // out of Et centroids (this is only used when
        // the 4-vectors are not in use).
        inline virtual bool useEtCentroid() {return true;}

        // The following function builds the output 4-vectors
        void buildOutput(std::vector<RecombinedJet<VectorLike> >* outputJets,
                         bool setRecoScaleRatiosToZero);

        // If the following function returns "true", we are done
        // with this event
        bool performPreliminaryProcessing(
            const std::vector<Peak>& peaks,
            const Grid2d<Real>& eventData,
            const BgData* bgData, unsigned nBgEta, unsigned nBgPhi,
            std::vector<RecombinedJet<VectorLike> >* outJets,
            VectorLike* unclustered, double* unused);

        // Input parameters
        ScaleSpaceKernel* const defaultKernel;
        const Functor2<double,double,BgData>* const bgWeightCalc;
        const double unlikelyBgWeight;
        const Real dataCutoff;
        const bool winnerTakesAll;
        const bool buildCorrelationMatrix;
        const bool buildClusterMask;
        const unsigned etaBinMin;
        const unsigned etaBinMax0;

        // Vector builder
        VBuilder vMaker;

        // Various useful buffers with the size
        // dependent on the number of input peaks
        double* weights;
        double* kernelScales;
        double* clusterScales;
        double* clusterScaleRatios;
        double* dEta;
        double* energySum;
        double* energyVariance;
        double* weightSum;
        double* etaEnergySum;
        double* phiEnergySum;
        double* etaPhiESum;
        double* etaSum;
        double* phiSum;
        unsigned* clusterMap;
        const Peak** peakPositions;
        AbsKernel2d** memFcns2d;
        AbsMembershipFunction** memFcns3d;
        VectorLike* jets;
        unsigned nWeights;

        // Cluster mask buffer and its size
        unsigned* clusterMask;
        unsigned maskMemSize;

        // Buffer for dphi values and its size
        double* dphiBuffer;
        unsigned dphiBufSize;

        // Using 2d or 3d membership functions?
        bool use3d;

        // Various variables from the last event
        unsigned lastNJets;
        unsigned lastNEta;
        unsigned lastNPhi;
        unsigned passedJets;
        unsigned bgDim;

    private:
        KernelRecombinationAlg();
        KernelRecombinationAlg(const KernelRecombinationAlg&);
        KernelRecombinationAlg& operator=(const KernelRecombinationAlg&);

        void remapClusterMask();
        void setupBuffers(unsigned njets, unsigned nEta, unsigned nPhi);
        void calculateDphi(const Grid2d<Real>& eventData);
        void allUnclustered(const Grid2d<Real>& evData,
                            VectorLike* unclus, double* unused);
        void processKernelScales(const std::vector<Peak>& peaks);
    };
}

#include "fftjet/KernelRecombinationAlg.icc"

#endif // FFTJET_KERNELRECOMBINATIONALG_HH_
