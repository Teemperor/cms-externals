//=========================================================================
// ClusteringSequencer.hh
//
// Generic sequencer of steps for constructing the clustering tree
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_CLUSTERINGSEQUENCER_HH_
#define FFTJET_CLUSTERINGSEQUENCER_HH_

#include <vector>

#include "fftjet/AbsConvolverBase.hh"
#include "fftjet/AbsClusteringTree.hh"
#include "fftjet/SimpleFunctors.hh"
#include "fftjet/PeakFinder.hh"
#include "fftjet/Grid2d.hh"

namespace fftjet {
    template<typename Real>
    class ClusteringSequencer
    {
    public:
        // The meaning of the constructor arguments is as follows:
        //
        // convolver     -- Object which calculates and manages Fourier
        //                  transforms and convolutions.
        //
        // selector      -- Object which defines whether the peak
        //                  should be included in the clustering tree.
        //                  The code will try to cast the selector
        //                  dynamically to type AbsPeakSelector<Real>.
        //                  If the cast succeeds, the "setEventData"
        //                  method of the selector will be called.
        //
        // peakFinder    -- Object which performs peak finding
        //
        // initialScales -- Initial set of scales for clustering
        //
        // maxAdaptiveScales -- Maximum number of scales to use for
        //                      growing the tree on top of the initial
        //                      set of scales
        //
        // minRatioLog    -- This argument will be passed to the
        //                   "nextBestScale" function of the clustering
        //                   tree. Used only when "maxAdaptiveScales"
        //                   is not 0.
        //
        // This class does not takes ownership of the "convolver" and
        // "selector" objects. It is a responsibility of the user of this
        // class to make sure that the "processEvent" method is called
        // only when these objects are alive.
        //
        ClusteringSequencer(
            AbsConvolverBase<Real>* convolver,
            Functor1<bool,Peak>* selector,
            const PeakFinder& peakFinder,
            const std::vector<double>& initialScales,
            unsigned maxAdaptiveScales=0, double minRatioLog=0.01);
        virtual ~ClusteringSequencer();

        // The main data processing function. It returns a status word.
        // 0 means everything is OK.
        virtual int run(const Grid2d<Real>& eventData,
                        AbsClusteringTree<Peak,long>* outputTree);

        // The following function inserts the complete event into
        // the clustering tree (but only those grid points whose
        // energy, pt, etc. is above the "dataCutoff" and also within
        // the eta bin limits of the peak finder). Status word is
        // returned (0 mean OK). The "scale" argument must be positive
        // but smaller than any of the scales already present in the tree.
        // This function will be useful in case the clustering tree
        // will be employed in the data analysis as a balltree.
        //
        // Note that Laplacian, Hessian, lifetime, etc. properties
        // will not be set in any meaningful way for this tree layer
        // because raw cells do not necessarily correspond to any peaks.
        //
        // This function causes calculation of the radii and separations
        // for the clusters in the tree if the flag "updateRadii" is true.
        //
        virtual int insertCompleteEvent(
            double scale,
            const Grid2d<Real>& eventData,
            AbsClusteringTree<Peak,long>* outputTree,
            double dataCutoff=0.0,
            bool updateRadii=true);

        // Various inspectors
        inline unsigned maxAdaptiveScales() const {return maxAdaptiveScales_;}
        inline double minRatioLog() const {return minRatioLog_;}
        inline unsigned nInitialScales() const {return initialScales.size();}
        inline const std::vector<double>& getInitialScales() const
            {return initialScales;}
        const PeakFinder& getPeakFinder() const {return peakFinder;}

    protected:
        virtual int processScale(double scale,
                                 const Grid2d<Real>& eventData,
                                 AbsClusteringTree<Peak,long>* outTree);
        void runPeakSelector(double scale);

        AbsConvolverBase<Real>* const convolver;
        Functor1<bool,Peak>* const peakSelector;
        PeakFinder peakFinder;
        std::vector<double> initialScales;
        const unsigned maxAdaptiveScales_;
        const double minRatioLog_;

        const unsigned nEta;
        const unsigned nPhi;
        std::vector<Peak> peaks;
        std::vector<Peak> selected;

        Real* convolvedData;

    private:
        ClusteringSequencer();
        ClusteringSequencer(const ClusteringSequencer&);
        ClusteringSequencer& operator=(const ClusteringSequencer&);
    };
}

#include "fftjet/ClusteringSequencer.icc"

#endif // FFTJET_CLUSTERINGSEQUENCER_HH_
