//=========================================================================
// PeakFinder.hh
//
// This class implements an algorithm for locating precluster positions
//
// I. Volobouev
// March 2008
//=========================================================================

#ifndef FFTJET_PEAKFINDER_HH_
#define FFTJET_PEAKFINDER_HH_

#include <vector>
#include <climits>

#include "fftjet/Peak.hh"

namespace fftjet {
    class PeakFinder
    {
    public:
        //
        // The "peakHeightCutoff" parameter of the constructor should be
        // used in case large plateaus are expected at the value of
        // "peakHeightCutoff" or just below. A typical use for this
        // parameter would be to raise the peak threshold above the
        // round-off noise of the FFT algorithm. If such plateaus are not
        // expected, it is best to set the "peakHeightCutoff" parameter to
        // some negative number of large magnitude.
        //
        // Other constructor arguments have the following meaning:
        //
        // subCellResolution -- Set this to "true" in order to improve
        //                      the peak position determination by fitting
        //                      the 3x3 region with the peak in the center
        //                      using a 2-d quadratic polynomial (and then
        //                      finding the peak of that polynomial).
        //
        // minEtaBin -- Minimum eta bin number for peak searching.
        //              Meaningful values of this parameter are 1 and higher.
        //              Value of 0 will be internally turned into 1.
        //
        // maxEtaBin -- Maximum eta bin number for peak searching.
        //              Meaningful values of this parameter are nEta-1
        //              and lower, where "nEta" is the size of the grid
        //              with which the "find" function will be called.
        //
        // printFitWarnings -- If "true", the peak finder will be verbose
        //                     about problems encountered during the peak
        //                     position fitting. This parameter is
        //                     meaningful only when the "subCellResolution"
        //                     value is "true".
        //
        explicit PeakFinder(double peakHeightCutoff,
                            bool subCellResolution=true,
                            unsigned minEtaBin=1, unsigned maxEtaBin=UINT_MAX,
                            bool printFitWarnings=false);
        PeakFinder(const PeakFinder&);
        ~PeakFinder();

        // Change the cutoff for the peak height 
        void setPeakHeightCutoff(double cutoff);

        // Trivial inspectors of the contents
        inline double peakHeightCutoff() const {return peakHeightCutoff_;}
        inline bool subCellResolution() const {return subCellResolution_;}
        inline unsigned smallestEtaBin() const {return minEtaBin_;}
        inline unsigned largestEtaBin() const {return maxEtaBin_;}

        // The following function will not look for peaks at
        // etaBin = 0 and etaBin = nEta-1. The peaks are searched
        // for only inside the grid, not on the boundaries.
        // The searched range can be limited even more (and the
        // code speed improved) by making the eta range smaller
        // with the "minEtaBin" and "maxEtaBin" parameters
        // in the constructor.
        //
        // It is assumed that the phi coordinate is periodic
        // (it has no boundary). The "data" array should be laid
        // out so that the fastest changing coordinate is phi.
        //
        // The peak coordinates (eta and phi) are returned
        // in the units of grid cell size. For example, if the peak
        // is in the middle of the input cell data[5*nPhi + 7] then 
        // the returned peak coordinate will be (5, 7).
        //
        // The function returns 0 in case everything is OK. When peaks
        // with subcell resultion are requested, it returns the number of
        // "unstable" peaks -- that is, peaks whose subcell position are
        // outside of the initial highest cell, the peaks for which the
        // fitted second derivatives are not negative, etc. Significant
        // number of such peaks usually indicates that the round-off errors
        // are too large and that the value of the "peakHeightCutoff"
        // constructor parameter should be increased.
        //
        template <typename Real>
        int find(const Real* data, unsigned nEta, unsigned nPhi,
                 std::vector<Peak>* result);

        // Helper function for fitting the peak position in a more precise
        // manner when curve values are known on a 3x3 grid. The values
        // are fitted to a quadratic polynomial in 2d and the extremum of
        // that polynomial is found. This code assumes that the grid
        // locations are at -1, 0, and 1 in both directions. Correct
        // shifting and scaling is up to the user of this function.
        //
        // The function returns "true" if the extremum found is actually
        // a maximum, "false" otherwise. The extremum coordinates are
        // placed into the locations pointed to by x and y.
        //
        template <typename Real>
        static bool find3by3(const Real gridValues[3][3], Real* x, Real* y,
                             double hessian[3]);

    private:
        PeakFinder();
        PeakFinder& operator=(const PeakFinder&);

        void allocateMaskBuffer(unsigned nEta, unsigned nPhi);
        void clearMaskBuffer(unsigned nEta, unsigned nPhi);
        template <typename Real>
        void fillMask(const Real* data, unsigned nEta, unsigned nPhi,
                      unsigned minEtaBin, unsigned maxEtaBin,
                      unsigned *nseeds, unsigned *npeaks);
        void processMask(unsigned nEta, unsigned nPhi,
                         unsigned minEtaBin, unsigned maxEtaBin);

        bool match3by3(unsigned iEta, unsigned iPhi,
                       unsigned nPhi, const int pattern[3][3]);
        void tag3by3(unsigned iEta, unsigned iPhi,
                     unsigned nPhi, int tagValue);
        template <typename Real>
        void print3by3(const Real gridValues[3][3]) const;

        double peakHeightCutoff_;
        const unsigned minEtaBin_;
        const unsigned maxEtaBin_;
        int *mask;
        unsigned nbins;
        const bool subCellResolution_;
        const bool printFitWarnings_;

        static const int plateau3x3[3][3];
        static const int ridges3x3[4][3][3];
    };
}

#include "fftjet/PeakFinder.icc"

#endif // FFTJET_PEAKFINDER_HH_
