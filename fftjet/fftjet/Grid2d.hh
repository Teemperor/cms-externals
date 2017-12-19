//=========================================================================
// Grid2d.hh
//
// 2d grid (histogram) for deposited energy (Pt, etc.) with cylindrical
// topology. The coverage in phi is 2 Pi.
//
// Note that the data from this grid will be fed directly into
// the FFT engine. This means that the number of bins in this grid
// should be appropriate for FFT and that there should be some
// additional space in eta to avoid energy leakage between high
// and low eta regions during convolution. It is better to add this
// additional space on both sides, as the peaks will not be found
// in the first and the last eta bin.
//
// This class is used inside some pretty tight loops. Don't create
// any virtual functions in it.
//
// I. Volobouev
// January 2008
//=========================================================================

#ifndef FFTJET_GRID2D_HH_
#define FFTJET_GRID2D_HH_

#include <string>
#include <iostream>

namespace fftjet {
    class StatAccumulator;

    template<typename Real>
    class Grid2d
    {
    public:
        // Constructors
        Grid2d(unsigned nEtaBins, Real etaMin, Real etaMax,
               unsigned nPhiBins, Real phiBin0Edge,
               const char* title = "");
        Grid2d(const Grid2d&);

        // Destructor
        ~Grid2d();

        // Assignment operator
        Grid2d& operator=(const Grid2d&);

        // Change title (could be useful after copying)
        void setTitle(const char* newtitle);

        // Comparison for equality
        bool operator==(const Grid2d& r) const;
        bool operator!=(const Grid2d& r) const;

        // Basic accessors
        unsigned nEta() const;
        unsigned nPhi() const;
        Real etaMin() const;
        Real etaMax() const;
        Real phiBin0Edge() const;
        const char* title() const;

        // Bin width quantities
        Real etaBinWidth() const;
        Real phiBinWidth() const;

        // Translators from coordinates to bin numbers.
        // Eta bin number can be negative or very large
        // which means that eta is out of histogram range.
        int getEtaBin(Real eta) const;
        unsigned getPhiBin(Real phi) const;

        // Translators from bin numbers to coordinates.
        Real etaBinCenter(int etaBin) const;
        Real phiBinCenter(unsigned phiBin) const;

        // Translators from real-valued bin numbers to coordinates.
        // These function perform the linear mapping from the coordinates
        // in the units of bin width to the normal coordinates.
        Real etaFromBinNum(Real etaBin) const;
        Real phiFromBinNum(Real phiBin) const;

        // The following function checks grid compatibility
        // for various calculations
        bool isCompatible(const Grid2d& g) const;

        // The following function checks whether all data
        // enstries are non-negative
        bool isNonNegative() const;

        // Finding minimum/maximum of the data values
        void findMinimum(unsigned* etaBin, unsigned* phiBin, Real* min) const;
        void findMaximum(unsigned* etaBin, unsigned* phiBin, Real* max) const;

        // Pointer to the data. To avoid any future incompatibilities,
        // do not manipulate the data "by hand". FFT engine needs
        // the direct data access for speed. It knows how the data
        // is laid out internally (and you are not supposed to).
        const Real* data() const;

        // Access to data values
        Real binValue(unsigned etaBin, unsigned phiBin) const;
        Real coordValue(Real eta, Real phi) const;

        // Values using unchecked array bounds
        Real uncheckedAt(unsigned etaBin, unsigned phiBin) const;

        // Sum of all data bins
        Real sum() const;

        // Integral of the data (sum of the data bins times bin area)
        Real integral() const;

        // Number of grid cells above some threshold value
        unsigned countAboveThreshold(Real threshold) const;

        // Integral of the data above some threshold value
        Real integralAboveThreshold(Real threshold) const;

        // Fuzzy integral of the data using the given threshold value.
        // Each datum is partitioned between the "signal" integral
        // (returned by this function) and "background" integral
        // according to the weights determined by datum/(threshold+datum)
        // and threshold/(threshold+datum), respectively. Note that this
        // partitioning makes sense only when both datum and threshold
        // are non-negative. The uncertainty of the integral is calculated
        // assuming that all cells are independent and that each cell
        // belongs completely to either signal or background.
        void fuzzyIntegral(Real thresh, Real *value, Real *uncertainty) const;

        // Accumulate statistics over the grid values (just the values,
        // not the coordinates). The accumulator will not be reset, so
        // the statistics can be accumulated from several grids (or from
        // multiple events).
        void accumulateDataStats(StatAccumulator* accumulator) const;

        // Fill function by coordinate. Overflows are ignored.
        void fill(Real eta, Real phi, Real energy);

        // The standard "fill" function fills the nearest four bins
        // in such a way that the center of mass of these fills
        // coincides with the coordinate given. Instead, the functions
        // below fill or set grid bins as in a histogram. Use these
        // functions (they are much faster) if you do not care about
        // loss of resolution due to binning.
        void fillFast(Real eta, Real phi, Real energy);
        void setFast(Real eta, Real phi, Real value);

        // Fill/set functions by bin number. Bin number
        // out of range will result in an assertion error.
        void fillBin(unsigned etaBin, unsigned phiBin, Real energy);
        void setBin(unsigned etaBin, unsigned phiBin, Real value);

        // Unchecked fill/set functions by bin number.
        // Can be used to improve the code speed in case
        // you are sure about your array bounds.
        void uncheckedFillBin(unsigned etaBin, unsigned phiBin, Real energy);
        void uncheckedSetBin(unsigned etaBin, unsigned phiBin, Real value);

        // Set everything. Assume that the data layout is correct.
        void blockSet(const Real *from, unsigned nEtaFrom, unsigned nPhiFrom);

        // Fill everything (add to internal data).
        // This function assumes that the data layout is correct.
        void blockFill(const Real *from, unsigned nEtaFrom, unsigned nPhiFrom);

        // Fill everything from another grid. The binning of the other
        // grid must be compatible
        Grid2d& operator+=(const Grid2d& r);

        // Subtract another grid. The binning of the other grid
        // must be compatible.
        Grid2d& operator-=(const Grid2d& r);

        // The following function sets all data values to the given constant
        void reset(Real value=0);

        // The following function multiplies all data values
        // by the given constant
        void scaleData(Real scaleFactor);

        // The following function multiplies data values by an eta-dependent
        // scale factor
        void scaleData(const Real* scaleFactorArray, unsigned arrayLength);

        // The following function sets all data below or equal to 
        // the given threshold to 0
        void applyThreshold(Real thresh);

        // The following function sets to 0 all data below or equal to
        // the ratio of the given threshold to cosh(eta). This is useful
        // for suppressing noise whose energy (not Et) distribution is
        // independent from eta while the grid itself is filled using Et.
        void applyCoshThreshold(Real thresh);

        // The following function sets all data below or equal to 
        // the given threshold to 0 and above the given threshold to 1
        void applyLogicalThreshold(Real thresh);

        // It may be interesting to cluster pow(Et, alpha) instead of
        // just Et. Apply the "power" function after filling the grid
        // to modify the data accordingly. Note that the "alpha" argument
        // cannot be negative, and that all grid entries less or equal
        // to 0 are set to 0 (even if "alpha" is 0).
        //
        // Doing this should be considered as an attempt to modify
        // the metric of the angular space (energies in real calorimeter
        // cells are added linearly no matter what we do here). One should
        // keep in mind that the occupancy away from the jet axis decays
        // approximately as r^-2 and that the average energy scales roughly
        // as r^-1, so that the angular energy profile decays as r^-3.
        // This means that, for example, using the "power" function with
        // argument 2 will cause the angular Et^2 profile to decay as r^-4.
        // The same asymptotic behavior of the Et profile can be achieved
        // by working with the variable r' = pow(r, k), k = 1/2. This is
        // because the occupancy behaves as r'^-2 no matter what k is 
        // (logarithmic infrared divergence in the number of soft particles),
        // while average energy now decays as r'^(-1/k). In a similar vein,
        // setting "alpha" to 1/2 results in k = 2.
        //
        // It is not clear how important these arguments are in the presence
        // of calorimeter thresholds.
        void power(double alpha);

        // The following function clears some eta range.
        // Can be useful in order to prepare the data for FFT.
        // Bin number which corresponds to minEta is cleared
        // while bin number which corresponds to maxEta
        // is not cleared (but the last bin is cleared
        // in case maxEta is beyound the grid range).
        void clearEtaRange(Real minEta, Real maxEta);

        // Same thing but indexing by bin number. Bin number
        // etaBinMin is cleared and etaBinMax is not.
        void clearEtaBins(unsigned etaBinMin, unsigned etaBinMax);

        // The "write" function returns "true" if the object
        // is successfully written out.
        bool write(std::ostream& of) const;

        // The following function reads the object in.
        // Returns NULL in case of a problem.
        static Grid2d* read(std::istream& in);

    private:
        Grid2d();

        void fillTwo(unsigned etaBin, Real phi, Real energy);

        Real* data_;
        std::string title_;
        Real phiBinWidth_;
        Real etaMin_;
        Real etaMax_;
        Real etaRange_;
        Real phiBin0Edge_;
        unsigned nEta_;
        unsigned nPhi_;

        // Version number for the I/O
        static inline unsigned version() {return 1;}

        // Class id for the I/O
        static inline unsigned classId() {return 7;}
    };
}

#include "fftjet/Grid2d.icc"

#endif // FFTJET_GRID2D_HH_
