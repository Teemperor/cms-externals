//=========================================================================
// EtaPhiEtInterpolator.hh
//
// A class for an efficient storage of eta-phi-energy histograms. It is
// assumed that for most eta-phi values only a small number of the lowest
// energy bins are populated. This class is useful for representing
// jet shapes.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_ETAPHIETINTERPOLATOR_HH_
#define FFTJET_ETAPHIETINTERPOLATOR_HH_

#include <string>
#include <iostream>

namespace fftjet {
    template<typename Real>
    class EtaPhiEtInterpolator
    {
    public:
        // Constructor from a regularly spaced data.
        // Note that we will be using it on some noisy
        // data (produced by FFT). Therefore, we need
        // to ignore bins which contribute very little
        // to the overall energy. The "energyFractionToIgnore"
        // argument tells us how much energy we are allowed
        // to ignore. The size of the "data" array should be
        // at least nEtaBins*nEtaBins*nEnergyBins, and the
        // ordering of the dimensions should be [eta][phi][energy].
        template <typename InputReal>
        EtaPhiEtInterpolator(const InputReal* data,
                     Real energyFractionToIgnore,
                     unsigned nEtaBins, Real etaMin, Real etaMax,
                     unsigned nPhiBins, Real phiBin0Edge,
                     unsigned nEnergyBins, Real energyMin, Real energyMax,
                     const char* title = "");

        // Copy constructor
        EtaPhiEtInterpolator(const EtaPhiEtInterpolator&);

        // Destructor
        ~EtaPhiEtInterpolator();

        // Assignment operator
        EtaPhiEtInterpolator& operator=(const EtaPhiEtInterpolator&);

        // Comparison (useful for testing)
        bool operator==(const EtaPhiEtInterpolator& r) const;
        inline bool operator!=(const EtaPhiEtInterpolator& r) const
            {return !(*this == r);}

        // Change title (could be useful after copying)
        void setTitle(const char* newtitle);

        // Basic accessors
        unsigned nEta() const;
        unsigned nPhi() const;
        unsigned nEnergy() const;
        Real ignoredFraction() const;
        Real etaMin() const;
        Real etaMax() const;
        Real phiBin0Edge() const;
        Real energyMin() const;
        Real energyMax() const;
        const char* title() const;

        // The 3-d bin size
        Real binSize() const;

        // The length of the zero-suppressed data buffer
        unsigned zeroSuppressedBufLen() const;

        // Integral of (bin contents times by energy),
        // after removing the ignored fraction
        Real totalEnergy() const;

        // Integral of all data bins (after removing the ignored fraction)
        Real multiplicity() const;

        // Data lookup functions
        Real binValue(unsigned etaBin, unsigned phiBin, unsigned eBin) const;
        Real uncheckedAt(unsigned etaBin, unsigned phiBin, unsigned eBin) const;

        // Linearly interpolated value of the energy density
        Real operator()(Real eta, Real phi, Real et) const;

        // Maximum energy value for which the density
        // is not set to 0. The value will be linearly
        // interpolated from the four closest energy bins.
        Real maxNon0Energy(Real eta, Real phi) const;

        // Bin number lookups. Eta and energy bin numbers
        // can be negative or very large which means that
        // the input argument is out of histogram range.
        int getEtaBin(Real eta) const;
        unsigned getPhiBin(Real phi) const;
        int getEnergyBin(Real et) const;

        // Translators from bin numbers to coordinates
        Real etaBinCenter(int etaBin) const;
        Real phiBinCenter(unsigned phiBin) const;
        Real energyBinCenter(int eBin) const;

        // The following function multiplies all data values
        // by the given constant
        void scaleData(Real scaleFactor);

        // The "write" function returns "true" if the object
        // is successfully written out.
        bool write(std::ostream& of) const;

        // The following function reads the object in.
        // Returns NULL in case of a problem.
        static EtaPhiEtInterpolator* read(std::istream& in);

    private:
        EtaPhiEtInterpolator();

        // Constructor from already sparsified data 
        // (for use in read/write procedures). Assumes
        // ownership of "sparseData" and "binBounds".
        EtaPhiEtInterpolator(Real* sparseData,
                     unsigned* binBounds,
                     Real energyFractionToIgnore,
                     Real totalEnergy, Real totalCount,
                     unsigned nEtaBins, Real etaMin, Real etaMax,
                     unsigned nPhiBins, Real phiBin0Edge,
                     unsigned nEnergyBins, Real energyMin, Real energyMax,
                     const char* title);

        Real* data_;
        unsigned* binBounds_;
        std::string title_;
        Real ignoredFraction_;
        Real totalEnergy_;
        Real totalCount_;
        Real phiBinWidth_;
        Real etaMin_;
        Real etaMax_;
        Real etaRange_;
        Real phiBin0Edge_;
        Real energyMin_;
        Real energyMax_;
        Real energyRange_;
        unsigned nEta_;
        unsigned nPhi_;
        unsigned nEnergy_;

        template <typename InputReal>
        static InputReal determineCutoff(const InputReal* data,
                                         unsigned angleBins, unsigned eBins,
                                         Real energyMin, Real energyMax,
                                         Real dataFractionToIgnore);
        // Version number for the I/O
        static inline unsigned version() {return 1;}

        // Class id for the I/O
        static inline unsigned classId() {return 1;}
    };
}

#include "fftjet/EtaPhiEtInterpolator.icc"

#endif // FFTJET_ETAPHIETINTERPOLATOR_HH_
