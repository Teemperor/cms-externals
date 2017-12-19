//=========================================================================
// InterpolatedMembershipFcn.hh
//
// Jet membership function whose data is tabulated on a 4-d grid
// (partially zero suppressed). In between, the data is interpolated
// linearly. For scales below the range of scales covered, the slice
// with the smallest scale is used. For scales above the range the slice
// with the largest scale is used.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_INTERPOLATEDMEMBERSHIPFCN_HH_
#define FFTJET_INTERPOLATEDMEMBERSHIPFCN_HH_

#include <vector>
#include <iostream>

#include "fftjet/AbsMembershipFunction.hh"
#include "fftjet/EtaPhiEtInterpolator.hh"

namespace fftjet {
    template<typename Real>
    class InterpolatedMembershipFcn : public AbsMembershipFunction
    {
    public:
        // The "scales" vector in the constructor must represent
        // the scales arranged in the increasing order (not
        // necessarily equidistant).
        //
        // The "useLogSpaceForScale" argument specifies whether
        // the code should use scale or log(scale) as the variable
        // in which the interpolation is linear.
        //
        // The "energyFractionToIgnore" argument specifies the
        // energy fraction of the whole distribution which can be
        // ignored in order to improve the compression. Numbers
        // about 1.0e-3 or 1.0e-4 seem to work well.
        //
        // "nEtaBins", "etaMin", and "etaMax" arguments specify
        // the binning in the eta direction. It is assumed that
        // the first bin is centered at
        // etaMin + (etaMax - etaMin)/(2 nEtaBins), just like
        // in a histogram. Same thing for the energy variable.
        //
        // The phi variable is assumed to have 2 Pi range,
        // with the first bin edge at "phiBin0Edge".
        //
        InterpolatedMembershipFcn(
            const std::vector<double>& scales,
            bool useLogSpaceForScale,
            Real energyFractionToIgnore,
            unsigned nEtaBins, Real etaMin, Real etaMax,
            unsigned nPhiBins, Real phiBin0Edge,
            unsigned nEnergyBins, Real energyMin, Real energyMax);

        virtual ~InterpolatedMembershipFcn();

        // Comparison (useful for testing)
        bool operator==(const InterpolatedMembershipFcn& r) const;
        inline bool operator!=(const InterpolatedMembershipFcn& r) const
            {return !(*this == r);}

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
        bool usesLogSpace() const;
        const std::vector<double>& scales() const;

        // Set the data for one of the scales. The "data" array should
        // have nEtaBins*nPhiBins*nEnergyBins elements. It is assumed
        // that the array index associated with the energy variable
        // changes most often, while the array index associated with the
        // eta variable changes least often.
        template <typename InputReal>
        void setScaleData(unsigned scaleBinNumber, const InputReal* data,
                          double dataScaleFactor=1.0);

        // The following function must be overriden from the base
        inline void setScaleRatio(double) {}

        // The following function returns "true" if the data were provided
        // for all the scales
        inline bool isReady() const {return ready_;}

        // Info for one of the scales
        Real totalEnergy(unsigned scaleBin) const;
        Real multiplicity(unsigned scaleBin) const;

        // AbsMembershipFunction methods to implement
        double operator()(double eta, double phi,
                          double et, double scale) const;
        double absorbableEnergy(double eta, double phi,
                                double scale) const;

        // The "write" function returns "true" if the object
        // is successfully written out.
        bool write(std::ostream& of) const;

        // The following function reads the object in.
        // Returns NULL in case of a problem.
        static InterpolatedMembershipFcn* read(std::istream& in);

    private:
        InterpolatedMembershipFcn();
        InterpolatedMembershipFcn(const InterpolatedMembershipFcn&);
        InterpolatedMembershipFcn& operator=(const InterpolatedMembershipFcn&);

        const std::vector<double> scales_;
        const Real ignoredFraction_;
        const Real etaMin_;
        const Real etaMax_;
        const Real phiBin0Edge_;
        const Real energyMin_;
        const Real energyMax_;
        const unsigned nEta_;
        const unsigned nPhi_;
        const unsigned nEnergy_;
        const bool useLogSpace_;

        bool ready_;
        std::vector<EtaPhiEtInterpolator<Real>*> interpols_;

        // The "setScale" function for use with I/O operators.
        // This class will assume ownership of the "interp" object.
        void setScale(unsigned scaleBin, EtaPhiEtInterpolator<Real>* interp);

        // Version number for the I/O
        static inline unsigned version() {return 1;}

        // Class id for the I/O
        static inline unsigned classId() {return 2;}
    };
}

#include "fftjet/InterpolatedMembershipFcn.icc"

#endif // FFTJET_INTERPOLATEDMEMBERSHIPFCN_HH_
