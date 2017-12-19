//=========================================================================
// MagneticSmearingKernel.hh
//
// This function represents smearing of a jet in phi due to the
// presence of magnetic field parallel to the beam axis, and it is
// really a 1-d function (no smearing in eta).
//
// The scale for this kernel should be set to 1/(jet Pt).
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_MAGNETICSMEARINGKERNEL_HH_
#define FFTJET_MAGNETICSMEARINGKERNEL_HH_

#include "fftjet/AbsPhiKernel.hh"
#include "fftjet/LinearInterpolator1d.hh"

namespace fftjet {
    // We will need the jet fragmentation function in order to build
    // the kernel. We will assume that this function is implemented
    // as a functor which has "double operator()(double x) const".
    // We will also assume that the lifetime of the fragmentation
    // function is longer than the lifetime of the kernel.
    //
    template <class FragmentationFunction>
    class MagneticSmearingKernel : public AbsPhiKernel
    {
    public:
        // Charged particle deviation in phi due to the presence
        // of magnetic field can be calculated as follows:
        //
        //            numeratorConst
        // sin(phi) = --------------
        //                  Pt
        //
        // Where "numeratorConst" depends on the units in which Pt,
        // magnetic field, and calorimeter radius are measured.
        // The "nominal" value of this constant is determined from
        // the relationship
        //
        // sin(phi) = R/(2*a)
        // a = Pt*1.0e9/(c*B), so that numeratorConst = c*B*R/2.0e9
        //
        // where
        //
        // Pt -- transverse particle momentum, in GeV/c
        // R  -- "typical" calorimeter radius, in m
        // a  -- gyration radius, in m, for a particle with given Pt
        //         and charge equal in magnitude to the charge of
        //         the electron (that is, |Z| == 1).
        // c  -- speed of light in m/s (c = 299792458.0)
        // B  -- magnetic field, in Tesla
        //
        // Of course, if particle Pt is below "numeratorConst" then
        // it simply never reaches the calorimeter.
        //
        // The sign of phi depends on the sign of charge and direction
        // of the magnetic field. We will assume that the number of
        // negatively charged particles is the same, on average, as
        // the number of positively charged particles, and that the
        // smearing is symmetric in phi. We will further assume that
        // the fragmentation function for particles with |Z| == 1
        // is the same as for |Z| == 2.
        //
        // The constructor arguments are as follows:
        //
        // fcn             -- the fragmentation function functor
        //
        // numeratorConst  -- explained above
        //
        // charge1Fraction -- average energy fraction deposited by
        //                    particles with |Z| == 1
        //
        // charge2Fraction -- average energy fraction deposited by
        //                    particles with |Z| == 2
        //
        // samplesPerBin   -- number of points per each phi bin
        //                    in which the fragmentation function
        //                    gets evaluated (and then integrated over)
        //
        // objectOwnsFcn   -- if "true", the fragmentation function functor
        //                    will be deleted in the destructor
        //
        // It is assumed that the energy fraction deposited by neutral
        // particles is (1 - charge1Fraction - charge2Fraction).
        //
        MagneticSmearingKernel(const FragmentationFunction* fcn,
                            double numeratorConst, double charge1Fraction,
                            double charge2Fraction, unsigned samplesPerBin=100,
                            bool objectOwnsFcn=false);
        virtual ~MagneticSmearingKernel();

        inline double numeratorConst() const {return numeratorConst_;}
        inline double charge1Fraction() const {return charge1Fraction_;}
        inline double charge2Fraction() const {return charge2Fraction_;}
        inline unsigned samplesPerBin() const {return samplesPerBin_;}

        inline bool isDensity() const {return true;}
        double rectangleAverage(double x, double y, double scale,
                                double dx, double dy) const;
    private:
        MagneticSmearingKernel();
        MagneticSmearingKernel(const MagneticSmearingKernel&);
        MagneticSmearingKernel& operator=(const MagneticSmearingKernel&);

        double phiFcn(double phi, double scale) const;
        void phiSupport(double scale,
                        double *phimin, double *phimax) const;
        double phiRandom(double rnd, double scale) const;

        double deltaFunAverage(double x, double y,
                               double dx, double dy) const;
        void buildRandomizer();

        const FragmentationFunction* fcn_;
        const double numeratorConst_;
        const double charge0Fraction_;
        const double charge1Fraction_;
        const double charge2Fraction_;
        const unsigned samplesPerBin_;
        const bool objectOwnsFcn_;

        LinearInterpolator1d* randomizer_;
    };
}

#include "fftjet/MagneticSmearingKernel.icc"

#endif // FFTJET_MAGNETICSMEARINGKERNEL_HH_
