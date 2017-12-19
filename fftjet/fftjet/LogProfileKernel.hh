//=========================================================================
// LogProfileKernel.hh
//
// Interpolation kernel using data-defined radial profile
// curve in the log(density) space. In the log profile,
// the points are assumed to be _in the center_ of uniformly
// spaced bins between 0 and rmax. Beyond rmax the kernel
// is decaying exponentially, with the decay time determined
// from the last two points in the profile.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_LOGPROFILEKERNEL_HH_
#define FFTJET_LOGPROFILEKERNEL_HH_

#include <vector>
#include "fftjet/AbsSymmetricKernel.hh"

namespace fftjet {
    class LogProfileKernel : public AbsSymmetricKernel
    {
    public:
        // The "profile" vector must have at least two points,
        // and the point before last must be larger than the
        // last (to ensure exponential decay at large distances)
        LogProfileKernel(double sx, double sy, int scalePower,
                         const std::vector<double>& profile,
                         double rmax=1.0);
        inline virtual ~LogProfileKernel() {delete [] cdf_;}

        inline bool isDensity() const {return true;}
        double axialWeight(double r, unsigned nsteps=0) const;

        inline static int nParameters() {return -1;}

    private:
        LogProfileKernel(const LogProfileKernel&);
        LogProfileKernel& operator=(const LogProfileKernel&);

        virtual double eval(double rsquared) const;
        virtual double randomRadius(double rnd) const;
        inline double supportDistance() const {return support_;}

        double annulusIntegral(double middleR, double width) const;
        double estimateSupport() const;
        void normalize();
        void buildCdf();

        const std::vector<double> params_;
        const double maxRadius_;
        double normfactor_;
        double weightBeforeMR_;
        double support_;
        double* cdf_;
    };
}

#endif // FFTJET_LOGPROFILEKERNEL_HH_
