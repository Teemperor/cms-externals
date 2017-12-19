//=========================================================================
// ProfileKernel.hh
//
// Interpolation kernel using data-defined radial profile
// curve. It is assumed that profile values are calculated
// on the radius interval [0, 1] with uniform step in radius.
// The first value corresponds to r = 0.0, and the last one
// corresponds to 1.0 - 1.0/profile.size(). Each profile value
// must be non-negative.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_PROFILEKERNEL_HH_
#define FFTJET_PROFILEKERNEL_HH_

#include <vector>
#include "fftjet/AbsSymmetricKernel.hh"

namespace fftjet {
    class ProfileKernel : public AbsSymmetricKernel
    {
    public:
        ProfileKernel(double sx, double sy, int scalePower,
                      const std::vector<double>& profile);
        inline virtual ~ProfileKernel() {delete [] cdf_;}

        inline bool isDensity() const {return true;}
        double axialWeight(double r, unsigned nsteps=0) const;

        inline static int nParameters() {return -1;}

    protected:
        void normalizeProfile();
        double rawInterpolatedValue(double r) const;
        double gaussIntegral() const;

    private:
        ProfileKernel(const ProfileKernel&);
        ProfileKernel& operator=(const ProfileKernel&);

        virtual double eval(double rsquared) const;
        virtual double randomRadius(double rnd) const;
        inline virtual double supportDistance() const {return 1.0;}

        void buildCdf();
        unsigned rCellNumber(double r, double *delta) const;

        const std::vector<double> params_;
        double normfactor_;
        double* cdf_;
    };
}

#endif // FFTJET_PROFILEKERNEL_HH_
