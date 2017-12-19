//=========================================================================
// PhiKernels.hh
//
// Kernel functions in phi
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_PHIKERNELS_HH_
#define FFTJET_PHIKERNELS_HH_

#include <vector>

#include "fftjet/AbsScalablePhiKernel.hh"

namespace fftjet {
    class PhiGauss : public AbsScalablePhiKernel
    {
    public:
        inline PhiGauss(const double sy, const int scalePow)
            : AbsScalablePhiKernel(sy, scalePow) {}
        inline PhiGauss(double, double sy, int scalePow,
                        const std::vector<double>&)
            : AbsScalablePhiKernel(sy, scalePow) {}

        inline bool isDensity() const {return true;}

        inline static int nParameters() {return 0;}

    private:
        double unscaledPhiFcn(double phi) const;
        void unscaledPhiSupport(double *phimin, double *phimax) const;
        double unscaledPhiRandom(double rnd) const;
    };

    // Interpolation kernel for smearing in the phi direction only.
    // It is assumed that the kernel is symmetric w.r.t. the change of
    // sign of the phi angle, and the values are calculated on the phi
    // interval [0, Pi/2] with uniform step in phi. The first value
    // corresponds to r = 0.0, and the last one corresponds to
    // Pi/2*(1.0 - 1.0/profile.size()). Each profile value must be
    // non-negative.
    class PhiProfileKernel : public AbsScalablePhiKernel
    {
    public:
        PhiProfileKernel(double sx, double sy, int scalePower,
                         const std::vector<double>& profile);
        inline virtual ~PhiProfileKernel() {delete [] cdf_;}

        inline bool isDensity() const {return true;}

        inline static int nParameters() {return -1;}

    private:
        PhiProfileKernel(const PhiProfileKernel&);
        PhiProfileKernel& operator=(const PhiProfileKernel&);

        double unscaledPhiFcn(double phi) const;
        void unscaledPhiSupport(double *phimin, double *phimax) const;
        double unscaledPhiRandom(double rnd) const;
        unsigned rCellNumber(const double r, double *diff) const;
        double randomDistance(double rnd) const;

        const std::vector<double> params_;
        double* cdf_;
        double normfactor_;
    };
}

#endif // FFTJET_PHIKERNELS_HH_
