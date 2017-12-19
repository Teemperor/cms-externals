//=========================================================================
// Kernels.hh
//
// Simple 2d kernel functions. More complicated kernels are implemented
// in their own separate files.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_KERNELS_HH_
#define FFTJET_KERNELS_HH_

#include <vector>
#include "fftjet/AbsSymmetricKernel.hh"

namespace fftjet {
    // Symmetric beta kernel. Constant kernel can be obtained
    // by using symmetric beta kernel with power 0. Epanechnikov
    // kernel is the symmetric beta kernel with power 1, etc.
    class SymmetricBeta : public AbsSymmetricKernel
    {
    public:
        inline SymmetricBeta(double sx, double sy, int scalePow, double power)
            : AbsSymmetricKernel(sx, sy, scalePow), n_(power) {}
        SymmetricBeta(double sx, double sy, int scalePow,
                      const std::vector<double>& params);
        virtual ~SymmetricBeta() {}

        inline bool isDensity() const {return true;}
        inline double power() const {return n_;}

        inline static int nParameters() {return 1;}

    private:
        double eval(double rsquared) const;
        double randomRadius(double rnd) const;
        inline double supportDistance() const {return 1.0;}
        double calculateResolution() const;

        double n_;
    };

    // Linear kernel. The radial profile descends linearly
    // from a constant value at 0 to zero value at 1.
    class Linear2d : public AbsSymmetricKernel
    {
    public:
        inline Linear2d(double sx, double sy, int scalePow)
            : AbsSymmetricKernel(sx, sy, scalePow) {}
        inline Linear2d(double sx, double sy, int scalePow,
                        const std::vector<double>&)
            : AbsSymmetricKernel(sx, sy, scalePow) {}

        inline bool isDensity() const {return true;}

        inline static int nParameters() {return 0;}

    private:
        double eval(double rsquared) const;
        double randomRadius(double rnd) const;
        inline double supportDistance() const {return 1.0;}
        inline double calculateHalfWeightRadius() const {return 0.5;}
        inline double calculateResolution() const {return 1.0;}
    };

    // Gaussian kernel. Standard deviation is 1 in any direction.
    class Gauss2d : public AbsSymmetricKernel
    {
    public:
        inline Gauss2d(double sx, double sy, int scalePow)
            : AbsSymmetricKernel(sx, sy, scalePow) {}
        inline Gauss2d(double sx, double sy, int scalePow,
                       const std::vector<double>&)
            : AbsSymmetricKernel(sx, sy, scalePow) {}

        inline bool isDensity() const {return true;}
        double axialWeight(double r, unsigned nsteps=0) const;

        inline static int nParameters() {return 0;}

    private:
        double eval(double rsquared) const;
        double randomRadius(double rnd) const;
        double supportDistance() const;
        inline double calculateHalfWeightRadius() const
            {return 1.177410022515474691;}
        inline double calculateResolution() const
            {return 1.414213562373095;}
    };

    // The following functor, together with RealFrequencyKernel class,
    // can be used to calculate stable rotationally symmetric Sub-Gaussian
    // distributions. If you need to know what they are see, for example,
    // section 6 (page 40) of
    // http://academic2.american.edu/~jpnolan/stable/DataAnalysis.ps 
    class SubGauss : public AbsSymmetricKernel
    {
    public:
        SubGauss(double sx, double sy, int scalePow, double alpha);
        SubGauss(double sx, double sy, int scalePow,
                 const std::vector<double>& params);
        virtual ~SubGauss() {}

        inline bool isDensity() const {return true;}
        inline double alpha() const {return alpha_;}

        inline static int nParameters() {return 1;}

    private:
        double eval(double rsquared) const;
        double supportDistance() const;
        double randomRadius(double rnd) const;

        double alpha_;
        double normfactor_;

        static double calculateNormfactor(double alpha);
    };

    // Huber kernel. Weight of the tail must be >= 0 and < 1.
    class Huber2d : public AbsSymmetricKernel
    {
    public:
        Huber2d(double sx, double sy, int scalePow, double tailWeight);
        Huber2d(double sx, double sy, int scalePow,
                const std::vector<double>& params);
        virtual ~Huber2d() {}

        inline bool isDensity() const {return true;}
        inline double tailWeight() const {return tWeight_;}
        inline double transition() const {return a_;}

        inline static int nParameters() {return 1;}

    private:
        double eval(double rsquared) const;
        double randomRadius(double rnd) const;
        double supportDistance() const;
        double calculateHalfWeightRadius() const;

        void setup();
        double weight_(double asquared) const;

        double tWeight_;
        double a_;
        double normfactor_;
    };

    class DeltaFunctionKernel : public AbsKernel2d
    {
    public:
        inline explicit DeltaFunctionKernel(double volume) : value_(volume) {}
        DeltaFunctionKernel(double sx, double sy, int scalePow,
                            const std::vector<double>& params);
        virtual ~DeltaFunctionKernel() {}

        inline void setScaleRatio(double) {}

        inline bool isDensity() const {return true;}
        inline double volume() const {return value_;}

        double operator()(double x, double y, double scale) const;
        void supportRectangle(double scale,
                              KernelSupportRectangle *r) const;
        double rectangleAverage(double x, double y, double scale,
                                double dx, double dy) const;
        void random(double r1, double r2, double scale,
                            double* px, double* py) const;

        inline static int nParameters() {return 1;}

    private:
        double value_;
    };

    // Inverse power kernel proportional to 1/(1 + (r^2)^n).
    // n must be >= 1.05. For n <= 1 the integral diverges,
    // and for n between 1 and 1.05 we encounter various
    // numerical difficulties.
    class InvPower2d : public AbsSymmetricKernel
    {
    public:
        InvPower2d(double sx, double sy, int scalePow, double n);
        InvPower2d(double sx, double sy, int scalePow,
                   const std::vector<double>& params);
        inline virtual ~InvPower2d() {delete [] cdf_;}

        inline bool isDensity() const {return true;}
        inline double power() const {return n_;}
        double axialWeight(double r, unsigned nsteps=0) const;

        inline static int nParameters() {return 1;}

    private:
        InvPower2d(const InvPower2d&);
        InvPower2d& operator=(const InvPower2d&);

        inline double calculateResolution() const {return 1.0;}
        double eval(double rsquared) const;
        double randomRadius(double rnd) const;
        double supportDistance() const;

        void setup();
        void buildCdf();
        double tailIntegral(double r) const;
        double headIntegral(double r) const;
        double annulusIntegral(double middleR, double width) const;
        double calculateHalfWeightRadius() const;

        double n_;
        double normfactor_;
        double lowbound_;
        double lowweight_;
        double upbound_;
        double upweight_;

        double* cdf_;
    };
}

#endif // FFTJET_KERNELS_HH_
