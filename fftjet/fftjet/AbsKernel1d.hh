//=========================================================================
// AbsKernel1d.hh
//
// Interface class for 1d kernel functions in scale space
//
// I. Volobouev
// November 2009
//=========================================================================

#ifndef FFTJET_ABSKERNEL1D_HH_
#define FFTJET_ABSKERNEL1D_HH_

namespace fftjet {
    // Utility class for kernel-related calculations
    class KernelSupportInterval
    {
    public:
        inline KernelSupportInterval(double minval, double maxval)
            : xmin_(minval), xmax_(maxval) {}

        inline double xmin() const {return xmin_;}
        inline double xmax() const {return xmax_;}

        // Multiplication by a scale factor
        inline KernelSupportInterval operator*(const double& scale) const
        {
            return KernelSupportInterval(xmin_*scale, xmax_*scale);
        }

    private:
        KernelSupportInterval();
        double xmin_;
        double xmax_;
    };

    // Interface class for 1d kernel functions. In all calculations,
    // it should be assumed that the "scale" parameter is positive.
    //
    // Subsequently, derived classes will be used with a templated
    // factory. In order to play nice with that factory, each kernel
    // class must have a constructor which looks like
    //
    // ClassName(double s, int n, const std::vector<double>& v);
    //
    // "s" is typically (but not necessarily) the scale factor
    // which maps the scale argument into bandwidth values. n is normally
    // the power of the scale in the bandwidth calculations: 1 means
    // bandwidth is proportional to the scale, 0 means the bandwidth
    // is independent from the scale, and -1 means that bandwidth is
    // inversely proportional to the scale. Vector "v" is a vector of
    // parameters for the kernels. In addition, each class must have
    // a static function
    //
    // static int nParameters();
    //
    // which defines the number of parameters the class expects to get
    // in "v". If this function returns a negative number, it means that
    // the number of parameters is arbitrary (for example, this would be
    // the case for kernels defined by interpolation tables).
    //
    // Of course, more complicated constructors may be needed for some
    // kernels. Such constructors can not be used with the common factory.
    // These kernels will, consequently, require special parsers in
    // a configuration file or in an interpretive language interface.
    //
    class AbsKernel1d
    {
    public:
        virtual ~AbsKernel1d() {}

        // The function value
        virtual double operator()(double x, double scale) const = 0;

        // The following function should return the smallest interval
        // which bounds the region in which the kernel value is not 0.
        virtual KernelSupportInterval supportInterval(double scale) const = 0;

        // An average over a small interval centered at x
        // with side dx. This is more precise (but slower)
        // than calculating the value of the function at
        // the interval center.
        //
        // Degenerate kernels (e.g., delta function) must override
        // the default implementation of this function.
        virtual double intervalAverage(double x, double sc, double dx) const;

        // The following function tells us whether the kernel
        // can be treated as a probability density. Normally,
        // it means that the kernel must be non-negative everywhere
        // and that its integral over the whole space must be positive.
        virtual bool isDensity() const = 0;

        // The following function returns random numbers distributed
        // according to the kernel treated as a probability density.
        // The argument "rnd" should be a random number between 0 and 1.
        // Calling this function on a kernel which is not a density
        // should result in a run-time error.
        virtual double random(double rnd, double scale) const = 0;

        // The following function scans the kernel for use on intervals.
        // The function returns the area under the scan before normalization.
        // It is assumed that the coordinate of the leftmost point in the
        // grid is xmin + xstep/2, where xstep = (xmax - xmin)/nx.
        // x0 is the position of the kernel center.
        template<typename Real>
        double scanInterval(Real* to, unsigned nx,
                            double x0, double scale,
                            double xmin, double xmax,
                            double normalizationArea=1.0,
                            bool normalize=true) const;

        // The following function scans the kernel for use with
        // periodic topologies. The function returns the
        // area under the scan before normalization.
        template<typename Real>
        double scanCircle(Real* to, unsigned nPhi,
                          double phi0, double scale,
                          double phiBin0Edge,
                          double normalizationArea=1.0,
                          bool normalize=true) const;

        // The following function assumes periodic 1d region
        // with the 2*Pi period. to[0][0] corresponds to kernel
        // argument (0, 0). The function returns the area under
        // the scan before normalization.
        template<typename Real>
        double scanFFT(Real* to, unsigned nPhi,
                       double scale, double normalizationArea=1.0,
                       bool normalize=true) const;
    };
}

#include "fftjet/AbsKernel1d.icc"

#endif // FFTJET_ABSKERNEL1D_HH_
