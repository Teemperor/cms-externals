//=========================================================================
// AbsKernel2d.hh
//
// Interface class for 2d kernel functions in scale space
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_ABSKERNEL2D_HH_
#define FFTJET_ABSKERNEL2D_HH_

#include "fftjet/ScaleSpaceKernel.hh"

namespace fftjet {
    // Utility class for kernel-related calculations
    struct KernelSupportRectangle
    {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
    };

    // Interface class for 2d kernel functions. In all calculations,
    // it should be assumed that the "scale" parameter is positive.
    //
    // Subsequently, derived classes will be used with a templated
    // factory. In order to play nice with that factory, each kernel
    // class must have a constructor which looks like
    //
    // ClassName(double sx, double sy, int n, const std::vector<double>& v);
    //
    // "sx" and "sy" are typically (but not necessarily) scale factors
    // for the x and y directions which map the scale argument into
    // bandwidth values. n is normally the power of the scale in the
    // bandwidth calculations: 1 means bandwidth is proportional to
    // the scale, 0 means the bandwidth is independent from the scale, and
    // -1 means that bandwidth is inversely proportional to the scale.
    // Vector "v" is a vector of parameters for the kernels. In addition,
    // each class must have a static function
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
    class AbsKernel2d : public ScaleSpaceKernel
    {
    public:
        virtual ~AbsKernel2d() {}

        // The following member sets eta to phi (or x to y) scale ratio
        virtual void setScaleRatio(double r) = 0;

        // The function value
        virtual double operator()(double x, double y, double scale) const = 0;

        // The following function should return the smallest rectangle
        // which bounds the region in which the kernel value is not 0.
        //
        // The size of the rectangle can strongly affect the speed
        // of various scanning algorithms. Therefore, resist the
        // temptation to return DBL_MAX for kernels with infinite support
        // and do something more computationally meaningful instead.
        // For example, the Gaussian kernel, while having infinite
        // support theoretically, turns to 0 in double precision
        // representation outside the square with corners at (-39,-39)
        // and (39,39).
        virtual void supportRectangle(double scale,
                                      KernelSupportRectangle *r) const = 0;

        // An average over a small rectangle centered at (x, y)
        // with sides dx and dy. This is more precise (but slower)
        // than calculating the value of the function at the rectangle
        // center.
        //
        // Degenerate kernels (e.g., delta functions in one or two
        // dimensions) must override the default implementation
        // of this function.
        virtual double rectangleAverage(double x, double y, double scale,
                                        double dx, double dy) const;

        // The following function tells us whether the kernel
        // can be treated as a probability density. Normally,
        // it means that the kernel must be non-negative everywhere
        // and that its integral over the whole space must be positive.
        virtual bool isDensity() const = 0;

        // The following function returns random numbers distributed
        // according to the kernel treated as a probability density. The
        // arguments "r1" and "r2" should be random numbers between 0 and 1.
        // Calling this function on a kernel which is not a density
        // should result in a run-time error.
        virtual void random(double r1, double r2, double scale,
                            double* px, double* py) const = 0;

        // The following function scans the kernel for use with planar
        // geometries. The function returns the area under the scan before
        // normalization. It is assumed that the coordinate of the
        // bottom left point in the grid is (xmin + xstep/2, ymin + ystep/2),
        // where xstep = (xmax - xmin)/nx and ystep = (ymax - ymin)/ny.
        template<typename Real>
        double scanPlane(Real* to, unsigned nx, unsigned ny,
                         double x0, double y0, double scale,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         double normalizationArea=1.0,
                         bool normalize=true) const;

        // The following function scans the kernel for use with
        // cylindrical geometries. The function returns the
        // area under the scan before normalization.
        template<typename Real>
        double scanCylinder(Real* to, unsigned nEta, unsigned nPhi,
                            double eta0, double phi0, double scale,
                            double etaMin, double etaMax, double phiBin0Edge,
                            double normalizationArea=1.0,
                            bool normalize=true) const;

        // The following function assumes periodic 2d region
        // with the 2*Pi periods in each direction. to[0][0]
        // corresponds to kernel argument (0, 0). The function
        // returns the area under the scan before normalization.
        template<typename Real>
        double scanFFT(Real* to, unsigned nx, unsigned ny,
                       double scale, double normalizationArea=1.0,
                       bool normalize=true) const;
    private:
        // Various utility functions follow
        template<typename Real>
        long double fillFFTy_(Real* to, unsigned nx, unsigned ny,
                              double x, double scale) const;

        void optimalLineScan(unsigned nx, double xmin, double xmax,
                             double rmin, double rmax,
                             unsigned* nxmin, unsigned* nxmax) const;

        // Clear out the whole scanned area if the following functions
        // return "true". In this case some optimal scan may be
        // available.
        bool optimalPlaneScan(unsigned nx, unsigned ny,
                              double x0, double y0, double scale,
                              double xmin, double xmax,
                              double ymin, double ymax,
                              unsigned* nxmin, unsigned* nxmax,
                              unsigned* nymin, unsigned* nymax) const;
        bool optimalCylinderScan(unsigned nEta, unsigned nPhi,
                                 double eta0, double phi0, double scale,
                                 double etaMin, double etaMax,
                                 double phiBin0Edge,
                                 unsigned* netaMin, unsigned* netaMax,
                                 unsigned* nphiMin, unsigned* nphiMax) const;
    };
}

#include "fftjet/AbsKernel2d.icc"

#endif // FFTJET_ABSKERNEL2D_HH_
