//=========================================================================
// InterpolatedKernel3d.hh
//
// Kernel whose data is tabulated on a 3d grid. In between, the data
// is interpolated linearly. For scales below the range of scales covered,
// the slice with the smallest scale is used. For scales above the range
// the slice with the largest scale is used.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_INTERPOLATEDKERNEL3D_HH_
#define FFTJET_INTERPOLATEDKERNEL3D_HH_

#include <vector>
#include <iostream>

#include "fftjet/AbsKernel2d.hh"
#include "fftjet/LinearInterpolator2d.hh"

namespace fftjet {
    class InterpolatedKernel3d : public AbsKernel2d
    {
    public:
        // The "scales" vector in the constructor must represent
        // the scales arranged in the increasing order
        // (not necessarily equidistant). The "data" array
        // can be NULL. In that case it is expected that the
        // data for each scale will be provided by calling the
        // "setScaleData" function. If the "data" array in the
        // constructor is not NULL, the number of elements there
        // should be at least scales.size()*nx*ny. It is assumed
        // that the array index which corresponds to the y variable
        // changes most often and that the index which corresponds
        // to the scale variable changes least often.
        //
        // The "useLogSpaceForScale" argument specifies whether
        // the code should use scale or log(scale) as the variable
        // in which the interpolation is linear.
        //
        // "nx", "xmin", and "xmax" arguments specify the binning
        // in the x direction. It is assumed that the first bin 
        // is centered at xmin + (xmax - xmin)/(2 nx), just like
        // in a histogram. Same thing for y.
        //
        template <typename Real>
        InterpolatedKernel3d(const Real* data,
                             const std::vector<double>& scales,
                             bool useLogSpaceForScale,
                             unsigned nx, double xmin, double xmax,
                             unsigned ny, double ymin, double ymax);
        virtual ~InterpolatedKernel3d() {}

        bool operator==(const InterpolatedKernel3d& r) const;
        inline bool operator!=(const InterpolatedKernel3d& r) const
            {return !(*this == r);}

        inline double xMin() const {return xmin_;}
        inline double xMax() const {return xmax_;}
        inline double yMin() const {return ymin_;}
        inline double yMax() const {return ymax_;}
        inline unsigned nx() const {return nxpoints_;}
        inline unsigned ny() const {return nypoints_;}
        inline bool usesLogSpace() const {return useLogSpace_;}
        inline const std::vector<double>& scales() const {return scales_;}

        template <typename Real>
        void setScaleData(unsigned scaleBinNumber, const Real* data);

        // AbsKernel2d functions which must be implemented
        void setScaleRatio(double r);
        bool isDensity() const;
        double operator()(double x, double y, double scale) const;
        void supportRectangle(double scale,
                              KernelSupportRectangle *r) const;
        void random(double r1, double r2, double scale,
                            double* px, double* py) const;

        // The "write" function returns "true" if the object
        // is successfully written out.
        bool write(std::ostream& of) const;

        // The following function reads the object in.
        // Returns NULL in case of a problem.
        static InterpolatedKernel3d* read(std::istream& in);

    private:
        InterpolatedKernel3d();
        void buildRandomizer(unsigned lowerBound, double scale) const;

        const double xmin_;
        const double xmax_;
        const double ymin_;
        const double ymax_;
        const unsigned nxpoints_;
        const unsigned nypoints_;
        const bool useLogSpace_;
        const std::vector<double> scales_;

        std::vector<LinearInterpolator2d> interpols_;

        mutable LinearInterpolator2d randomizer_;
        mutable double randomizerScale_;

        // isDensity_: 1 - yes, 0 - unknown, -1 - no
        mutable int isDensity_;

        // Version number for the I/O
        static inline unsigned version() {return 1;}

        // Class id for the I/O
        static inline unsigned classId() {return 3;}
    };
}

#include "fftjet/InterpolatedKernel3d.icc"

#endif // FFTJET_INTERPOLATEDKERNEL3D_HH_
