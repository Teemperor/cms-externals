//=========================================================================
// InterpolatedKernel.hh
//
// Kernel whose data is tabulated on a 2d grid (histogram bin centers).
// In between, the data is interpolated linearly.
//
// I. Volobouev
// August 2008
//=========================================================================

#ifndef FFTJET_INTERPOLATEDKERNEL_HH_
#define FFTJET_INTERPOLATEDKERNEL_HH_

#include <cassert>
#include <iostream>

#include "fftjet/AbsScalableKernel.hh"
#include "fftjet/LinearInterpolator2d.hh"

namespace fftjet {
    class InterpolatedKernel : public AbsScalableKernel
    {
    public:
        template <typename Real>
        inline InterpolatedKernel(double xScaleFactor, double yScaleFactor,
                                  int scalePow, const Real* data,
                                  unsigned nx, double xmin, double xmax,
                                  unsigned ny, double ymin, double ymax)
            : AbsScalableKernel(xScaleFactor, yScaleFactor, scalePow),
              in(data, nx, xmin, xmax, ny, ymin, ymax)
        {
            assert(data);
            in.normalize(1.0);
            isDensity_ = in.isNonNegative();
        }
        virtual ~InterpolatedKernel() {}

        bool operator==(const InterpolatedKernel& r) const;
        inline bool operator!=(const InterpolatedKernel& r) const
            {return !(*this == r);}

        inline bool isDensity() const {return isDensity_;}

        // The "write" function returns "true" if the object
        // is successfully written out.
        bool write(std::ostream& of) const;

        // The following function reads the object in.
        // Returns NULL in case of a problem.
        static InterpolatedKernel* read(std::istream& in);

    private:
        // Constructor for the binary IO
        inline InterpolatedKernel(double xScaleFactor, double yScaleFactor,
                                  int scalePow, const LinearInterpolator2d& i)
            : AbsScalableKernel(xScaleFactor, yScaleFactor, scalePow),
              in(i), isDensity_(in.isNonNegative()) {}

        inline double calculate(double x, double y) const {return in(x, y);}
        void unscaledRectangle(KernelSupportRectangle *r) const;
        void unscaledRandom(double r1, double r2,
                            double* px, double* py) const;

        LinearInterpolator2d in;
        bool isDensity_;

        // Version number for the I/O
        static inline unsigned version() {return 1;}

        // Class id for the I/O
        static inline unsigned classId() {return 4;}
    };
}

#endif // FFTJET_INTERPOLATEDKERNEL_HH_
