//=========================================================================
// CompositeKernel.hh
//
// A linear combination of kernels. The first member of the pair
// is the factor for a kernel component in the combination. The
// pointer to the kernel component itself is the second member
// of the pair. Note that, if the factors do not sum up to 1,
// the resulting kernel is not likely to be a density.
//
// This class is mainly intended for modeling different smearing
// for charged and neutral parts of a jet (e.g., widening by
// a magnetic field). It is up to the user of this class to make
// sure that the combined result makes sense. In particular,
// use of negative component factors should normally be avoided.
//
// If the "takePointerOwnership" argument is true, the object
// will call the "delete" operator on the component kernel
// pointers in the destructor -- that is, the object will own
// the pointers. You have to decide on the ownership policy
// at the construction time and follow it throughout the use
// of the object.
//
// Note that there is no way to use different bandwidth values
// for the component kernels.
//
// I. Volobouev
// May 2008
//=========================================================================

#ifndef FFTJET_COMPOSITEKERNEL_HH_
#define FFTJET_COMPOSITEKERNEL_HH_

#include <utility>
#include <vector>

#include "fftjet/AbsKernel2d.hh"

namespace fftjet {
    class CompositeKernel :
        public AbsKernel2d,
        public std::vector<std::pair<double,AbsKernel2d*> >
    {
    public:
        // An empty CompositeKernel can be created.
        // The components can be subsequently added
        // by using the "push_back" function.
        explicit CompositeKernel(bool takePointerOwnership=false);
        CompositeKernel(const std::vector<std::pair<double,AbsKernel2d*> >&,
                        bool takePointerOwnership=false);
        virtual ~CompositeKernel();

        void setScaleRatio(double r);
        bool isDensity() const;
        double operator()(double x, double y, double scale) const;
        void supportRectangle( double scale, KernelSupportRectangle *r) const;
        double rectangleAverage(double x, double y, double scale,
                                double dx, double dy) const;

        // The following function assumes that all component
        // kernels are normalized
        void random(double r1, double r2, double scale,
                    double* px, double* py) const;

        inline bool takesPointerOwnership() const {return destroyKernelSet_;}

    private:
        CompositeKernel(const CompositeKernel&);
        CompositeKernel& operator=(const CompositeKernel&);

        const bool destroyKernelSet_;
        mutable std::vector<double> cdf_;
    };
}

#endif // FFTJET_COMPOSITEKERNEL_HH_
