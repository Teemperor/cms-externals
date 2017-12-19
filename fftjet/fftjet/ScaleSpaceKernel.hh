//=========================================================================
// ScaleSpaceKernel.hh
//
// Generic base class for kernel functions in scale space.
// Application code is not supposed to derive directly from this class.
// Instead, use "AbsKernel2d" or "AbsMembershipFunction" as bases.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_SCALESPACEKERNEL_HH_
#define FFTJET_SCALESPACEKERNEL_HH_

namespace fftjet {
    class AbsKernel2d;
    struct AbsMembershipFunction;

    class ScaleSpaceKernel
    {
    public:
        virtual ~ScaleSpaceKernel() {}

        // The following member sets eta to phi (or x to y) scale ratio
        virtual void setScaleRatio(double r) = 0;

    private:
        friend class AbsKernel2d;
        friend struct AbsMembershipFunction;

        inline ScaleSpaceKernel() {}
    };
}

#endif // FFTJET_SCALESPACEKERNEL_KERNEL
