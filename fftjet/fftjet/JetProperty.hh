//=========================================================================
// JetProperty.hh
//
// Selectors of properties for Peak/RecombinedJet classes. For use
// with OpenDXTreeFormatter and in other similar contexts.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_JETPROPERTY_HH_
#define FFTJET_JETPROPERTY_HH_

#include <cmath>
#include <cassert>

#include "fftjet/SimpleFunctors.hh"

namespace fftjet {
    // Generic wrapper for all properties which can be
    // expressed as jet member functions with no arguments
    // returning double
    template<class Jet>
    class JetProperty : public Functor1<double,Jet>
    {
    public:
        typedef double (Jet::* JetMemberFunction)() const;

        inline explicit JetProperty(JetMemberFunction fn) : fn_(fn) {}
        inline double operator()(const Jet& j) const {return (j.*fn_)();}

    private:
        JetProperty();
        JetMemberFunction fn_;
    };

    // Generic wrapper for taking logarithms of other properties
    template<class Jet>
    class LogProperty : public Functor1<double,Jet>
    {
    public:
        inline LogProperty(const Functor1<double,Jet>* f, bool own=false)
            : func_(f), ownsPointer_(own) {}
        inline ~LogProperty() {if (ownsPointer_) delete func_;}

        inline double operator()(const Jet& j) const
        {
            const double value = (*func_)(j);
            assert(value > 0.0);
            return log(value);
        }

    private:
        LogProperty();
        const Functor1<double,Jet>* func_;
        const bool ownsPointer_;
    };

    // We need couple extra constants to calculate the Laplacian
    template <class Jet>
    struct MinusScaledLaplacian : public Functor1<double,Jet>
    {
        inline MinusScaledLaplacian(const double bwx, const double bwy)
            : bwx_(bwx), bwy_(bwy) {}

        inline double operator()(const Jet& c) const
        {
            // The following normally makes the Laplacian non-negative
            return -c.scaledLaplacian(bwx_, bwy_);
        }

    private:
        MinusScaledLaplacian();
        const double bwx_;
        const double bwy_;
    };

    // Some potentially useful combinations of jet variables
    template <class Jet>
    struct ScaledHessianDet : public Functor1<double,Jet>
    {
        inline double operator()(const Jet& c) const
        {
            const double s = c.scale();
            const double ssquared = s*s;
            return ssquared*ssquared*c.hessianDeterminant();
        }
    };

    template <class Jet>
    struct ScaledMagnitude : public Functor1<double,Jet>
    {
        inline double operator()(const Jet& c) const
        {
            return c.scale()*c.magnitude();
        }
    };

    template <class Jet>
    struct ScaledMagnitude2 : public Functor1<double,Jet>
    {
        inline double operator()(const Jet& c) const
        {
            const double s = c.scale();
            return s*s*c.magnitude();
        }
    };
}

#endif // FFTJET_JETPROPERTY_HH_
