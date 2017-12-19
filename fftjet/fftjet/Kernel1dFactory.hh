//=========================================================================
// Kernel1dFactory.hh
//
// Factories for 1d kernel functions
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_KERNEL1DFACTORY_HH_
#define FFTJET_KERNEL1DFACTORY_HH_

#include <string>
#include <vector>
#include <map>
#include <cassert>

#include "fftjet/AbsKernel1d.hh"

namespace fftjet {
    // Abstract factory for Kernel1d types
    class AbsKernel1dFactory
    {
    public:
        virtual ~AbsKernel1dFactory() {}
        virtual AbsKernel1d* create(double sx, int scalePower,
                                    const std::vector<double>&) const = 0;
        virtual int nParameters() const = 0;
    };

    // Concrete factory for Kernel1d types
    template <typename T>
    class Kernel1dFactory : public AbsKernel1dFactory
    {
    public:
        virtual ~Kernel1dFactory() {}

        inline int nParameters() const {return T::nParameters();}

        T* create(double sx, int scalePower,
                  const std::vector<double>& params) const
        {
            if (nParameters() >= 0)
                assert(params.size() == static_cast<unsigned>(nParameters()));
            return new T(sx, scalePower, params);
        }
    };

    // Default kernel factory. It is intended for use inside some
    // interpretive language or testing framework and can iterate
    // over kernel classes. Usage is like this:
    //
    // DefaultKernel1dFactory factory;
    // AbsKernel1d *kernel = factory[userName]->create(sx, p, userParameters);
    //
    // where "userName" is one of the known kernel names, sx is
    // double, p is an integer which determines the dependence of the
    // kernel bandwidth values on the scale, and "userParameters" is
    // a vector of doubles. At some point after this one has to call
    // "delete" operator on the kernel pointer. Known kernel names
    // correspond to the names of classes derived from "AbsKernel1d".
    // For example, to create the Gaussian kernel, one can use
    //
    // AbsKernel1d *kernel = factory["Gauss1d"]->create(1, 1,
    //                                                  std::vector<double>());
    //
    // You can add your own kernel to this factory like this:
    //
    // factory["YourClass"] = new Kernel1dFactory<YourClass>();
    //
    class DefaultKernel1dFactory :
        public std::map<std::string, AbsKernel1dFactory *>
    {
    public:
        DefaultKernel1dFactory();
        virtual ~DefaultKernel1dFactory();

    private:
        DefaultKernel1dFactory(const DefaultKernel1dFactory&);
        DefaultKernel1dFactory& operator=(const DefaultKernel1dFactory&);
    };
}

#endif // FFTJET_KERNEL1DFACTORY_HH_
