//=========================================================================
// Kernel2dFactory.hh
//
// Factories for 2d kernel functions
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_KERNEL2DFACTORY_HH_
#define FFTJET_KERNEL2DFACTORY_HH_

#include <string>
#include <vector>
#include <map>
#include <cassert>

#include "fftjet/AbsKernel2d.hh"

namespace fftjet {
    // Abstract factory for Kernel2d types
    class AbsKernel2dFactory
    {
    public:
        virtual ~AbsKernel2dFactory() {}
        virtual AbsKernel2d* create(double sx, double sy, int scalePower,
                                    const std::vector<double>&) const = 0;
        virtual int nParameters() const = 0;
    };

    // Concrete factory for Kernel2d types
    template <typename T>
    class Kernel2dFactory : public AbsKernel2dFactory
    {
    public:
        virtual ~Kernel2dFactory() {}

        inline int nParameters() const {return T::nParameters();}

        T* create(double sx, double sy, int scalePower,
                  const std::vector<double>& params) const
        {
            if (nParameters() >= 0)
                assert(params.size() == static_cast<unsigned>(nParameters()));
            return new T(sx, sy, scalePower, params);
        }
    };

    // Default kernel factory. It is intended for use inside some
    // interpretive language or testing framework and can iterate
    // over kernel classes. Usage is like this:
    //
    // DefaultKernel2dFactory factory;
    // AbsKernel2d *kernel = factory[userName]->create(sx, sy, p,
    //                                                 userParameters);
    //
    // where "userName" is one of the known kernel names, sx and sy are
    // doubles, p is an integer which determines the dependence of the
    // kernel bandwidth values on the scale, and "userParameters" is
    // a vector of doubles. At some point after this one has to call
    // "delete" operator on the kernel pointer. Known kernel names
    // correspond to the names of classes derived from "AbsKernel2d".
    // For example, to create the Gaussian kernel, one can use
    //
    // AbsKernel2d *kernel = factory["Gauss2d"]->create(1, 1, 1,
    //                                                  std::vector<double>());
    //
    // You can add your own kernel to this factory like this:
    //
    // factory["YourClass"] = new Kernel2dFactory<YourClass>;
    //
    class DefaultKernel2dFactory :
        public std::map<std::string, AbsKernel2dFactory *>
    {
    public:
        DefaultKernel2dFactory();
        virtual ~DefaultKernel2dFactory();

    private:
        DefaultKernel2dFactory(const DefaultKernel2dFactory&);
        DefaultKernel2dFactory& operator=(const DefaultKernel2dFactory&);
    };
}

#endif // FFTJET_KERNEL2DFACTORY_HH_
