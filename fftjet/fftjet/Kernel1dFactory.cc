#include "fftjet/Kernel1dFactory.hh"
#include "fftjet/Kernels1d.hh"

namespace fftjet {
    DefaultKernel1dFactory::DefaultKernel1dFactory()
    {
        (*this)["SymmetricBeta1d"] = new Kernel1dFactory<SymmetricBeta1d>();
        (*this)["Gauss1d"] = new Kernel1dFactory<Gauss1d>();
        (*this)["DeltaFunction1d"] = new Kernel1dFactory<DeltaFunction1d>();
    }

    DefaultKernel1dFactory::~DefaultKernel1dFactory()
    {
        for (std::map<std::string, AbsKernel1dFactory *>::iterator it = begin();
             it != end(); ++it)
            delete it->second;
    }
}
