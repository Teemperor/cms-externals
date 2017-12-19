#include "fftjet/Kernel2dFactory.hh"
#include "fftjet/Kernels.hh"
#include "fftjet/ProfileKernel.hh"
#include "fftjet/LogProfileKernel.hh"
#include "fftjet/PythiaKernel_30_100_v0.hh"
#include "fftjet/PhiKernels.hh"

namespace fftjet {
    DefaultKernel2dFactory::DefaultKernel2dFactory()
    {
        // "Normal" kernels
        (*this)["SymmetricBeta"] = new Kernel2dFactory<SymmetricBeta>;
        (*this)["Gauss2d"] = new Kernel2dFactory<Gauss2d>;
        (*this)["SubGauss"] = new Kernel2dFactory<SubGauss>;
        (*this)["Huber2d"] = new Kernel2dFactory<Huber2d>;
        (*this)["Linear2d"] = new Kernel2dFactory<Linear2d>;
        (*this)["InvPower2d"] = new Kernel2dFactory<InvPower2d>;
        (*this)["ProfileKernel"] = new Kernel2dFactory<ProfileKernel>;
        (*this)["LogProfileKernel"] = new Kernel2dFactory<LogProfileKernel>;
        (*this)["PythiaKernel_30_100_v0"] = new Kernel2dFactory<PythiaKernel_30_100_v0>;

        // Phi kernels
        (*this)["PhiGauss"] = new Kernel2dFactory<PhiGauss>;
        (*this)["PhiProfileKernel"] = new Kernel2dFactory<PhiProfileKernel>;

        // Delta function
        (*this)["DeltaFunctionKernel"] = new Kernel2dFactory<DeltaFunctionKernel>;
    }

    DefaultKernel2dFactory::~DefaultKernel2dFactory()
    {
        for (std::map<std::string, AbsKernel2dFactory *>::iterator it = begin();
             it != end(); ++it)
            delete it->second;
    }
}
