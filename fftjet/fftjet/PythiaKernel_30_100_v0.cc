#include <cassert>

#include "fftjet/PythiaKernel_30_100_v0.hh"
#include "fftjet/pythia_kernel_30_100_v0.scc"

namespace fftjet {
    PythiaKernel_30_100_v0::PythiaKernel_30_100_v0(int scalePower)
        : LogProfileKernel(1.0, 1.0, scalePower,
                           log_profile_vector(), log_profile_range)
    {
    }

    PythiaKernel_30_100_v0::PythiaKernel_30_100_v0(
        double, double, int scalePower, const std::vector<double>& params)
        : LogProfileKernel(1.0, 1.0, scalePower,
                           log_profile_vector(), log_profile_range)
    {
        assert(params.size() == static_cast<unsigned>(nParameters()));
    }    
}
