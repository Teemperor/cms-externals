//=========================================================================
// PythiaKernel_30_100_v0.hh
//
// Angular energy distribution for jets produced by light quarks
// with momenta between 30 and 100 GeV/c, generated using Pythia
// jet gun. The natural scale for this kernel is 1.0/(parton Pt),
// (if scalePower == 1) where parton Pt is in units GeV/c.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_PYTHIAKERNEL_30_100_V0_HH_
#define FFTJET_PYTHIAKERNEL_30_100_V0_HH_

#include "fftjet/LogProfileKernel.hh"

namespace fftjet {
    class PythiaKernel_30_100_v0 : public LogProfileKernel
    {
    public:
        explicit PythiaKernel_30_100_v0(int scalePower = 1);
        PythiaKernel_30_100_v0(double sx, double sy, int scalePower,
                               const std::vector<double>& params);

        inline static int nParameters() {return 0;}
    };
}

#endif // FFTJET_PYTHIAKERNEL_30_100_V0_HH_
