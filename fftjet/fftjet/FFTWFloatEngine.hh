//=========================================================================
// FFTWFloatEngine.hh
//
// Single precision FFT engine for jet reconstruction which uses
// the FFTW library from http://www.fftw.org
//
// I. Volobouev
// January 2008
//=========================================================================

#ifndef FFTJET_FFTWFLOATENGINE_HH_
#define FFTJET_FFTWFLOATENGINE_HH_

#include "fftjet/FFTWEngine.hh"
#include "fftw3.h"

namespace fftjet {
    class FFTWFloatEngine : public FFTWEngine<float, fftwf_complex, fftwf_plan>
    {
    public:
        FFTWFloatEngine(unsigned nEta, unsigned nPhi,
                        unsigned optimization = FFTW_PATIENT);
        virtual ~FFTWFloatEngine();
    };
}

#include "fftjet/FFTWFloatEngine.icc"

#endif // FFTJET_FFTWFLOATENGINE_HH_
