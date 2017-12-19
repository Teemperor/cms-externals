//=========================================================================
// FFTWDoubleEngine.hh
//
// Double precision FFT engine for jet reconstruction which uses
// the FFTW library from http://www.fftw.org
//
// I. Volobouev
// January 2008
//=========================================================================

#ifndef FFTJET_FFTWDOUBLEENGINE_HH_
#define FFTJET_FFTWDOUBLEENGINE_HH_

#include "fftjet/FFTWEngine.hh"
#include "fftw3.h"

namespace fftjet {
    class FFTWDoubleEngine : public FFTWEngine<double, fftw_complex, fftw_plan>
    {
    public:
        FFTWDoubleEngine(unsigned nEta, unsigned nPhi,
                         unsigned optimization = FFTW_PATIENT);
        virtual ~FFTWDoubleEngine();
    };
}

#include "fftjet/FFTWDoubleEngine.icc"

#endif // FFTJET_FFTWDOUBLEENGINE_HH_
