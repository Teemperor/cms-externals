//========================================================================
// scanFrequencyKernel.hh
//
// A standard way to scan frequency kernels
//
// I. Volobouev
// May 2009
//========================================================================

#ifndef FFTJET_SCANFREQUENCYKERNEL_HH_
#define FFTJET_SCANFREQUENCYKERNEL_HH_

#include <cassert>

#include "fftjet/AbsFFTEngine.hh"
#include "fftjet/DiscreteGauss1d.hh"
#include "fftjet/DiscreteGauss2d.hh"

namespace fftjet {
    template<typename Real, typename Complex>
    void scanFrequencyKernel(
        const AbsFFTEngine<Real,Complex>* const engine,
        const AbsFrequencyKernel* kernel,
        const double scale, Complex* to)
    {
        assert(engine);
        assert(kernel);
        assert(to);

        // The "DiscreteGauss2d" kernel is too important
        // to allow for random screw-ups. Force binning
        // compatibility here.
        {
            const DiscreteGauss2d* dg = 
                dynamic_cast<const DiscreteGauss2d*>(kernel);
            if (dg)
                assert(dg->isCompatible(engine->nEta(), engine->nPhi()));
        }

        const int nEta = static_cast<int>(engine->nEta());
        const int nPhi = static_cast<int>(engine->nPhi());
        for (int ieta=0; ieta<nEta; ++ieta)
        {
            const int ie = nEta - ieta < ieta ? ieta - nEta : ieta;
            for (int iphi=0; iphi<nPhi; ++iphi)
            {
                const int ip = nPhi - iphi < iphi ? iphi - nPhi : iphi;
                engine->setTransformPoint(to, ieta, iphi,
                                          (*kernel)(ie, ip, scale));
            }
        }
    }

    template<typename Real, typename Complex>
    void scanFrequencyKernel1d(
        const AbsFFTEngine<Real,Complex>* const engine,
        const AbsFrequencyKernel1d* kernel,
        const double scale, Complex* to)
    {
        assert(engine);
        assert(kernel);
        assert(to);

        // The "DiscreteGauss1d" kernel is too important
        // to allow for random screw-ups. Force binning
        // compatibility here.
        {
            const DiscreteGauss1d* dg = 
                dynamic_cast<const DiscreteGauss1d*>(kernel);
            if (dg)
                assert(dg->nPhi() == engine->nPhi());
        }

        const int nPhi = static_cast<int>(engine->nPhi());
        for (int iphi=0; iphi<nPhi; ++iphi)
        {
            const int ip = nPhi - iphi < iphi ? iphi - nPhi : iphi;
            engine->setTransformPoint(to, 0, iphi, (*kernel)(ip, scale));
        }
    }

    // The following function can be used to normalize
    // the Fourier image, so that its power corresponds
    // to the power of the given spatial image treated
    // as a density
    template<typename Real, typename Complex>
    Real normalizeFrequencyAsDensity(
        const AbsFFTEngine<Real,Complex>* const engine,
        Complex* imageBuf, const Real* timeRep)
    {
        assert(engine);
        assert(imageBuf);
        assert(timeRep);

        const unsigned nEta = engine->nEta();
        const unsigned nPhi = engine->nPhi();
        const unsigned nbins = nEta*nPhi;

        long double sum = 0.0L, sumsq = 0.0L;
        for (unsigned i=0; i<nbins; ++i)
        {
            sum += timeRep[i];
            sumsq += timeRep[i]*timeRep[i];
        }
        assert(sum != 0.0L);

        const double timePower = sumsq/sum/sum;
        const double imagePower = engine->totalPower(imageBuf);
        assert(imagePower > 0.0);

        // The image power should be equal to timePower*nEta*nPhi
        engine->scaleTransform(imageBuf, sqrt(timePower*nbins/imagePower),
                               imageBuf);
        return static_cast<Real>(sum);
    }
}

#endif // FFTJET_SCANFREQUENCYKERNEL_HH_
