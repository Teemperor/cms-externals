#include <iostream>
#include <cassert>
#include <cstdlib>

#include "fftjet/CUFFTFloatEngine.hh"
#include "cuda_runtime_api.h"

namespace fftjet {
    void CUFFTFloatEngine::checkFFTStatus(const cufftResult status) const
    {
        switch (status)
        {
        case CUFFT_SUCCESS:
            break;

        case CUFFT_SETUP_FAILED:
            std::cerr << "CUFFT library failed to initialize. Aborting."
                      << std::endl;
            std::cerr.flush();
            abort();
            break;

        case CUFFT_INVALID_SIZE:
            std::cerr << "Dimensions " << nEta_ << " x " << nPhi_
                      << " are not supported by CUFFT. Aborting."
                      << std::endl;
            std::cerr.flush();
            abort();
            break;

        case CUFFT_INVALID_TYPE:
            // This should never happen
            assert(!"Invalid CUFFT transform type");
            break;

        case CUFFT_ALLOC_FAILED:
            std::cerr << "CUFFT can't allocate GPU resources. Aborting."
                      << std::endl;
            std::cerr.flush();
            abort();
            break;

        case CUFFT_EXEC_FAILED:
            std::cerr << "CUFFT failed to execute transform on GPU. Aborting."
                      << std::endl;
            std::cerr.flush();
            abort();
            break;

        case CUFFT_INVALID_PLAN:
            // This should never happen
            assert(!"Invalid CUFFT plan");
            break;

        case CUFFT_INVALID_VALUE:
            // This should never happen
            assert(!"Invalid CUFFT argument");
            break;

        default:
            // This should never happen
            assert(!"Unexpected CUFFT status");
            break;
        }
    }

    CUFFTFloatEngine::CUFFTFloatEngine(const unsigned nEta,
                                       const unsigned nPhi)
        : AbsFFTEngine<float, cufftComplex>(nEta, nPhi),
          in(0),
          out(0)
    {
        cudaError_t status = cudaMalloc((void**)&in, nEta*nPhi*sizeof(float));
        assert(status == 0);
        assert(in);

        status = cudaMalloc((void**)&out,nEta*(nPhi/2+1)*sizeof(cufftComplex));
        assert(status == 0);
        assert(out);

        if (nEta == 1)
        {
            checkFFTStatus(cufftPlan1d(&pf, nPhi, CUFFT_R2C, 1));
            checkFFTStatus(cufftPlan1d(&pb, nPhi, CUFFT_C2R, 1));
        }
        else
        {
            checkFFTStatus(cufftPlan2d(&pf, nEta, nPhi, CUFFT_R2C));
            checkFFTStatus(cufftPlan2d(&pb, nEta, nPhi, CUFFT_C2R));
        }
    }

    CUFFTFloatEngine::~CUFFTFloatEngine()
    {
        cufftDestroy(pb);
        cufftDestroy(pf);
        cudaFree(out);
        cudaFree(in);
    }

    void CUFFTFloatEngine::transformForward(
        const float* from, cufftComplex* to) const
    {
        assert(from);
        assert(to);
        cudaError_t status = cudaMemcpy(
            in, from, nEta_ * nPhi_ * sizeof(float),
            cudaMemcpyHostToDevice);
        assert(status == 0);
        checkFFTStatus(cufftExecR2C(pf, in, out));
        status = cudaMemcpy(to, out, nEta_*(nPhi_/2+1)*sizeof(cufftComplex),
                            cudaMemcpyDeviceToHost);
        assert(status == 0);
    }

    void CUFFTFloatEngine::transformBack(
        const cufftComplex* from, float* to) const
    {
        assert(from);
        assert(to);
        cudaError_t status = cudaMemcpy(
            out, from, nEta_*(nPhi_/2+1)*sizeof(cufftComplex),
            cudaMemcpyHostToDevice);
        assert(status == 0);
        checkFFTStatus(cufftExecC2R(pb, out, in));
        status = cudaMemcpy(to, in, nEta_*nPhi_*sizeof(float),
                            cudaMemcpyDeviceToHost);
        assert(status == 0);
        const unsigned n2 = nEta_ * nPhi_;
        const float fn2 = n2;
        for (unsigned i=0; i<n2; ++i)
            to[i] /= fn2;
    }

    void CUFFTFloatEngine::multiplyTransforms(
        const cufftComplex* l, const cufftComplex* r, cufftComplex* res) const
    {
        assert(l);
        assert(r);
        assert(res);
        const unsigned nc = nEta_ * (nPhi_/2 + 1);
        for (unsigned i=0; i<nc; ++i)
            res[i] = cuCmulf(l[i], r[i]);
    }

    void CUFFTFloatEngine::addTransforms(const cufftComplex* l,
                                         const cufftComplex* r,
                                         cufftComplex* res) const
    {
        assert(l);
        assert(r);
        assert(res);
        const unsigned nc = nEta_ * (nPhi_/2 + 1);
        for (unsigned i=0; i<nc; ++i)
            res[i] = cuCaddf(l[i], r[i]);
    }

    void CUFFTFloatEngine::divideTransforms(const cufftComplex* l,
                                            const cufftComplex* r,
                                            cufftComplex* res) const
    {
        assert(l);
        assert(r);
        assert(res);
        const unsigned nc = nEta_ * (nPhi_/2 + 1);
        for (unsigned i=0; i<nc; ++i)
            res[i] = cuCdivf(l[i], r[i]);
    }

    void CUFFTFloatEngine::deconvolutionRatio(const cufftComplex* l,
                                              const cufftComplex* r,
                                              cufftComplex* res) const
    {
        assert(l);
        assert(r);
        assert(res);
        const unsigned nc = nEta_ * (nPhi_/2 + 1);
        const cufftComplex zero = make_cuFloatComplex(0.f, 0.f);
        for (unsigned i=0; i<nc; ++i)
        {
            if (cuCrealf(r[i]) == 0.f && cuCimagf(r[i]) == 0.f)
            {
                if (cuCrealf(l[i]) == 0.f && cuCimagf(l[i]) == 0.f)
                    res[i] = zero;
                else
                    assert(!"Complex division by zero");
            }
            else
                res[i] = cuCdivf(l[i], r[i]);
        }
    }

    void CUFFTFloatEngine::subtractTransforms(const cufftComplex* l,
                                              const cufftComplex* r,
                                              cufftComplex* res) const
    {
        assert(l);
        assert(r);
        assert(res);
        const unsigned nc = nEta_ * (nPhi_/2 + 1);
        for (unsigned i=0; i<nc; ++i)
            res[i] = cuCsubf(l[i], r[i]);
    }

    void CUFFTFloatEngine::scaleTransform(const cufftComplex* a,
                                          const float scale,
                                          cufftComplex* res) const
    {
        assert(a);
        assert(res);
        const unsigned nc = nEta_ * (nPhi_/2 + 1);
        const cufftComplex cScale = make_cuFloatComplex(scale, 0.f);
        for (unsigned i=0; i<nc; ++i)
            res[i] = cuCmulf(a[i], cScale);
    }

    void CUFFTFloatEngine::amplitudeSquared(const cufftComplex* a,
                                            cufftComplex* res) const
    {
        assert(a);
        assert(res);
        const unsigned nc = nEta_ * (nPhi_/2 + 1);
        for (unsigned i=0; i<nc; ++i)
            res[i] = cuCmulf(a[i], cuConjf(a[i]));
    }

    double CUFFTFloatEngine::totalPower(const cufftComplex* a) const
    {
        assert(a);
        const unsigned physNPhi(nPhi_/2 + 1);
        double sum = 0.0;
        for (unsigned iEta=0; iEta<nEta_; ++iEta)
            for (unsigned iph=0; iph<nPhi_; ++iph)
            {
                const unsigned iPhi = iph >= physNPhi ? nPhi_ - iph : iph;
                const unsigned idx(iEta*physNPhi + iPhi);
                const float re = cuCrealf(a[idx]);
                const float im = cuCimagf(a[idx]);
                sum += (re*re + im*im);
            }
        return sum;
    }

    void CUFFTFloatEngine::setToReal(cufftComplex* res, float v) const
    {
        assert(res);
        const unsigned nc = nEta_ * (nPhi_/2 + 1);
        for (unsigned i=0; i<nc; ++i)
            res[i] = make_cuFloatComplex(v, 0.f);
    }

    void CUFFTFloatEngine::setTransformPoint(
        cufftComplex* points,
        unsigned iEta, unsigned iPhi,
        const std::complex<double>& value) const
    {
        const unsigned physNPhi(nPhi_/2 + 1);
        iPhi %= nPhi_;
        if (iPhi >= physNPhi)
            iPhi = nPhi_ - iPhi;
        points[(iEta % nEta_)*physNPhi + iPhi] = 
            make_cuFloatComplex(static_cast<float>(value.real()),
                                static_cast<float>(value.imag()));
    }

    std::complex<double> CUFFTFloatEngine::getTransformPoint(
        const cufftComplex* points,
        unsigned iEta, unsigned iPhi) const
    {
        const unsigned physNPhi(nPhi_/2 + 1);
        iPhi %= nPhi_;
        if (iPhi >= physNPhi)
            iPhi = nPhi_ - iPhi;
        const unsigned idx((iEta % nEta_)*physNPhi + iPhi);
        return std::complex<double>(cuCrealf(points[idx]),
                                    cuCimagf(points[idx]));
    }

    cufftComplex* CUFFTFloatEngine::allocateComplex() const
    {
        return reinterpret_cast<cufftComplex *>(
            malloc(nEta_*(nPhi_/2+1)*sizeof(cufftComplex)));
    }

    void CUFFTFloatEngine::destroyComplex(cufftComplex* p) const
    {
        if (p) free(p);
    }
}
