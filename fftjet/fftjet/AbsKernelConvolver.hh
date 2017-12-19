//=========================================================================
// AbsKernelConvolver.hh
//
// Base class for performing kernel convolutions by FFT. This class manages
// a collection of kernel scans, images and, if requested, the data
// normalization curves. The actual building of kernel images is up to
// derived classes.
//
// Normally, all scans should be performed using the same kernel or
// set of kernels, and with the same grid dimensions. The scans should
// differ only by the scale.
//
// I. Volobouev
// May 2008
//=========================================================================

#ifndef FFTJET_ABSKERNELCONVOLVER_HH_
#define FFTJET_ABSKERNELCONVOLVER_HH_

#include <map>

#include "fftjet/AbsConvolverBase.hh"
#include "fftjet/AbsFFTEngine.hh"
#include "fftjet/KernelData.hh"

namespace fftjet {
    // Scan information for multiple scales.
    // This class will destroy all KernelData
    // data buffers in its own destructor.
    template<typename Real, typename Complex>
    class AbsKernelConvolver : public AbsConvolverBase<Real>
    {
    public:
        // This AbsKernelConvolver object will not own the AbsFFTEngine.
        // It is the responsibility of the user of this class to ensure
        // that this object is not destroyed before the AbsKernelConvolver
        // itself is destroyed.
        AbsKernelConvolver(const AbsFFTEngine<Real,Complex>* fftEngine,
                           unsigned minFixBin, unsigned maxFixBin);
        virtual ~AbsKernelConvolver();

        // Some trivial inspectors
        inline bool fixingEfficiency() const {return fixEfficiency;}
        inline unsigned minFixedBin() const {return minFixBin;}
        inline unsigned maxFixedBin() const {return maxFixBin;}
        inline unsigned nEtaFFT() const {return fftEngine->nEta();}
        inline unsigned nPhiFFT() const {return fftEngine->nPhi();}

        // In the next two methods, the "nEta" and "nPhi" arguments 
        // which define array dimensionalities must correspond to the
        // grid dimensions of the FFT engine used to construct this object.
        // "setEventData" must be called at least once before any call
        // to "convolveWithKernel".
        void setEventData(const Real* data, unsigned nEta, unsigned nPhi);
        void convolveWithKernel(double scale, Real* result, unsigned nEta,
                                unsigned nPhi);

        // There is no accessor method to get the event data back.
        // This is because there is no particular reason for this
        // class to keep the original event data in case the efficiency
        // needs no fixing. In this case only the Fourier image
        // of the data is stored internally.

        // The next method simply requests construction of kernel data
        // (scan, image, and, possibly, normalization) for the given scale.
        // After this, "isScaleProcessed" method will return "true" for
        // the given scale, and scan/image/normalization accessors will
        // return valid pointers. It is not necessary to call this function
        // explicitly if all you want is to convolve the data with
        // the kernel at the given scale.
        void processScale(double scale);

        // Find out how many scales were processed so far
        unsigned nProcessedScales() const;

        // Return the processed scales in the reverse sorted order
        void getProcessedScales(std::vector<double>* scales) const;

        // Check whether a certain scale is already processed
        bool isScaleProcessed(double scale) const;

        // The following methods allow the user to examine the kernel
        // scans/images/normalization curves for the given scale.
        // These methods return NULL if the scale was never processed.
        //
        // Note that the pointer to KernelData<Real,Complex> can be
        // invalidated whenever a new scale is added to the system,
        // but pointers to the scan, image, and normalization curve
        // will remain valid.
        //
        const KernelData<Real,Complex>* getKernelData(double scale) const;
        const Real* getKernelScan(double scale) const;
        const Complex* getKernelImage(double scale) const;
        const Real* getNormalization(double scale) const;

        // The following function returns 0.0 in case the scale
        // was never processed
        double getScanArea(double scale) const;

    protected:
        // The following helper function can be used to build
        // the normalization curve from the scanned kernel
        Real* buildEffNorm(const Real* scan) const;

        // The following function must be implemented by the derived classes.
        // It must build the kernel scan, image and, if requested, the inverse
        // efficiency curve for data normalization. The image and scan buffers
        // must be allocated on the heap using fftEngine->allocateComplex()
        // and "new Real[]", respectively. Construction of the normalization
        // curve can be performed by the "buildEffNorm" function if the kernel
        // image is available in the normal space (not Fourier transformed).
        virtual KernelData<Real,Complex> buildKernelImage(double scale) = 0;

        // Input arguments
        const AbsFFTEngine<Real,Complex>* const fftEngine;
        const unsigned minFixBin;
        const unsigned maxFixBin;
        const bool fixEfficiency;

        // Complex buffer which can be used by the derived classes
        // inside the "buildKernelImage" function. It will be allocated
        // before "buildKernelImage" is called.
        Complex* complexBuffer;

    private:
        AbsKernelConvolver();
        AbsKernelConvolver(const AbsKernelConvolver&);
        AbsKernelConvolver& operator=(const AbsKernelConvolver&);

        // Collection of kernel images and normalization curves
        std::map<double, KernelData<Real,Complex> > kernelProperties;

        // Various memory buffers
        Complex* dataImage;
        const Real* dataBuffer;
        mutable Real* effBuf;
    };
}

#include "fftjet/AbsKernelConvolver.icc"

#endif // FFTJET_ABSKERNELCONVOLVER_HH_
