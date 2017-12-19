//=========================================================================
// AbsSequentialConvolver.hh
//
// Base class for performing kernel convolutions sequentially (first,
// convolve each eta bin, then each phi bin). This allows for different
// kernel scales in different eta bins.
//
// This class manages a collection of kernel scans, images and, if
// requested, the data normalization curves. The actual building of
// kernel images is up to derived classes.
//
// Normally, all scans should be performed using the same kernel or
// set of kernels, and with the same grid dimensions. The scans should
// differ only by scale.
//
// I. Volobouev
// November 2009
//=========================================================================

#ifndef FFTJET_ABSSEQUENTIALCONVOLVER_HH_
#define FFTJET_ABSSEQUENTIALCONVOLVER_HH_

#include <map>

#include "fftjet/AbsConvolverBase.hh"
#include "fftjet/AbsFFTEngine.hh"
#include "fftjet/KernelData.hh"

namespace fftjet {
    // There is one KernelData<Real,Complex> element for each eta bin
    // and one eta kernel for all phi bins
    template<typename Real, typename Complex>
    class SequentialKernelData : public std::vector<KernelData<Real,Complex> >
    {
    public:
        inline SequentialKernelData(const KernelData<Real,Complex>& etaKernel)
            : etaKernel_(etaKernel) {}
        inline const KernelData<Real,Complex>& etaKernelData() const
            {return etaKernel_;}

    private:
        SequentialKernelData();
        KernelData<Real,Complex> etaKernel_;
    };

    // Scan information for multiple scales.
    // This class will destroy all KernelData
    // data buffers in its own destructor.
    template<typename Real, typename Complex>
    class AbsSequentialConvolver : public AbsConvolverBase<Real>
    {
    public:
        // This AbsSequentialConvolver object will not own the AbsFFTEngine
        // objects. It is the responsibility of the user of this class
        // to ensure that these objects are not destroyed before the
        // AbsSequentialConvolver itself is destroyed.
        //
        // Note that in both engines the nEta() function should return 1
        // and nPhi() should return the number of bins. The number of
        // elements in the "phiScales" vector should be equal to nPhi()
        // of the "etaEngine". The scale for each eta bin will be obtained
        // by multiplying the overall scale by the corresponding bin scale.
        //
        AbsSequentialConvolver(const AbsFFTEngine<Real,Complex>* etaEngine,
                               const AbsFFTEngine<Real,Complex>* phiEngine,
                               const std::vector<double>& phiScales,
                               unsigned minEtaBin, unsigned maxEtaBin,
                               bool fixEfficiency);
        virtual ~AbsSequentialConvolver();

        // Some trivial inspectors
        inline bool fixingEfficiency() const {return fixEfficiency;}
        inline unsigned minFixedBin() const {return minEtaBin;}
        inline unsigned maxFixedBin() const {return maxEtaBin;}
        inline unsigned nEtaFFT() const {return etaEngine->nPhi();}
        inline unsigned nPhiFFT() const {return phiEngine->nPhi();}

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
        // class to keep the original event data. Only the Fourier
        // image of the data is stored internally.

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
        // Note that the pointer to SequentialKernelData can be
        // invalidated whenever a new scale is added to the system,
        // but pointers to the scan and image curves will remain valid.
        //
        const SequentialKernelData<Real,Complex>* getKernelData(
            double scale) const;
        const Real* getKernelScanPhi(double scale, unsigned etaBin) const;
        const Complex* getKernelImagePhi(double scale, unsigned etaBin) const;
        double getScanAreaPhi(double scale, unsigned etaBin) const;

        const Real* getKernelScanEta(double scale) const;
        const Complex* getKernelImageEta(double scale) const;
        const Real* getNormalization(double scale) const;
        double getScanAreaEta(double scale) const;

    protected:
        // The following helper function can be used to build
        // the normalization curve from the scanned kernel
        Real* buildEffNorm(const Real* scan, unsigned scanLength) const;

        // The following functions must be implemented by the derived classes.
        // It must build the kernel scan, image and, if requested, the inverse
        // efficiency curve for data normalization (the latter must be done
        // in case minFixBin < maxFixBin).
        virtual KernelData<Real,Complex> buildKernelImageEta(
            double scale, const AbsFFTEngine<Real,Complex>* engine,
            unsigned minFixBin, unsigned maxFixBin) = 0;
        virtual KernelData<Real,Complex> buildKernelImagePhi(
            double scale, const AbsFFTEngine<Real,Complex>* engine) = 0;

        // Complex buffer which can be used by the derived classes.
        // The buffer will be sufficiently large to accomodate the image
        // of either phi or eta kernel.
        Complex* complexBuffer;

    private:
        AbsSequentialConvolver();
        AbsSequentialConvolver(const AbsSequentialConvolver&);
        AbsSequentialConvolver& operator=(const AbsSequentialConvolver&);

        // Collection of kernel images and normalization curves
        std::map<double, SequentialKernelData<Real,Complex> > kernelProperties;
        typedef typename std::map<double, SequentialKernelData<Real,Complex> 
            >::iterator MapIterator;
        typedef typename std::map<double, SequentialKernelData<Real,Complex> 
            >::const_iterator ConstMapIterator;

        // Unique kernel images for the phi convolutions
        std::map<double, KernelData<Real,Complex> > uniquePhiKernels;

        // Collection of buffers for the phi images in each eta bin
        std::vector<Complex*> dataImages;

        // Work buffer for eta convolutions
        Real* workbuf;

        // Input arguments
        const AbsFFTEngine<Real,Complex>* const etaEngine;
        const AbsFFTEngine<Real,Complex>* const phiEngine;
        const std::vector<double> phiScales;
        const unsigned minEtaBin;
        const unsigned maxEtaBin;
        const bool fixEfficiency;

        // See if we have data
        bool haveData;
    };
}

#include "fftjet/AbsSequentialConvolver.icc"

#endif // FFTJET_ABSSEQUENTIALCONVOLVER_HH_
