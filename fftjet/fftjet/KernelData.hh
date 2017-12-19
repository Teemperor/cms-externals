//=========================================================================
// KernelData.hh
//
// Scan information for a single scale. This class will not own
// the pointers. "normalization()" method returns the normalization
// curve (or NULL pointer in case there is no normalization data).
//
// I. Volobouev
// November 2009
//=========================================================================

#ifndef FFTJET_KERNELDATA_HH_
#define FFTJET_KERNELDATA_HH_

namespace fftjet {
    template<typename Real, typename Complex>
    class KernelData
    {
    public:
        // The "fftScanArea" parameter is not the area under the scan
        // curve (which is usually normalized). Instead, it should
        // represent the area returned by the "scanFFT" function
        // of the kernel. If it is very different from 1.0, it means
        // there are significant distortions of the kernel shape
        // due to binning.
        KernelData(const Real* scan, const Complex* image,
                   const Real* norm, double fftScanArea);

        inline const Real* scan() const {return scan_;}
        inline const Complex* image() const {return image_;}
        inline const Real* normalization() const {return normalization_;}
        inline double fftScanArea() const {return scanArea_;}

    private:
        KernelData();

        const Real* scan_;
        const Complex* image_;
        const Real* normalization_;
        double scanArea_;
    };
}

#include "fftjet/KernelData.icc"

#endif // FFTJET_KERNELDATA_HH_
