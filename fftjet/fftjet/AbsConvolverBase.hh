//=========================================================================
// AbsConvolverBase.hh
//
// Abstract base class for performing kernel convolutions
//
// I. Volobouev
// November 2009
//=========================================================================

#ifndef FFTJET_ABSCONVOLVERBASE_HH_
#define FFTJET_ABSCONVOLVERBASE_HH_

#include <vector>

namespace fftjet {
    template<typename Real>
    class AbsConvolverBase
    {
    public:
        virtual ~AbsConvolverBase() {}

        // Some inspectors of the functionality
        virtual bool fixingEfficiency() const = 0;
        virtual unsigned minFixedBin() const = 0;
        virtual unsigned maxFixedBin() const = 0;
        virtual unsigned nEtaFFT() const = 0;
        virtual unsigned nPhiFFT() const = 0;

        // In the next two methods, the "nEta" and "nPhi" arguments 
        // which define array dimensionalities must correspond to the
        // grid dimensions of the FFT engine used to construct this object.
        // "setEventData" must be called at least once before any call
        // to "convolveWithKernel".
        virtual void setEventData(const Real* data,
                                  unsigned nEta, unsigned nPhi) = 0;
        virtual void convolveWithKernel(double scale, Real* result,
                                        unsigned nEta, unsigned nPhi) = 0;

        // Functions related to scale management
        virtual void processScale(double scale) = 0;

        // Find out how many scales were processed so far
        virtual unsigned nProcessedScales() const = 0;

        // Return the processed scales in the reverse sorted order
        virtual void getProcessedScales(std::vector<double>* scales) const = 0;

        // Check whether a certain scale is already processed
        virtual bool isScaleProcessed(double scale) const = 0;
    };
}

#endif // FFTJET_ABSCONVOLVERBASE_HH_
