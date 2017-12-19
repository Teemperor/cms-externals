//=========================================================================
// DiscreteGauss1d.hh
//
// The Fourier transform of the Gaussian kernel, corrected for the grid
// discretization effects
//
// I. Volobouev
// November 2009
//=========================================================================

#ifndef FFTJET_DISCRETEGAUSS1D_HH_
#define FFTJET_DISCRETEGAUSS1D_HH_

#include "fftjet/AbsFrequencyKernel1d.hh"

namespace fftjet {
    class DiscreteGauss1d : public AbsFrequencyKernel1d
    {
    public:
        // Parameter "sPhi" is the scale factors which has
        // the same meaning as the corresponding "sx" parameter
        // for the "Gauss1d" kernel.
        //
        // Parameter "nPhi" represents the number of cells
        // in the discretization grid with which this kernel
        // will be used.
        //
        DiscreteGauss1d(double sPhi, unsigned nPhi);
        inline ~DiscreteGauss1d() {}

        inline double sPhi() const {return sPhi_;}
        inline unsigned nPhi() const {return nPhi_;}

        // The following inherited function must be overriden
        std::complex<double> operator()(int ix, double scale) const;

    private:
        DiscreteGauss1d();

        const double sPhi_;
        const unsigned nPhi_;
    };
}

#endif // FFTJET_DISCRETEGAUSS1D_HH_
