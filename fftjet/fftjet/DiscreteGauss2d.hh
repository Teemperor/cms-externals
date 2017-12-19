//=========================================================================
// DiscreteGauss2d.hh
//
// The Fourier transform of the Gaussian kernel, corrected for the grid
// discretization effects
//
// I. Volobouev
// May 2009
//=========================================================================

#ifndef FFTJET_DISCRETEGAUSS2D_HH_
#define FFTJET_DISCRETEGAUSS2D_HH_

#include "fftjet/AbsFrequencyKernel.hh"

namespace fftjet {
    class DiscreteGauss2d : public AbsFrequencyKernel
    {
    public:
        // Parameters "sEta" and "sPhi" are scale factors for
        // eta and phi directions, respectively. They have the same
        // meaning as the corresponding "sx", "sy" parameters for
        // the "Gauss2d" kernel.
        //
        // Parameters "nEta" and "nPhi" represent the number of cells
        // in the eta-phi discretization grid with which this kernel
        // will be used.
        //
        DiscreteGauss2d(double sEta, double sPhi,
                        unsigned nEta, unsigned nPhi);
        inline ~DiscreteGauss2d() {}

        inline double sEta() const {return sEta_;}
        inline double sPhi() const {return sPhi_;}
        inline unsigned nEta() const {return nEta_;}
        inline unsigned nPhi() const {return nPhi_;}

        // The following function checks whether a grid with
        // given numbers of eta and phi bins is compatible
        // with this kernel
        inline bool isCompatible(const unsigned etaBins,
                                 const unsigned phiBins) const
            {return nEta_ == etaBins && nPhi_ == phiBins;}

        // The following inherited function must be overriden
        std::complex<double> operator()(int ix, int iy, double scale) const;

    private:
        DiscreteGauss2d();

        const double sEta_;
        const double sPhi_;
        const unsigned nEta_;
        const unsigned nPhi_;
    };
}

#endif // FFTJET_DISCRETEGAUSS2D_HH_
