//========================================================================
// convertPeakCoords.hh
//
// The following function will convert the peak positions from the units
// of cell size into absolute eta/phi units. For input grids with
// reasonable phiBin0Edge values the phi angle will end up in the
// interval [0, 2*Pi).
//
// I. Volobouev
// May 2008
//========================================================================

#ifndef FFTJET_CONVERTPEAKCOORDS_HH_
#define FFTJET_CONVERTPEAKCOORDS_HH_

#include <vector>
#include <cmath>

#include "fftjet/Grid2d.hh"
#include "fftjet/Peak.hh"

namespace fftjet {
    template<typename Real>
    void convertPeakCoords(const Grid2d<Real>& grid,
                           std::vector<Peak>* peaks)
    {
        const double twopi = 2.0*M_PI;
        const unsigned npeaks = peaks->size();
        const double etaWidth = grid.etaBinWidth();
        const double phiWidth = grid.phiBinWidth();
        const double h0 = etaWidth*etaWidth;
        const double h1 = etaWidth*phiWidth;
        const double h2 = phiWidth*phiWidth;
        double hessian[3];

        for (unsigned ijet=0; ijet<npeaks; ++ijet)
        {
            Peak& peak((*peaks)[ijet]);

            // Convert peak position
            double phi(grid.phiFromBinNum(peak.phi()));
            if (phi < 0.0)
                phi += twopi;
            else if (phi >= twopi)
                phi -= twopi;
            peak.setEtaPhi(grid.etaFromBinNum(peak.eta()), phi);

            // Convert peak Hessian
            peak.hessian(hessian);
            hessian[0] /= h0;
            hessian[1] /= h1;
            hessian[2] /= h2;
            peak.setHessian(hessian);
        }
    }
}

#endif // FFTJET_CONVERTPEAKCOORDS_HH_
