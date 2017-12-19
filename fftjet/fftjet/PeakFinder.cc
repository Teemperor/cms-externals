#include <cassert>
#include "fftjet/PeakFinder.hh"

namespace fftjet {
    const int PeakFinder::plateau3x3[3][3] = {        
        {2, 2, 2},
        {2, 1, 2},
        {2, 2, 2}
    };

    const int PeakFinder::ridges3x3[4][3][3] = {        
        {
            {0, 0, 0},
            {2, 2, 2},
            {0, 0, 0}
        },{
            {0, 2, 0},
            {0, 2, 0},
            {0, 2, 0}
        },{
            {2, 0, 0},
            {0, 2, 0},
            {0, 0, 2}
        },{
            {0, 0, 2},
            {0, 2, 0},
            {2, 0, 0}
        }
    };

    PeakFinder::PeakFinder(const double peakHeightCut,
                           const bool subCellRes,
                           const unsigned minEtaBin,
                           const unsigned maxEtaBin,
                           const bool printFitWarnings)
        : peakHeightCutoff_(peakHeightCut),
          minEtaBin_(minEtaBin),
          maxEtaBin_(maxEtaBin),
          mask(0),
          nbins(0),
          subCellResolution_(subCellRes),
          printFitWarnings_(printFitWarnings)
    {
    }

    PeakFinder::~PeakFinder()
    {
        delete [] mask;
    }

    PeakFinder::PeakFinder(const PeakFinder& r)
        : peakHeightCutoff_(r.peakHeightCutoff_),
          minEtaBin_(r.minEtaBin_),
          maxEtaBin_(r.maxEtaBin_),
          mask(0),
          nbins(0),
          subCellResolution_(r.subCellResolution_),
          printFitWarnings_(r.printFitWarnings_)
    {
    }

    void PeakFinder::allocateMaskBuffer(const unsigned nEta,
                                        const unsigned nPhi)
    {
        if (nEta*nPhi > nbins)
        {
            delete [] mask;
            mask = 0;
            nbins = 0;
            mask = new int[nEta*nPhi];
            nbins = nEta*nPhi;
        }
    }

    void PeakFinder::clearMaskBuffer(const unsigned nEta,
                                     const unsigned nPhi)
    {
        const unsigned n = nEta*nPhi;
        assert(n <= nbins);
        for (unsigned i=0; i<n; ++i)
            mask[i] = 0;
    }

    bool PeakFinder::match3by3(const unsigned ieta, const unsigned iphi,
                               const unsigned nPhi, const int pattern[3][3])
    {
        const int* mm1(mask + (ieta - 1)*nPhi);
        const int* m(mask + ieta*nPhi);
        const int* mp1(mask + (ieta + 1)*nPhi);
        const unsigned iphim1(iphi ? iphi - 1 : nPhi - 1);
        const unsigned iphip1((iphi + 1) % nPhi);
        return 
            pattern[0][0] == mm1[iphim1] &&
            pattern[0][1] == mm1[iphi] &&
            pattern[0][2] == mm1[iphip1] &&
            pattern[1][0] == m[iphim1] &&
            pattern[1][1] == m[iphi] &&
            pattern[1][2] == m[iphip1] &&
            pattern[2][0] == mp1[iphim1] &&
            pattern[2][1] == mp1[iphi] &&
            pattern[2][2] == mp1[iphip1];
    }

    void PeakFinder::tag3by3(const unsigned ieta, const unsigned iphi,
                             const unsigned nPhi, const int tagValue)
    {
        int* mm1(mask + (ieta - 1)*nPhi);
        int* m(mask + ieta*nPhi);
        int* mp1(mask + (ieta + 1)*nPhi);
        const unsigned iphim1(iphi ? iphi - 1 : nPhi - 1);
        const unsigned iphip1((iphi + 1) % nPhi);

        mm1[iphim1] = 0;
        mm1[iphi] = 0;
        mm1[iphip1] = 0;
        m[iphim1] = 0;
        m[iphi] = tagValue;
        m[iphip1] = 0;
        mp1[iphim1] = 0;
        mp1[iphi] = 0;
        mp1[iphip1] = 0;
    }

    void PeakFinder::processMask(const unsigned nEta, const unsigned nPhi,
                                 const unsigned minEtaBin,
                                 const unsigned maxEtaBin)
    {
        assert(nEta*nPhi <= nbins);

        // Recognize plateaus/ridges in a 3x3 window.
        // We assume that these patterns are infrequent.
        for (unsigned ieta = minEtaBin; ieta < maxEtaBin; ++ieta)
        {
            const int* m = mask + ieta*nPhi;
            for (unsigned iphi = 0; iphi < nPhi; ++iphi)
            {
                switch (m[iphi])
                {
                case 1:
                    // Recognize the 3x3 plateau
                    if (m[(iphi + 1) % nPhi] == 2)
                        if (match3by3(ieta, iphi, nPhi, plateau3x3))
                            tag3by3(ieta, iphi, nPhi, 4);
                    break;
                case 2:
                    // Recognize various 3x3 ridges
                    for (unsigned ir=0;
                         ir<sizeof(ridges3x3)/sizeof(ridges3x3[0]); ++ir)
                        if (match3by3(ieta, iphi, nPhi, ridges3x3[ir]))
                        {
                            tag3by3(ieta, iphi, nPhi, 3);
                            break;
                        }
                    break;
                default:
                    break;
                }
            }
        }

        // Do something with ridges which remain in a 2x2 window
        for (unsigned ieta = minEtaBin; ieta < maxEtaBin; ++ieta)
        {
            const int* m = mask + ieta*nPhi;
            for (unsigned iphi = 0; iphi < nPhi; ++iphi)
                if (m[iphi] == 2)
                    tag3by3(ieta, iphi, nPhi, 3);
        }
    }
}
