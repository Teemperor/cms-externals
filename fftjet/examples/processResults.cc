#include <cstdio>

#include "processResults.hh"
#include "SimpleEvent.hh"

#include "matchOneToOne.hh"

namespace {
    // The following functor calculates the distance in the eta-phi
    // space between a Monte Carlo jet and a reconstructed jet
    class EventToJetDistance
    {
    public:
        inline explicit EventToJetDistance(const double etaToPhiBandwidthRatio)
            : bwEta(sqrt(etaToPhiBandwidthRatio)), bwPhi(1.0/bwEta)
        {
            assert(bwEta > 0.0);
        }

        inline double operator()(const SimulatedJet& ev,
                      const fftjet::RecombinedJet<VectorLike>& jet) const
        {
            geom3::Vector3 p(ev.jetSum().momentum());
            const double deta = (jet.vec().eta() - p.eta())/bwEta;
            const double dphi = geom3::deltaPhi(jet.vec().phi(),p.phi())/bwPhi;
            return sqrt(deta*deta + dphi*dphi);
        }

    private:
        EventToJetDistance();
        const double bwEta;
        const double bwPhi;
    };
}


// Print the results. We will show the original Pythia jets
// together with the corresponding reconstructed jets.
void processResults(
    const unsigned eventNumber, const double etaToPhiBandwidthRatio,
    const SimpleEvent& mcEvent,
    const std::vector<fftjet::RecombinedJet<VectorLike> >& recoJets,
    const VectorLike& /* unclustered */, const double unused)
{
    static bool printHeader = true;
    if (printHeader)
    {
        printf("========================================================\n");
        printf("    Original Pythia Jet     |  FFTJet Reconstructed Jet \n");
        printf("   Pt        Eta      Phi   |    Pt        Eta      Phi \n");
        printf("========================================================\n");
        printHeader = false;
    }

    // Print some info about the event
    printf("Event %u, %u jets reconstructed, %7.3f GeV unused/noise\n",
           eventNumber, static_cast<unsigned>(recoJets.size()), unused);

    // Match the Monte Carlo jets to the reconstructed jets
    std::vector<int> matchTable;
    matchOneToOne(mcEvent, recoJets,
                  EventToJetDistance(etaToPhiBandwidthRatio), &matchTable);

    // Go over the MC jets and print them together with the matched
    // reconstructed jets
    const int nMCjets = mcEvent.size();
    for (int imc=0; imc<nMCjets; ++imc)
    {
        const rk::P4 mcJet(mcEvent[imc].jetSum());
        printf("%8.3f %8.3f %8.3f  | ", mcJet.pt(), mcJet.eta(), mcJet.phi());

        // Do we have a matched jet?
        const int i = matchTable[imc];
        if (i >= 0)
        {
            const rk::P4 recoJet(recoJets[i].vec());
            printf("%8.3f %8.3f %8.3f\n",
                   recoJet.pt(), recoJet.eta(), recoJet.phi());
        }
        else
            printf("%8.3f %8.3f %8.3f\n", 0.0, 0.0, 0.0);
    }

    // Print reconstructed jets which are not matched to MC jets (if any)
    const int nRecoJets = recoJets.size();
    if (nRecoJets > nMCjets)
    {
        for (int ireco=0; ireco<nRecoJets; ++ireco)
            if (std::find(matchTable.begin(), matchTable.end(), ireco) == 
                matchTable.end())
            {
                const rk::P4 recoJet(recoJets[ireco].vec());
                printf("%8.3f %8.3f %8.3f  | %8.3f %8.3f %8.3f\n", 0.0, 0.0,
                       0.0, recoJet.pt(), recoJet.eta(), recoJet.phi());
            }
    }
}
