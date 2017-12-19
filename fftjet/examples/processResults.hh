// This code pretty-prints the jets recostructed by FFTJet
//
// I. Volobouev
// May 2009

#ifndef PROCESSRESULTS_HH_
#define PROCESSRESULTS_HH_

#include <vector>

#include "fftjet/RecombinedJet.hh"
#include "fftjet_typedefs.hh"

struct SimpleEvent;

void processResults(unsigned eventNumber, double etaToPhiBandwidthRatio,
               const SimpleEvent& mcEvent,
               const std::vector<fftjet::RecombinedJet<VectorLike> >& recoJets,
               const VectorLike& unclustered, double unused);

#endif // PROCESSRESULTS_HH_
