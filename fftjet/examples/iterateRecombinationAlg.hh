//=========================================================================
// iterateRecombinationAlg.hh
//
// A convenience function for iterating FFTJet energy recombination
// algorithms.
//
// This function makes a few important assumptions on top of what
// FFTJet normally assumes:
//
// 1) The 4-vector class has px(), py(), eta(), and phi() methods.
//    The "eta()" method of the 4-vector and the eta coordinate used
//    for clustering have the same meaning.
//
// 2) The energy recombination scale is set to 1.0/(jet Pt).
//
// 3) The user knows how to obtain an unbiased estimate of the jet Pt
//    for the jets reconstructed by the algorithm (that is, the user
//    knows how to "correct" jets).
//
// I. Volobouev
// May 2009
//=========================================================================

#ifndef ITERATERECOMBINATIONALG_HH_
#define ITERATERECOMBINATIONALG_HH_

#include "fftjet/AbsRecombinationAlg.hh"
#include "fftjet/AbsJetCorrector.hh"

// The function arguments are as follows:
//
//  recoAlg        -- The energy recombination algorithm to iterate
//
//  peaks, eventData,       -- Inputs to the energy recombination algorithm.
//  bgData, nBgEta, nBgPhi     These arguments have the same names and
//                             meanings as the input arguments of the
//                             AbsRecombinationAlg "run" method.
//                             
//  jetCorr        -- "Jet corrector": this functor will be used to convert
//                     the Pt of the reconstructed jet into the "original"
//                     jet Pt.
//
//  maxIterations  -- Maximum number of iterations to perform. The jet
//                    recombination algorithm will be called just once
//                    in case this argument is 0 or 1.
//
//  etaToPhiBandwidthRatio,   -- These parameters define when we should
//  relativePtBandwidth,         conclude that the iterations have converged.
//  convergedDistanceSquared     The code calculates a "distance" between
//                               two consecutive iterations in the
//                               chi-squared sense. When the squared distance
//                               is less than "convergedDistanceSquared", the
//                               iterations are stopped.
//
//  jets, unclustered, unused -- Outputs of the energy recombination
//                               algorithm. These arguments have the same
//                               names and meanings as the output arguments
//                               of the AbsRecombinationAlg "run" method.
//
//  nIterationsPerformed  -- On exit, will be set to the actual number
//                           of times the algorithm was called and
//                           returned successfully.
//
//  finalSquaredDistance  -- The squared distance between the last and
//                           the next to last iterations. Will be set to
//                           -1.0 in case only one iteration was performed.
//
// The function returns the status word produced by the last invokation
// of the algorithm's "run" method. The output variables (jets, unclustered,
// unused, etc) are filled whenever at least one "run" was successful.
// In fact, the output will contain some meaningful results as long as
// *nIterationsPerformed is not 0 on exit, even if the function itself
// returns non-0 status.
//
template<typename Real, typename VectorLike, typename BgData>
int iterateRecombinationAlg(
    fftjet::AbsRecombinationAlg<Real,VectorLike,BgData>& recoAlg,
    const std::vector<fftjet::Peak>& peaks,
    const fftjet::Grid2d<Real>& eventData,
    const BgData* bgData, unsigned nBgEta, unsigned nBgPhi,
    const fftjet::AbsJetCorrector<fftjet::RecombinedJet<VectorLike> >& jetCorr,
    unsigned maxIterations, double etaToPhiBandwidthRatio,
    double relativePtBandwidth, double convergedDistanceSquared,
    std::vector<fftjet::RecombinedJet<VectorLike> >* jets,
    VectorLike* unclustered, double* unused, unsigned* nIterationsPerformed,
    double* finalSquaredDistance);

#include "iterateRecombinationAlg.icc"

#endif // ITERATERECOMBINATIONALG_HH_
