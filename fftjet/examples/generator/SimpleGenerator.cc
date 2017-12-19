#include <algorithm>
#include <cassert>
#include <cmath>

#include "SimpleGenerator.hh"
#include "SimpleEvent.hh"

using namespace PythiaJetGun;

SimpleGenerator::SimpleGenerator(unsigned njets, unsigned partonCode,
                                 double ptMin, double ptMax,
                                 double etaMin, double etaMax,
                                 double phiMin, double phiMax,
                                 bool useOriginalSingleParticleGun,
                                 bool generateFlatInLogPt)
    : njets_(njets),
      partonCode_(partonCode),
      ptMin_(ptMin),
      ptMax_(ptMax),
      etaMin_(etaMin),
      etaMax_(etaMax),
      phiMin_(phiMin),
      phiMax_(phiMax),
      single_(useOriginalSingleParticleGun),
      flatlog_(generateFlatInLogPt)
{
    assert(ptMin_ > 0.0);
    assert(ptMax_ > 0.0);
    logMin_ = log(ptMin_);
    logMax_ = log(ptMax_);
}

void SimpleGenerator::generate(SimpleEvent* ev) const
{
    ev->clear();
    ev->reserve(njets_);

    SimulatedJet simj;

    for (unsigned i=0; i<njets_; ++i)
    {
        double pt;
        if (flatlog_)
            pt = exp(logMin_ + (logMax_ - logMin_)*pyrandom());
        else
            pt = ptMin_ + (ptMax_ - ptMin_)*pyrandom();
        const double eta = etaMin_ + (etaMax_ - etaMin_)*pyrandom();
        const double phi = phiMin_ + (phiMax_ - phiMin_)*pyrandom();

        if (single_)
            shoot1(partonCode_, pt, eta, phi,
                   &simj.parton, &simj.jetParticles);
        else
            shoot2(partonCode_, pt, eta, phi,
                   &simj.parton, &simj.jetParticles);
        ev->push_back(simj);
    }

    std::sort(ev->begin(), ev->end(), std::greater<SimulatedJet>());
}
