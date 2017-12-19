#include <cassert>
#include <cfloat>
#include <cmath>

#include "SimpleEvent.hh"
#include "PythiaJetGun.hh"

using namespace geom3;
using namespace PythiaJetGun;

rk::P4 SimulatedJet::jetSum() const
{
    rk::P4 sum;
    const unsigned npart = jetParticles.size();
    for (unsigned j=0; j<npart; ++j)
        sum += jetParticles[j].p4();
    return sum;
}

rk::P4 SimulatedJet::chargedSum() const
{
    rk::P4 sum;
    const unsigned npart = jetParticles.size();
    for (unsigned j=0; j<npart; ++j)
        if (jetParticles[j].charge() != 0.0)
            sum += jetParticles[j].p4();
    return sum;
}

rk::P4 SimpleEvent::closestJet(const double eta, const double phi) const
{
    const unsigned njets = size();
    assert(njets);

    double bestDr = DBL_MAX;
    rk::P4 bestJet;

    const geom3::UnitVector3 direction(cos(phi), sin(phi), sinh(eta));
    for (unsigned i=0; i<njets; ++i)
    {
        const SimulatedJet& jet((*this)[i]);
        rk::P4 jetSum(jet.jetSum());
        const double dr = geom3::deltaR(direction, jetSum.momentum());
        if (dr < bestDr)
        {
            bestDr = dr;
            bestJet = jetSum;
        }
    }

    return bestJet;
}

unsigned SimpleEvent::closestJetNumber(const unsigned thisJetNumber) const
{
    rk::P4 thisSum(at(thisJetNumber).jetSum());
    const double eta = thisSum.eta();
    const double phi = thisSum.phi();
    const unsigned njets = size();

    double bestDr = DBL_MAX;
    unsigned bestJet = njets;

    const geom3::UnitVector3 direction(cos(phi), sin(phi), sinh(eta));
    for (unsigned i=0; i<njets; ++i)
        if (i != thisJetNumber)
        {
            const SimulatedJet& jet((*this)[i]);
            rk::P4 jetSum(jet.jetSum());
            const double dr = geom3::deltaR(direction, jetSum.momentum());
            if (dr < bestDr)
            {
                bestDr = dr;
                bestJet = i;
            }
        }

    return bestJet;
}

double SimpleEvent::minDr() const
{
    const unsigned njets = size();
    if (njets < 2)
        return 0.0;

    double dr, drmin = DBL_MAX;
    for (unsigned i=0; i<njets-1; ++i)
    {
        const geom3::Vector3 p((*this)[i].jetSum().momentum());
        for (unsigned j=i+1; j<njets; ++j)
        {
            dr = geom3::deltaR(p, (*this)[j].jetSum().momentum());
            if (dr < drmin)
                drmin = dr;
        }
    }
    return drmin;
}

bool phiAtRadius(const Particle& particle, const double B,
                 const double r, double *phi)
{
    assert(r >= 0.0);
    *phi = particle.p4().phi();
    if (B == 0.0 || r == 0.0)
        return true;
    if (particle.charge() == 0.0)
        return true;
    const double a = particle.gyrationRadius(B);
    const double twoa = 2.0*fabs(a);
    if (r > twoa)
    {
        *phi = DBL_MAX;
        return false;
    }
    if (r == twoa)
        *phi += a/twoa*M_PI;
    else
        *phi += asin(r/2.0/a);
    return true;
}
