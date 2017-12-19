#ifndef SIMPLEEVENT_HH_
#define SIMPLEEVENT_HH_

#include "PythiaJetGun.hh"

// The following structure represents a jet
struct SimulatedJet
{
    // Parton originating the jet
    PythiaJetGun::Particle parton;

    // The constituent jet particles
    std::vector<PythiaJetGun::Particle> jetParticles;

    // The sum of all constituent 4-vectors
    rk::P4 jetSum() const;

    // The sum of 4-vectors which includes charged particles only
    rk::P4 chargedSum() const;

    // Comparison operators for sorting vectors of jets using std::sort.
    // Jets will be sorted by their transverse momentum.
    inline bool operator<(const SimulatedJet& r) const
    {
        return jetSum().pt() < r.jetSum().pt();
    }
    inline bool operator>(const SimulatedJet& r) const
    {
        return jetSum().pt() > r.jetSum().pt();
    }
};

// The "event" is a collection of jets. In addition,
// there are some convenience function for studying
// event properties.
struct SimpleEvent : public std::vector<SimulatedJet>
{
    // Find the 4-momentum of a jet closest to the given direction.
    // This function can be used to find event jets matching to
    // the reconstructed jets.
    rk::P4 closestJet(double eta, double phi) const;

    // Find another jet which has closest direction
    // to the jet with the given number
    unsigned closestJetNumber(unsigned thisJetNumber) const;

    // Return the "delta R" for the two closest jets in the event
    double minDr() const;
};

// The following function figures out the angle phi at which
// particle track hits the calorimeter situated at radius r
// away from the production point, assuming simple circular
// trajectory in magnetic field B and neglecting the energy
// losses along the way. The return value indicates whether
// the particle actually reaches the calorimeter.
bool phiAtRadius(const PythiaJetGun::Particle& particle, double B,
                 double r, double *phi);

#endif // SIMPLEEVENT_HH_
