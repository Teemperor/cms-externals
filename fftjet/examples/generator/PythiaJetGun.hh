// Convenient wrapper for two simple PYTHIA jet production mechanisms

#ifndef PYTHIAJETGUN_HH_
#define PYTHIAJETGUN_HH_

#include <vector>

#include "rk.hh"

namespace PythiaJetGun
{
    // A simple particle class which carries the 4-momentum
    // of the particle and its PYTHIA particle code.
    // Jets will be made of collections of such particles.
    class Particle
    {
    public:
        inline Particle() : p4_(), code_(0) {}
        inline Particle(int code, const rk::P4& p4) : p4_(p4), code_(code) {}

        inline const rk::P4& p4() const {return p4_;}
        inline int code() const {return code_;}
        double charge() const;

        // The gyration radius is returned in meters for magnetic field B
        // in Tesla and particle momentum in GeV/c. The radius is positive
        // if signs of B and particle charge are the same, and negative
        // otherwise. DBL_MAX is returned for neutral particles or in case
        // B is 0.
        double gyrationRadius(double B) const;

    private:
        rk::P4 p4_;
        int code_;
    };

    // The following function uses the original Pythia single jet gun
    // function PY1ENT. Such a jet is not required to conserve
    // energy, momentum, or flavour.
    void shoot1(int partonCode, double pt, double eta, double phi,
                Particle* parton, std::vector<Particle>* jetParticles);

    // The following creates a q qbar pair and takes a jet from the
    // hemisphere of q or qbar (choice of q or qbar is random).
    // Chosen parton and jet particles are returned.
    void shoot2(int partonCode, double pt, double eta, double phi,
                Particle* parton, std::vector<Particle>* jetParticles);

    // Random number generators based on Pythia "pyr" function
    double pyrandom();
    double gaussRandom(double mean, double sigma);

    // Wrapper for "pygive"
    void pythiacmd(const char* cmd);
}

#endif // PYTHIAJETGUN_HH_
