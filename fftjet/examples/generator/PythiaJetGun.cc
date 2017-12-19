#include <cassert>
#include <cmath>
#include <cfloat>

#include "PythiaJetGun.hh"
#include "pythia6.h"

using namespace geom3;

namespace PythiaJetGun
{
    double Particle::charge() const
    {
        return pychge(code_)/3.0;
    }

    double Particle::gyrationRadius(const double B) const
    {
        const double qB = charge()*B;
        if (qB == 0.0)
            return DBL_MAX;
        else
        {
            const double c = 299792458.0;
            return 1.e9/c*p4_.pt()/qB;
        }
    }

    double pyrandom()
    {
        return pyr();
    }

    void pythiacmd(const char* cmd)
    {
        pygive(cmd);
    }

    double gaussRandom(const double mean, const double sigma)
    {
        double r;
        do {
            r = 1.0 - pyrandom();
        } while (r <= 0.0 || r > 1.0);
        return mean + sigma*sqrt(-2.0*log(r))*cos(2.0*M_PI*pyrandom());
    }

    void shoot1(const int code, const double pt, const double eta,
                const double phi, Particle* partonptr,
                std::vector<Particle>* output)
    {
        assert(pt > 0.0);
        assert(code >= 1 && code <= 500);
        int zcode = code;
        if (pyrandom() < 0.5)
            zcode = -code;

        rk::P4 parton(Vector3(pt*cos(phi), pt*sin(phi),
                              pt*sinh(eta)), pymass(zcode));

        // Generate the hadronization, shooting either particle
        // or antiparticle in the given direction
        py1ent(0, zcode, parton.e(), parton.theta(), phi);

        // Clear all decayed particles (this also removes neutrinos)
        pyedit(2);

        // Pick up all particles
        const int nparticles = pyjets_.n;
        if (nparticles > 0)
        {
            output->clear();
            output->reserve(nparticles);
            for (int i=0; i<nparticles; ++i)
            {
                Vector3 p(pyjets_.p[0][i],pyjets_.p[1][i],pyjets_.p[2][i]);
                output->push_back(Particle(pyjets_.k[1][i], 
                                           rk::P4(p, pyjets_.p[4][i])));
            }
            *partonptr = Particle(zcode, parton);
        }
        else
        {
            /* Hmm... Seems like something did not get generated. Retry. */
            shoot1(code, pt, eta, phi, partonptr, output);
        }
    }

    void shoot2(const int code, const double pt, const double eta,
                const double phi, Particle* partonptr,
                std::vector<Particle>* output)
    {
        assert(pt > 0.0);
        assert(code >= 1 && code <= 500);
        int zcode = code;
        if (pyrandom() < 0.5)
            zcode = -code;

        rk::P4 parton(Vector3(pt*cos(phi), pt*sin(phi),
                              pt*sinh(eta)), pymass(zcode));
        const UnitVector3 pdir(parton.momentum().direction());

        // Generate the hadronization, shooting either particle
        // or antiparticle in the Z direction
        py2ent(0, zcode, -zcode, 2.0*parton.e());

        // Clear all decayed particles (this also removes neutrinos)
        pyedit(2);

        // Create the rotation to the lab
        const UnitVector3 axis(
            UnitVector3::zAxis().cross(pdir).direction());
        const double angle = UnitVector3::zAxis().angle(pdir);
        Rotation3 rot(axis, angle);

        // Also, rotate the jet randomly around its own axis
        rot *= Rotation3(pdir, pyrandom()*2.0*M_PI);

        // Pick up only the particles going into positive Z direction
        const int nparticles = pyjets_.n;
        unsigned count = 0;
        for (int i=0; i<nparticles; ++i)
            if (pyjets_.p[2][i] > 0.0)
                ++count;

        if (count > 0)
        {
            output->clear();
            output->reserve(count);
            for (int i=0; i<nparticles; ++i)
                if (pyjets_.p[2][i] > 0.0)
                {
                    Vector3 p(pyjets_.p[0][i],pyjets_.p[1][i],pyjets_.p[2][i]);
                    output->push_back(Particle(pyjets_.k[1][i],
                                      rk::P4(rot.rotate(p), pyjets_.p[4][i])));
                }
            *partonptr = Particle(zcode, parton);
        }
        else
        {
            /* Hmm... Seems like something did not get generated. Retry. */
            shoot2(code, pt, eta, phi, partonptr, output);
        }
    }
}
