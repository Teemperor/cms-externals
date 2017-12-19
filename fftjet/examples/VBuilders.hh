//=========================================================================
// VBuilders.hh
//
// Functors for building various quantities from grid points. They use
// the "kinematics" package for 3- and 4-vector implementations. They are
// intended for use with "KernelRecombinationAlg" or similar templated
// classes.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef VBUILDERS_HH_
#define VBUILDERS_HH_

#include <cmath>

#include "rk.hh"

struct PtEtaP4Builder
{
    inline rk::P4 operator()(double pt, double eta, double phi) const
    {
        return rk::P4(geom3::Vector3(cos(phi), sin(phi),
                                     sinh(eta)), 0.0)*pt;
    }
};

struct EnergyEtaP4Builder
{
    inline rk::P4 operator()(double energy, double eta, double phi) const
    {
        // There is no mass associated with this energy... We will
        // assume that the mass is 0 and proceed as if the energy
        // is the momentum.
        const double pt = energy/cosh(eta);
        return rk::P4(geom3::Vector3(cos(phi), sin(phi),
                                     sinh(eta)), 0.0)*pt;
    }
};

struct CentroidBuilder
{
    inline geom3::Vector3 operator()(double mag, double eta, double phi) const
    {
        return geom3::Vector3(mag*eta, mag*phi, mag);
    }
};

#endif // VBUILDERS_HH_
