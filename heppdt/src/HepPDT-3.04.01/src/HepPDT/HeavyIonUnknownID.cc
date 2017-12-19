// ----------------------------------------------------------------------
//
// HeavyIonUnknownID.cc
// Author: Lynn Garren
//
// ----------------------------------------------------------------------
#include "HepPDT/HeavyIonUnknownID.hh"
#include "HepPDT/ParticleDataTable.hh"

namespace HepPDT {

ParticleData * HeavyIonUnknownID::processUnknownID
              ( ParticleID key, const ParticleDataTable & pdt ) 
{ 
    if( ! key.isNucleus() ) return 0;
    // have to create a TempParticleData with all properties first
    TempParticleData tpd(key);
    // calculate approximate mass
    // WARNING: any calls to particle() from here MUST reference 
    //          a ParticleData which is already in the table
    // This convention is enforced.
    const ParticleData * proton = pdt.particle(2212);
    if( ! proton ) return 0;
    double protonMass = proton->mass();
    tpd.tempMass = Measurement(key.A()*protonMass, 0.);
    return new ParticleData(tpd);
}

}	// HepPDT
