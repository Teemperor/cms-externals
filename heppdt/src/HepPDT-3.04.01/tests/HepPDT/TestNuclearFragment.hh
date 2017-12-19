#ifndef TestNuclearFragment_HH
#define TestNuclearFragment_HH
// ----------------------------------------------------------------------
//
// TestNuclearFragment.hh
// Author: Lynn Garren
//
// ----------------------------------------------------------------------

#include <sstream>

#include "HepPDT/ProcessUnknownID.hh"
#include "HepPDT/ParticleData.hh"

namespace HepPDT {

//! The TestNuclearFragment class inherits from ProcessUnknownID

///
/// \class TestNuclearFragment
/// \author Lynn Garren
///
/// Create and return a pointer to a ParticleData object if 
/// invoked by a valid nuclear fragment.
/// TestNuclearFragment::processUnknownID( ParticleID ) calculates 
/// both charge and an approximate mass for the particle.
/// This implementation also sends output to cout for code checking
///
class TestNuclearFragment : public ProcessUnknownID {
public:
  TestNuclearFragment() {}

  virtual ParticleData  * processUnknownID( ParticleID,  const ParticleDataTable & pdt );
    
};


inline ParticleData * TestNuclearFragment::processUnknownID
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
    std::ostringstream sname;
    sname << "fragment-" << key.pid();
    tpd.tempParticleName = sname.str();
    std::cout << "TestNuclearFragment::processUnknownID: creating ParticleData for " 
	      << key.pid() << std::endl;
    return new ParticleData(tpd);
}

}	// HepPDT

#endif // TestNuclearFragment_HH
