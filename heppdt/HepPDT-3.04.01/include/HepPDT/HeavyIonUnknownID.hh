#ifndef HEAVYIONUNKNOWNID_HH
#define HEAVYIONUNKNOWNID_HH
// ----------------------------------------------------------------------
//
// HeavyIonUnknownID.hh
// Author: Lynn Garren
//
// ----------------------------------------------------------------------
#include "HepPDT/ProcessUnknownID.hh"
#include "HepPDT/ParticleData.hh"

namespace HepPDT {

//! The HeavyIonUnknownID class inherits from ProcessUnknownID

///
/// \class HeavyIonUnknownID
/// \author Lynn Garren
///
/// Create and return a pointer to a ParticleData object if 
/// invoked by a valid Heavy Ion.
/// HeavyIonUnknownID::processUnknownID( ParticleID ) calculates 
/// both charge and an approximate mass for the particle.
///
class HeavyIonUnknownID : public ProcessUnknownID {
public:
  HeavyIonUnknownID() {}

  virtual ParticleData  * processUnknownID( ParticleID,  const ParticleDataTable & pdt );
    
};

}	// HepPDT

#endif // HEAVYIONUNKNOWNID_HH
