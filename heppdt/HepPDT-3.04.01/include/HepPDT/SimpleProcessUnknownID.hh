#ifndef SIMPLEPROCESSUNKNOWNID_HH
#define SIMPLEPROCESSUNKNOWNID_HH
// ----------------------------------------------------------------------
//
// SimpleProcessUnknownID.hh
// Author: Lynn Garren
//
// ----------------------------------------------------------------------
#include "HepPDT/ProcessUnknownID.hh"
#include "HepPDT/ParticleData.hh"

namespace HepPDT {
  
//! The SimpleProcessUnknownID class inherits from ProcessUnknownID

///
/// \class SimpleProcessUnknownID
/// \author Lynn Garren
///
/// This provides the default HepPDT behaviour and will be invoked
/// by default unless you specify a different behaviour when creating
/// your ParticleDataTable object.
///
class SimpleProcessUnknownID : public ProcessUnknownID {
public:
  SimpleProcessUnknownID() {}

  virtual ParticleData  * processUnknownID
                      ( ParticleID, const ParticleDataTable & )
    { return 0; }
    
};

}	// HepPDT

#endif // SIMPLEPROCESSUNKNOWNID_HH
