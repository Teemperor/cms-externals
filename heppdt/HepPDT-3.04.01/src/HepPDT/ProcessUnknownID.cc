// ----------------------------------------------------------------------
//
// ProcessUnknownID.cc
// Author: Lynn Garren
//
// ----------------------------------------------------------------------
#include "HepPDT/ProcessUnknownID.hh"
#include "HepPDT/ParticleDataTable.hh"

namespace HepPDT {

ParticleData * ProcessUnknownID::callProcessUnknownID
              ( ParticleID key, const ParticleDataTable & pdt ) 
{ 
    if ( alreadyHere ) return 0;
    alreadyHere = true;
    ParticleData * pd = processUnknownID( key, pdt );
    alreadyHere = false;
    return pd;
}

}	// HepPDT
