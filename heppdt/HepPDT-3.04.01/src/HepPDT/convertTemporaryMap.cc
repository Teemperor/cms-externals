// ----------------------------------------------------------------------
//
// convertTemporaryMap.cc
// Author: Lynn Garren, Walter Brown
//
// ----------------------------------------------------------------------

#include "HepPDT/defs.h"
#include "HepPDT/ParticleID.hh"
#include "HepPDT/ParticleDataTable.hh"
#include "HepPDT/ParticleData.hh"
#include "HepPDT/TempParticleData.hh"
#include "HepPDT/Version.hh"

namespace HepPDT {

void ParticleDataTable::convertTemporaryMap( TempMap & tempPDT, std::ostream & )
{
  std::map<ParticleID,TempParticleData>::iterator cit;
  // first loop and create a PD map with empty DecayData
  for( cit=tempPDT.begin(); cit != tempPDT.end(); ++cit ) {
    TempParticleData & tempd = cit->second;
    // create PD and add it to the list
    ParticleData pd( tempd );
    // now create ParticleData and add it to the map
    addParticle( pd );
  }
}

}	// HepPDT
