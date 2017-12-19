// ----------------------------------------------------------------------
//
// addParticleTable.cc
// Author: Lynn Garren
//
// this has the functions used by addParticleTable
//
// ----------------------------------------------------------------------

#include <iostream>
#include <string>
#include <sstream>

#include "HepPDT/defs.h"
#include "HepPDT/TempParticleData.hh"
#include "HepPDT/TableBuilder.hh"
#include "HepPID/ParticleName.hh"

namespace HepPDT {

bool  addParticleTable( std::istream & pdfile, TableBuilder & tb,
                        bool validate )
{
  //    validate => verify that the ParticleID is valid
  std::string pdline, aname;
  int id;
  // read and parse each line
  while( std::getline( pdfile, pdline) ) {
    if( detail::getParticleID( id, pdline ) ) {
      ParticleID pid( id );
      if( validate ) {
	if( pid.isValid() ) {
          // this is a new particle definition
          TempParticleData& tpd = tb.getParticleData( pid );
          detail::parseParticleLine( tpd, pdline );
	}
      } else {
        // this is a new particle definition
        TempParticleData& tpd = tb.getParticleData( pid );
        detail::parseParticleLine( tpd, pdline );
      }
    }
  }
  std::cout << "found " << tb.size() << " particles" << std::endl;
  return true;
}


namespace detail {
bool getParticleID( int & id, const std::string & pdline )
{
    int sl = pdline.length();
    id = 0;
    // line is too short
    if( sl < 30 ) return false;
    // now check for possible comments 
    std::string firstc = pdline.substr(0,1);
    if( firstc == "#" ) return false;
    std::string first2c = pdline.substr(0,2);
    if( first2c == "//" ) return false;
    // hope that this is now a valid line
    //std::istringstream var1(pdline.substr(21,12).c_str());
    std::istringstream var1(pdline.c_str());
    var1 >> id;
    if( id == 0 ) return false;
    // have non-zero ID
    return true;
}

void parseParticleLine( TempParticleData & tpd, const std::string & pdline )
{
    // this line defines a particle
    std::string name1;
    int id, chg;
    double mass, width, lifet;
    
    // check for valid TempParticleData
    if( tpd.tempID.pid() == 0 ) { return; }
    // have a valid PID, so proceed
    std::istringstream particle( pdline.c_str() );
    particle >> id >> name1 >> chg >> mass >> width >> lifet ;
    // allow for Q-balls
    if( tpd.tempID.isQBall() ) {
       // 10x the charge
       tpd.tempCharge = double(chg)/10.0;
    } else {
       // 3x the charge
       tpd.tempCharge = double(chg)/3.0;
    }
    tpd.tempParticleName = name1;
    tpd.tempSource = "ParticleTable";
    tpd.tempOriginalID = id;
    tpd.tempMass = Measurement( mass, 0.0 );
    // either width or lifetime is defined - not both
    if( width > 0. ) {
        tpd.tempWidth = Measurement( width, 0.0 );
    } else if( width == -1. ) { 
        tpd.tempWidth = Measurement( -1., 0.0 );
    } else if( lifet > 0. ) { 
        tpd.tempWidth = Measurement( calculateWidthFromLifetime( lifet ), 0.0 );
    } else {
        tpd.tempWidth = Measurement( 0.0, 0.0 );
    }
}
}  // namespace detail

}	// HepPDT
