// ----------------------------------------------------------------------
//
// addPythiaParticles.cc
// Author: Lynn Garren
//
// this has the functions used by addPythiaParticles
//
// ----------------------------------------------------------------------

#include <iostream>
#include <string>
#include <sstream>

#include "HepPDT/defs.h"
#include "HepPDT/TempParticleData.hh"
#include "HepPDT/TableBuilder.hh"

namespace HepPDT {

bool  addPythiaParticles( std::istream & pdfile, TableBuilder & tb )
{
  std::string pdline, aname;
  int id, kf;
  int saveid=0;
  int anti=0;
  // read and parse each line
  while( std::getline( pdfile, pdline) ) {
    if( detail::getPythiaid( kf, pdline ) ) {
      if( kf != 0 ) {
          // this is a new particle definition
          saveid = id = HepPID::translatePythiatoPDT( kf );
          TempParticleData& tpd = tb.getParticleData( ParticleID( id ) );
          detail::parsePythiaLine( tpd, anti, aname, pdline );
	  if( anti > 0 ) {
	      // code here to define antiparticles
	      TempParticleData& atpd = tb.getAntiParticle( ParticleID( id ), aname );
              // use this variable (fake out the compiler)
              atpd.tempSource = tpd.tempSource;
              atpd.tempOriginalID = -tpd.tempOriginalID;
              atpd.tempMass = tpd.tempMass;
	  }
      } else if( saveid != 0 ) {
          TempParticleData& tpd = tb.getParticleData( ParticleID( saveid ) );
          detail::parsePythiaDecayLine( tpd, pdline );
	  if( anti > 0 ) {
	      // code here to append antiparticle decays
	  }
      }
    }
  }
  std::cout << "found " << tb.size() << " particles" << std::endl;
  return true;
}

namespace detail {
void parsePythiaLine( TempParticleData & tpd, int & anti, std::string & aname, const std::string & pdline )
{
    // this line defines a particle
    std::string name1;
    int kf, kc, chg, col, decay, blank;
    double mass, width, wcut, lifet;
    
    // check for valid TempParticleData
    anti = 0;
    int sl = pdline.length();
    if( tpd.tempID.pid() == 0 ) { return; }
    // have a valid PID, so proceed
    std::istringstream pids( pdline.substr(0,17).c_str() );
    pids >> kf >> kc;		// must access this before getting another stream
    // unfortunately, the istrstream trick does not work with strings in KCC
    aname = "";
    //std::istringstream pnames( pdline.substr(21,32).c_str() );
    //pnames >> name1 >> aname;
    blank = pdline.substr(21,16).find(" ");
    name1 = pdline.substr(21,blank);
    blank = pdline.substr(37,16).find(" ");
    aname = pdline.substr(37,blank);
    std::istringstream particle( pdline.substr(54,sl-53).c_str() );
    particle >> chg >> col >> anti >> mass >> width >> wcut >> lifet >> decay ;
    // std::cout << kf << " " << kc << " " << tpd.tempID.pid() << " " 
    // 	    << name1 << " " << aname << " " << chg 
    // 	    << " " << mass << " " << width << std::endl; 
    tpd.tempParticleName = name1;
    tpd.tempSource = "Pythia";
    tpd.tempOriginalID = kf;
    tpd.tempCharge = double(chg)/3.0;
    tpd.tempColorCharge = col;
    tpd.tempMass = Measurement( mass, 0.0 );
    tpd.tempHighCutoff = wcut;
    // either width or lifetime is defined - not both
    if( width > 0. ) {
        tpd.tempWidth = Measurement( width, 0.0 );
    } else if( lifet > 0. ) { 
        tpd.tempWidth = Measurement( calculateWidthFromLifetime( lifet ), 0.0 );
    }
}

void parsePythiaDecayLine( TempParticleData & tpd, const std::string & )
{
    // check for valid TempParticleData
    if( tpd.tempID.pid() == 0 ) { return; }
}

}  // namespace detail

}	// HepPDT
