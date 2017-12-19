// ----------------------------------------------------------------------
//
// addQQParticles.cc
// Author: Lynn Garren
//
// this has the functions used by addQQParticles
//
// ----------------------------------------------------------------------

#include <sstream>
#include <string>
#include <iostream>

#include "HepPDT/defs.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/TempParticleData.hh"
#include "HepPID/ParticleIDTranslations.hh"

namespace HepPDT {

// helper functions
bool getQQLineType( std::string & ltype, int & id, std::string & name, 
                     const std::string & pdline );
bool parseQQDecayLine( const std::string & pdline );
void parseQQParticle( TempParticleData & tpd, const std::string & pdline );

bool  addQQParticles( std::istream & pdfile, TableBuilder & tb )
{ 
  std::string pdline, ltype, def, name, cname;
  int id, eid;
  std::string genName = "QQ";
  // read and parse each line
  while( std::getline( pdfile, pdline) ) {
    if( getQQLineType( ltype, eid, name, pdline ) ) {
      if( ltype == "PARTICLE" ) {
          // this is a new particle definition
          id = HepPID::translateQQtoPDT( eid );
	  if( id != 0 ) { 
              TempParticleData& tpd = tb.getParticleData( ParticleID( id ) );
              parseQQParticle( tpd, pdline );
              tpd.tempSource = genName;
	      tb.addParticle( tpd );
	  }
      } else if( ltype=="QQBAR" ) {
          // these are psuedo-particles used in the decay tables
          id = HepPID::translateQQbar( eid );
          TempParticleData& tpd = tb.getParticleData( ParticleID( id ) );
          tpd.tempParticleName = name;
          tpd.tempOriginalID = eid;
          tpd.tempSource = genName;
	  tb.addParticle( tpd );
      } else if( ltype=="VERSION" ) {
          genName += ' ';
	  genName += name;
      } else if( ltype=="HIDE" ) {
      } else if( ltype=="PDG" ) {
      } else if( ltype=="PARITY" ) {
      } else if( ltype=="DECAY" ) {
          // QQ assumes that all particles have been defined before
          // they are used in the decay table
	  if( tb.hasParticleData( name ) ) {
	      // read lines until we reach an "Enddecay"
	      while( parseQQDecayLine( pdline ) ) {
        	  std::getline( pdfile, pdline);
	      }
	  } else {
	      std::cout << "HepPDT TableBuilder: could not match QQ name " << name << std::endl;
	  }
      }
    }
  }
  std::cout << "found " << tb.size() << " particles" << std::endl;
    return true;
}

bool getQQLineType( std::string & ltype, int & id, std::string & name, 
                    const std::string & pdline )
{
    std::string s1, s2;
    int sl = pdline.length();
    std::string firstc = pdline.substr(0,1);
    id = 0;
    name = "";
    ltype = "";
    // check for comments
    if( firstc == "*" || firstc == ";" ) { return false; }
    // a ; may be anywhere on the line
    // check for empty lines
    if( sl > 5 ){
	std::istringstream thisline( pdline.substr(0,sl).c_str() );
	thisline >> ltype;
	if( ltype == "PARTICLE" ) {
            thisline >> name >>  id;
        } else if( ltype=="QQBAR" ) {
            thisline >> name >>  id;
        } else if( ltype=="PDG" ) {
            thisline >> name >>  id;
        } else if( ltype=="PARITY" ) {
            thisline >> name >>  id;
        } else if( ltype=="CPARITY" ) {
            thisline >> name >>  id;
	} else if( ltype == "HIDE" ) {
            thisline >> name;
	} else if( ltype == "DECAY" ) {
            thisline >> name;
	} else if( ltype == "MIXING" ) {
            thisline >> name;
	} else if( ltype == "VERSION" ) {
            thisline >> name;
	}
	return true;
    } else {
        return false;
    }
}

void parseQQParticle( TempParticleData & tpd, const std::string & pdline )
{
    // this line defines a particle
    std::string s1, s2, ltype, name;
    int id, sid;
    double chg, spin;
    double mass, width, mmax, lifet, mmin;

    // check for valid TempParticleData
    if( tpd.tempID.pid() == 0 ) { return; }
    // have a valid PID, so proceed
    int sl = pdline.length();
    std::istringstream thisline( pdline.substr(0,sl).c_str() );
    thisline >> ltype >> name >> id >> sid >>  mass >> chg >> spin >> lifet;
    width = mmin = mmax = 0.;
    if( !thisline.eof() ) {
        thisline >> width >> mmin >> mmax; 
    }
    if( ltype != "PARTICLE" ) { 
        std::cout << "called parseQQParticle with wrong line type" << std::endl;
        return; 
    }
    tpd.tempParticleName = name;
    tpd.tempOriginalID = id;
    tpd.tempCharge = chg;
    tpd.tempMass = Measurement( mass, 0.0 );
    tpd.tempLowCutoff = mmin;
    tpd.tempHighCutoff = mmax;
    if( width > 0. ) {	// by default, only lifetime is specified
        tpd.tempWidth = Measurement( width, 0.0 );
    } else {
        tpd.tempWidth = Measurement( calculateWidthFromLifetime(lifet), 0.0 );
    }
    if( tpd.tempSpin.totalSpin() != spin ) {
         //std::cout << "QQ spin " << spin << " does not match ParticleID spin " 
	 //     << tpd.tempSpin.totalSpin() 
	 //     << " for " << name << std::endl;
	tpd.tempSpin.setTotalSpin( spin );
    }
}

bool parseQQDecayLine( const std::string & pdline)
{
    std::string ltype;
    int sl = pdline.length();
    std::istringstream thisline( pdline.substr(0,sl).c_str() );
    // better check for comments and blank lines.....
    if( sl < 5 ) { return true; }
    std::string firstc = pdline.substr(0,1);
    if( firstc == "*" || firstc == ";" ) { return true; }
    // check line type
    thisline >> ltype;
    if( ltype == "DECAY" ) {
        // begin decay
	return true;
    } else if( ltype == "ENDDECAY" ) {
        // end decay
	return false;
    }
    // this is a good decay line
    // accumulate information in temporary QQ structures
    return true;
}


}	// HepPDT
