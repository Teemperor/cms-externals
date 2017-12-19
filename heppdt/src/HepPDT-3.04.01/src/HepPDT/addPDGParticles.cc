// ----------------------------------------------------------------------
//
// addPDGParticles.cc
// Author:  Lynn Garren
//
// this has the functions used by addPDGParticles
//
// ----------------------------------------------------------------------

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>

#include "HepPDT/defs.h"
#include "HepPDT/TempParticleData.hh"
#include "HepPDT/TableBuilder.hh"

namespace HepPDT {

bool  addPDGParticles( std::istream & pdfile, TableBuilder & tb )
{
  // mass and width lines can be in any order
  std::vector<int> idlist;
  std::vector<std::string> names;
  std::string pdline;
  // read and parse each line
  while( std::getline( pdfile, pdline) ) {
    detail::getPDGpid( idlist, pdline );
    detail::getPDGnames( names, pdline );
    if ( idlist.size() != names.size() ) {
        std::cout << "addPDGParticles ERROR: cannot make sense of line:" << std::endl;
	std::cout << "   " << pdline << std::endl;
    } else {
	for( unsigned int i = 0; i < idlist.size(); ++i )
	{
            TempParticleData& tpd = tb.getParticleData( ParticleID( idlist[i] ) );
            tpd.tempSource = "PDG table";
            tpd.tempOriginalID = idlist[i];
	    tpd.tempParticleName = names[i];
            detail::parsePDGline( tpd, pdline );
	}
    }
  }
  std::cout << "found " << tb.size() << " particles" << std::endl;
  return true;
}

namespace detail {
void parsePDGline( TempParticleData & tpd,  std::string & pdline )
{
   double v, e1, e2, err;
   std::string name, ckey, charges, fullname, chg;
   //std::string bigname;
   //int sl = pdline.length() - 1;		// <cr> at ends of lines
   // we already know that this is a valid line
   ckey = pdline.substr(0,1);
   name = charges = fullname = "";
   v = e1 = e2 = 0.0;
   std::istringstream val( pdline.substr(34,33).c_str() );
   val >> v >> e1 >> e2;
   err = sqrt( (e1*e1 + e2*e2)/2.0 );
   CheckPDGEntry( tpd, ckey, v, err );
}

bool CheckPDGEntry( TempParticleData & tpd, const std::string & ckey, 
                    double val, double err  )
{
   // now add appropriate properties
   if( ckey.find("M") == 0 ) {
      tpd.tempMass = Measurement( val, err );
   } else if( ckey.find("W") == 0 ) {
      tpd.tempWidth = Measurement( val, err );
   } else {
     std::cout << "unrecognized key " << ckey << std::endl;
     return false;
   }
   return true;
}
}  // namespace detail

}	// HepPDT
