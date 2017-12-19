// ----------------------------------------------------------------------
//
// addIsajetParticles.cc
// Author: Lynn Garren
//
// ----------------------------------------------------------------------

#include "HepPDT/defs.h"
#include "HepPDT/TempParticleData.hh"
#include "HepPDT/TableBuilder.hh"

namespace HepPDT {

bool  addIsajetParticles( std::istream & pdfile, TableBuilder & tb )
{
  // this function should only be called on isaparticles.dat, not on isadecay.dat
  std::string pdline;
  int id, isaid;
  // read and parse each line
  while( std::getline( pdfile, pdline) ) {
    if( detail::getIsajetID( isaid, pdline ) ) {
          // this is a new particle definition
          id = HepPID::translateIsajettoPDT( isaid );
          TempParticleData& tpd = tb.getParticleData( ParticleID( id ) );
          detail::parseIsajetLine( tpd, pdline );
    }
  }
  std::cout << "found " << tb.size() << " particles" << std::endl;
  return true;
}

namespace detail {
void parseIsajetLine( TempParticleData & tpd, const std::string & pdline )
{
    // this line defines a particle
    std::string isaname;
    int isaid, fl1, fl2, fl3, spin, idx, blank;
    double mass, chg;
    
    // check for valid TempParticleData
    int sl = pdline.length();
    if( tpd.tempID.pid() == 0 ) { return; }
    // have a valid PID, so proceed
    std::istringstream pids( pdline.substr(0,17).c_str() );
    pids >> isaid ;		// must access this before getting another stream
    // unfortunately, the istrstream trick does not work with strings in KCC
    //std::istringstream pnames( pdline.substr(21,32).c_str() );
    //pnames >> isaname >> aname;
    blank = pdline.substr(11,10).find(" ");
    isaname = pdline.substr(11,blank);
    std::istringstream particle( pdline.substr(21,sl-20).c_str() );
    particle >> mass >> chg >> fl1 >> fl2 >> fl3 >> spin >> idx ;
    // std::cout << isaid << " " << tpd.tempID.pid() << " " 
    // 	    << isaname << " " << chg 
    // 	    << " " << mass << " " << fl1 << " " << fl2 << " " << fl3
    //	    << " " << spin << " " << idx << std::endl; 
    tpd.tempParticleName = isaname;
    tpd.tempCharge = chg;
    tpd.tempSpin = SpinState(0.,double(spin),0.);
    tpd.tempMass = Measurement( mass, 0.0 );
}
}  // namespace detail

}  // namespace HepPDT
