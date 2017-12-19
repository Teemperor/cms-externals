// ----------------------------------------------------------------------
//
// write.cc
// Author: Lynn Garren
//
// ----------------------------------------------------------------------

#include <iomanip>
#include <sstream>

#include "HepPDT/defs.h"
#include "HepPDT/ParticleData.hh"

using namespace std;

namespace HepPDT {

void ParticleData::write( std::ostream & ostr ) const
{
    if( !itsID.isValid() ) {
        std::cout << " WARNING: invalid particle " <<  itsID.pid() << std::endl;
        return;
    }
    Measurement mass = itsResonance.mass();
    Measurement width = itsResonance.totalWidth();
    Measurement ctau = lifetime();
    int il = itsParticleName.length();
    ostr << itsParticleName ;
    // pad the name
    int pad;
    for( pad = 0; pad < 20 - il; ++pad ) {
        ostr << " " ;
    }
    ostr.width(12);
    ostr <<  itsID.pid();
    ostr.width(8);
    ostr.precision(4);
    ostr << itsCharge;
    ostr.width(5);
    ostr.precision(2);
    ostr << " " << itsColorCharge;
    ostr.width(5);
    ostr << itsSpin.totalSpin() ;
    ostr.width(5);
    ostr << itsSpin.spin() ;
    ostr.width(5);
    ostr << itsSpin.orbAngMom() ;
    ostr.width(14);
    ostr.precision(5);
    ostr.setf(std::ios::scientific, std::ios::floatfield);
    ostr << mass.value() << "+/-" << mass.sigma() ;
    ostr.width(13);
    ostr.precision(4);
    ostr << width.value() << "+/-" << width.sigma() ;
    ostr.width(13);
    ostr.precision(4);
    ostr << lowerCutoff() << " " << upperCutoff();
    ostr.width(15);
    ostr.precision(4);
    ostr << ctau.value() << "+/-" << ctau.sigma();
    // return to default settings
    ostr.precision(6);
    ostr.setf(std::ios::fmtflags(0), std::ios::floatfield);
    // end the line
    ostr << std::endl;
}

void ParticleData::writeParticleInfo( std::ostream & os ) const
{
    write(os);
}

void ParticleData::writeParticleTranslation( std::ostream & os ) const
{
    os << " " << itsSource << ": ";
    os.width(20);
    os << itsParticleName;
    os.width(12);
    os << itsOriginalID;
    os << "   HepPDT: ";
    os.width(20);
    os << PDTname();
    os.width(12);
    os << pid() << std::endl;
    
    return;
}

}	// HepPDT
