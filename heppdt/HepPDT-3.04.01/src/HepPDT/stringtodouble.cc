// ----------------------------------------------------------------------
//
// stringtodouble.cc
// Author: Lynn Garren
//
// ----------------------------------------------------------------------

#include <sstream>

#include "HepPDT/defs.h"
#include "HepPDT/stringtodouble.hh"

namespace HepPDT {

double stringtodouble( std::string & numb )
{
    double dn;
    int sl = numb.length();
    std::istringstream snum( numb.substr(0,sl).c_str() );
    snum >> dn;
    return dn;
}

}  // namespace HepPDT
