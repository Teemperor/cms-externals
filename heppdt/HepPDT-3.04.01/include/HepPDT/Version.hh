// ----------------------------------------------------------------------
//
// Version.hh
// Author:  Lynn Garren
//
//  for now, this is a free function
//
// ----------------------------------------------------------------------
#ifndef HEPPDTVERSION_HH
#define HEPPDTVERSION_HH

#include <string>
#include <iostream>

namespace HepPDT {

void version( );			//!< print HepPDT version
void writeVersion( std::ostream & os );	//!< write HepPDT version to os
std::string versionName( );	//!< return HepPDT version

}	// HepPDT

#endif // HEPPDTVERSION_HH
