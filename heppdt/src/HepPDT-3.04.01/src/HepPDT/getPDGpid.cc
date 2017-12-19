// ----------------------------------------------------------------------
//
// getPDGpid.cc
// Author:  Lynn Garren
//
// parse a line in the PDG table and return the particle ID numbers
//
// ----------------------------------------------------------------------

#include <string>
#include <vector>
#include <sstream>

#include "HepPDT/defs.h"
#include "HepPDT/TableBuilder.hh"

namespace HepPDT {

namespace detail {
void getPDGpid( std::vector<int> & idlist,  std::string & pdline )
{
    int sl, cl, id1, id2, id3, id4;
    sl = pdline.length();
    cl = pdline.find('*');
    // make sure the vector is empty
    idlist.clear();
    if( cl != 0 && sl > 80 ){
      // this is a valid line, so parse it
      // initialize possible ID's to zero
      id1 = id2 = id3 = id4 = 0;
      std::istringstream idnet( pdline.substr(1,32).c_str() );
      idnet >> id1 >> id2 >> id3 >> id4;
      if( id1 > 0 ) { idlist.push_back( id1 ); }
      if( id2 > 0 ) { idlist.push_back( id2 ); }
      if( id3 > 0 ) { idlist.push_back( id3 ); }
      if( id4 > 0 ) { idlist.push_back( id4 ); }
    }
}

void getPDGnames( std::vector<std::string> & namelst, std::string & pdline )
{
    int sl, cl;
    std::string name, charges, fullname;
    sl = pdline.length();
    cl = pdline.find('*');
    // make sure the vector is empty
    namelst.clear();
    if( cl != 0 && sl > 80 ){
       // this is a valid line, so parse it
       std::istringstream namelist( pdline.substr(68,21).c_str() );
       namelist >> name >> charges;
       // parse the charge list and add charge to base name
       char buf[20];
       std::istringstream chglst(charges);
       while( chglst ) {
	  chglst.getline(buf,20,',');
	  if( chglst.gcount() > 0 ) {
	      fullname = name + buf;
	      namelst.push_back( fullname );
	  }
       }
    }
}


}  // namespace detail

} // namespace HepPDT

