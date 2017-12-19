// ----------------------------------------------------------------------
// listIsajetTranslation.cc
// Author: Lynn Garren
//
// Usage:  listIsajetTranslation 
//
// ----------------------------------------------------------------------

#include <fstream>
#include <iostream>
#include <cstdlib>     // for exit

#include "HepPID/ParticleIDTranslations.hh"

int main()
{
    const char outfile[] = "listIsajetTranslation.out";
    // open the output file
    std::ofstream wpdfile( outfile );
    if( !wpdfile ) { 
      std::cerr << "cannot open " << outfile << std::endl;
      exit(-1);
    }
    // write the particle names
    HepPID::writeIsajetTranslation( wpdfile );
}
