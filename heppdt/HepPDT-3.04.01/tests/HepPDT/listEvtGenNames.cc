// ----------------------------------------------------------------------
// listEvtGenNames.cc
// Author: Lynn Garren
//
// read EvtGen table and write out translation from EvtGen to HepPDT
//
// ----------------------------------------------------------------------

#include <fstream>
#include <iostream>

#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"

int main()
{
    const char infile1[] = "../../examples/data/evt.pdl";
    const char outfile[] = "listEvtGenNames.out";
    // open input files
    std::ifstream pdfile1( infile1 );
    if( !pdfile1 ) { 
      std::cerr << "cannot open " << infile1 << std::endl;
      exit(-1);
    }
    // construct PDT
    HepPDT::ParticleDataTable datacol( "EvtGen Table" );
    {
        // Construct table builder
        HepPDT::TableBuilder  tb(datacol);
	// read the input - put as many here as you want
        if( !addEvtGenParticles( pdfile1, tb ) ) { std::cout << "error reading EvtGen pdt file " << std::endl; }
    }	// the tb destructor fills datacol
    // open output file
    std::ofstream wpdfile( outfile );
    if( !wpdfile ) { 
      std::cerr << "cannot open " << outfile << std::endl;
      exit(-1);
    }
    // write a translation list
    datacol.writeParticleTranslation( wpdfile );

    return 0;
}
