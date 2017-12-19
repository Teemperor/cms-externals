// $Id: testReadEvtGen.cc.in,v 1.5 2009/01/09 20:28:34 garren Exp $
// ----------------------------------------------------------------------
// testReadEvtGen.cc
//
// read EvtGen table and write it out
//
// ----------------------------------------------------------------------

#include <fstream>

#include "HepPDT/defs.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"

int main()
{
    const char infile1[] = "../../examples/data/evt.pdl";
    const char infile2[] = "../../examples/data/DECAY.DEC";
    const char outfile[] = "testReadEvtGen.out";
    // open input files
    std::ifstream pdfile1( infile1 );
    if( !pdfile1 ) { 
      std::cerr << "cannot open " << infile1 << std::endl;
      exit(-1);
    }
    // construct empty PDT
    std::ifstream pdfile2( infile2 );
    if( !pdfile2 ) { 
      std::cerr << "cannot open " << infile2 << std::endl;
      exit(-1);
    }
    HepPDT::ParticleDataTable datacol( "EvtGen Table" );
    {
        // Construct table builder
        HepPDT::TableBuilder  tb(datacol);
	// read the input - put as many here as you want
        if( !addEvtGenParticles( pdfile1, tb ) ) { std::cout << "error reading EvtGen pdt file " << std::endl; }
        if( !addEvtGenParticles( pdfile2, tb ) ) { std::cout << "error reading EvtGen decay file " << std::endl; }
    }	// the tb destructor fills datacol
    std::ofstream wfile( outfile );
    if( !wfile ) { 
      std::cerr << "cannot open " << outfile << std::endl;
      exit(-1);
    }
    datacol.writeParticleData(wfile);
   
    return 0;
}