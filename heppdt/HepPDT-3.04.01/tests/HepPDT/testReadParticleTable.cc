// $Id: testReadParticleTable.cc.in,v 1.6 2009/11/25 02:20:37 garren Exp $
// ----------------------------------------------------------------------
// testReadParticleTable.cc
//
// read particle.tbl and write it out
//
// ----------------------------------------------------------------------

#include <fstream>

#include "HepPDT/defs.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"

int main()
{
    const char infile[] = "../../data/particle.tbl";
    const char infile2[] = "../../tests/HepPDT/extras.tbl";
    const char outfile[] = "testReadParticleTable.out";
    // open input files
    std::ifstream pdfile( infile );
    if( !pdfile ) { 
      std::cerr << "cannot open " << infile << std::endl;
      exit(-1);
    }
    std::ifstream pdfile2( infile2 );
    if( !pdfile2 ) { 
      std::cerr << "cannot open " << infile2 << std::endl;
      exit(-1);
    }
    // construct empty PDT
    HepPDT::ParticleDataTable datacol( "Generic Particle Table" );
    {
        // Construct table builder
        HepPDT::TableBuilder  tb(datacol);
	// read the input - put as many here as you want
        // bool  addParticleTable( std::istream&, TableBuilder&,
        //                         bool validate = false );
	// where:  validate=true => verify that the ParticleID is valid
        if( !addParticleTable( pdfile, tb, true ) ) { 
	    std::cout << "error reading PDG pdt file " << std::endl; 
	}
        if( !addParticleTable( pdfile2, tb, true ) ) { 
	    std::cout << "error reading extra pdt file " << std::endl; 
	}
    }	// the tb destructor fills datacol
    // open the output stream
    std::ofstream wfile( outfile );
    if( !wfile ) { 
      std::cerr << "cannot open " << outfile << std::endl;
      exit(-1);
    }
    // write the data table
    datacol.writeParticleData(wfile);
    // try some heavy ions
    wfile << std::endl;
    wfile << std::endl;
    HepPDT::ParticleData * pd;
    pd=datacol.particle(HepPDT::ParticleID(1000020040));
    if(pd) pd->write(wfile);
    pd=datacol.particle(HepPDT::ParticleID(1000050110));
    if(pd) pd->write(wfile);

    // check isStable
    const char outfile3[] = "testReadParticleTableStatus.out";
    std::ofstream wpdt3( outfile3 );
    if( !wpdt3 ) { 
      std::cerr << "cannot open " << outfile3 << std::endl;
      exit(-1);
    }
    datacol.writeParticleStatus(wpdt3);

    return 0;
}
