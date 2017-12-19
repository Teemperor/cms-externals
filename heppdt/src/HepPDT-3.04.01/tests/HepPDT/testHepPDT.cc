// ----------------------------------------------------------------------
// testHepPDT.cc
// Author: Lynn Garren
//
// test by reading the PDG table
// get filename and location of PDG table from input stream
//
// Usage:  testHepPDT 
//
// ----------------------------------------------------------------------

#include <fstream>
#include <iomanip>

#include "HepPDT/defs.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"
#include "HepPDT/HeavyIonUnknownID.hh"
// local include
#include "TestNuclearFragment.hh"

void pdtSimpleTest( char[300], std::ofstream & );
void pdtFragmentTest( char[300], std::ofstream & );
void duplicateFragmentTest( char[300], std::ofstream & );
void testPDMethods( HepPDT::ParticleDataTable&, std::ofstream & );

int main()
{
    char pdgfile[300] = "";
    const char outfile[] = "testHepPDT.out";
    std::cin >> pdgfile;
    // open output file
    std::ofstream wpdfile( outfile );
    if( !wpdfile ) { 
      std::cerr << "cannot open " << outfile << std::endl;
      exit(-1);
    }

    // construct a default PDT
    pdtSimpleTest( pdgfile, wpdfile );
    // now test the nuclear fragment option
    pdtFragmentTest( pdgfile, wpdfile );
    // check how we deal with duplicate fragments
    duplicateFragmentTest( pdgfile, wpdfile );

    return 0;
}

void pdtSimpleTest( char pdgfile[300], std::ofstream & wpdfile )
{
    // open input file
    std::ifstream pdfile( pdgfile );
    if( !pdfile ) { 
      std::cerr << "cannot open " << pdgfile << std::endl;
      exit(-1);
    }
    // construct empty PDT
    HepPDT::ParticleDataTable datacol( "2006 PDG Table" );
    {
        // Construct table builder
        HepPDT::TableBuilder  tb(datacol);
	// read the input - put as many here as you want
	if( !HepPDT::addPDGParticles( pdfile, tb ) ) { 
	    std::cout << "error reading PDG file " << std::endl; 
	}
    }	// the tb destructor fills datacol
    // done with pdfile, so close it
    pdfile.close();

    const char outfile1[] = "testHepPDTtable.out";
    // open output file
    std::ofstream wpdt1( outfile1 );
    if( !wpdt1 ) { 
      std::cerr << "cannot open " << outfile1 << std::endl;
      exit(-1);
    }
    datacol.writeParticleData(wpdt1);

    wpdfile << std::endl;

    // output some pion information
    HepPDT::ParticleData * pd;
    pd=datacol.particle(HepPDT::ParticleID(111));
    // test the ResonanceStructure cutoff methods here
    if(pd) { 
	pd->write(wpdfile);
	wpdfile << "Resonance info for 111 " 
        	<< pd->totalWidth().value() << " " 
        	<< pd->totalWidth().sigma() << " " 
		<< pd->lowerCutoff() << " " 
		<< pd->upperCutoff() << std::endl;
    }
    //  -111 is an illegal particle, no info will be written
    pd=datacol.particle(HepPDT::ParticleID(-111));
    if(pd) pd->write(wpdfile);
    pd=datacol.particle(HepPDT::ParticleID(211));
    if(pd) pd->write(wpdfile);
    // string lookup
    pd=datacol.particle(std::string("pi0"));
    if(pd) pd->write(wpdfile);

    // particle info
    datacol.writeParticleInfo(wpdfile);
       
    testPDMethods( datacol, wpdfile );
}

void pdtFragmentTest( char pdgfile[300], std::ofstream & wpdfile )
{
    wpdfile << std::endl;
    wpdfile << " Begin test of HeavyIonUnknownID " << std::endl;
    // reopen input file
    std::ifstream pdfile2( pdgfile );
    if( !pdfile2 ) { 
      std::cerr << "cannot open " << pdgfile << std::endl;
      exit(-1);
    }
    // construct another PDT instance that knows how to deal with unknown heavy ions
    // NOTE: normally you would construct a single ParticleDataTable with this option
    HepPDT::ParticleDataTable pdt2( "Handle Heavy Ions", 
                                    new HepPDT::HeavyIonUnknownID );
    {
        // Construct table builder
        HepPDT::TableBuilder  tb2(pdt2);
	// read the input - put as many here as you want
	if( !HepPDT::addPDGParticles( pdfile2, tb2 ) ) { 
	    std::cout << "error reading PDG file " << std::endl; 
	}
    }
    // done with pdfile, so close it
    pdfile2.close();

    // try a heavy ion
    HepPDT::ParticleData * pd=pdt2.particle(HepPDT::ParticleID(1000020040));
    wpdfile << " Printing information for unknown nuclear fragment " 
            << std::endl;
    if(pd) pd->write(wpdfile);
}

void duplicateFragmentTest( char pdgfile[300], std::ofstream & wpdfile )
{
    wpdfile << std::endl;
    wpdfile << " Begin test of duplicate nuclear fragments " << std::endl;
    // reopen input file
    std::ifstream pdfile2( pdgfile );
    if( !pdfile2 ) { 
      std::cerr << "cannot open " << pdgfile << std::endl;
      exit(-1);
    }
    // this test checks to see if we have actually added a fragment to the table
    HepPDT::ParticleDataTable pdt( "Duplicate Nuclear Fragments", 
                                    new HepPDT::TestNuclearFragment );
    {
        // Construct table builder
        HepPDT::TableBuilder  tb2(pdt);
	// read the input - put as many here as you want
	if( !HepPDT::addPDGParticles( pdfile2, tb2 ) ) { 
	    std::cout << "error reading PDG file " << std::endl; 
	}
    }
    // done with pdfile, so close it
    pdfile2.close();

    // try a heavy ion
    wpdfile << " Printing information for unknown nuclear fragments " 
            << std::endl;
    HepPDT::ParticleData * pd=pdt.particle(HepPDT::ParticleID(1000020040));
    if(pd) pd->write(wpdfile);
    pd=pdt.particle(HepPDT::ParticleID(1000020040));
    if(pd) pd->write(wpdfile);
    pd=pdt.particle(HepPDT::ParticleID(1000010040));
    if(pd) pd->write(wpdfile);
    pd=pdt.particle(HepPDT::ParticleID(1000020040));
    if(pd) pd->write(wpdfile);
    // what is the state of the table?
    const char outfile2[] = "testHepPDTfragment.out";
    std::ofstream wpdt( outfile2 );
    if( !wpdt ) { 
      std::cerr << "cannot open " << outfile2 << std::endl;
      exit(-1);
    }
    pdt.writeParticleData(wpdt);

    // check isStable
    const char outfile3[] = "testHepPDTstatus.out";
    std::ofstream wpdt3( outfile3 );
    if( !wpdt3 ) { 
      std::cerr << "cannot open " << outfile3 << std::endl;
      exit(-1);
    }
    pdt.writeParticleStatus(wpdt3);
}

void testPDMethods( HepPDT::ParticleDataTable& datacol, std::ofstream & wpdfile )
{
    wpdfile << std::endl;
    wpdfile << "Begin test of ParticleData methods " << std::endl;
    HepPDT::ParticleData * pd;
    int id[27] = { 5, 24, 15, 213, 3214, 10213, 9050225, 541, 129050225,
                   2000025, 3101, 3301, 2212, 1000020040, 1000060120, 555,
		   5000040, 5100005, 24, 5100024, 5100025, 9221132, 
                   4111370, -4120240, 4110050, 10013730, 1000612 };
    int it;
    for( it=0; it < 27; it++ ) {
        pd=datacol.particle(HepPDT::ParticleID(id[it]));
	if(pd) {
	    if( pd->hasUp() ) {
		wpdfile << "Particle " << pd->name()  << " " << pd->pid() 
        		  << " has an up quark" << std::endl;
	    }
	    if( pd->hasDown() ) {
		wpdfile << "Particle " << pd->name()  << " " << pd->pid() 
        		  << " has a down quark" << std::endl;
	    }
	    if( pd->hasStrange() ) {
		wpdfile << "Particle " << pd->name()  << " " << pd->pid() 
        		  << " has a strange quark" << std::endl;
	    }
	    if( pd->hasCharm() ) {
		wpdfile << "Particle " << pd->name() << " " << pd->pid() 
        		  << " has a charmed quark" << std::endl;
	    }
	    if( pd->hasBottom() ) {
		wpdfile << "Particle " << pd->name() << " " << pd->pid() 
        		  << " has a bottom quark" << std::endl;
	    }
	    if( pd->hasTop() ) {
		wpdfile << "Particle " << pd->name() << " " << pd->pid() 
        		  << " has a top quark" << std::endl;
	    }
	}
    }
}
