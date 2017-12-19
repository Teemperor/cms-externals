// -------------------------------------------------------------------
// examListIsajet.cc
// 
// ----------------------------------------------------------------------
///
/// \author Lynn Garren
///
/// List the isajet particle ID translations
/// Get ID list directly from isajet
///
/// When mixing C++ and Fortran, the main program must be C++
///
 
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>

#include "HepPID/ParticleIDTranslations.hh"
#include "HepPID/ParticleName.hh"
#include "HepPID/Version.hh"

extern "C" {
    void list_isajet_init_ ( );
    void flavor_( int *, int *, int *, int *, int *, int * );
    void get_label_( int * id, char *name );
}

int main()
{
    static char outfile[] = "examListIsajet.out";
    std::string title = "HepPID listing of Isajet translations";

    // initialize isajet 
    list_isajet_init_ ( );

    // open the output stream
    std::ofstream os( outfile );
    if( !os ) {
        std::cout << "error opening output file"  << std::endl;
        exit(1);
    }
    HepPID::writeVersion(os);

    os << "      " << title << std::endl;
    
    int i, j;
    int id, aid, fl1, fl2, fl3, js, indx;
    int pid;
    char cname[10];
    char acname[10];
    std::string hname;
    std::string pn;
    for( i=1; i<100005; ++i ) {
	// make sure names are empty
	for( j=0; j<10; ++j) { cname[j] = '\0'; }
	for( j=0; j<10; ++j) { acname[j] = '\0'; }
	// get info from isajet
	id = i;
	aid = 0;
	flavor_(&id, &fl1, &fl2, &fl3, &js, &indx );
	// we need both a valid index and a valid label
	// check the label only if there is a valid translation
	if ( indx > 0 ) {
	   get_label_(&id,cname);
	   aid = -id;
	   get_label_(&aid,acname);
	} else {
	   id = aid = 0;
	}
	
	// print particle
	if( id != 0 ) {
	    pid = HepPID::translateIsajettoPDT( id );
	    hname = std::string( cname );
	    if ( pid != 0 ) {
        	pn = HepPID::particleName( pid );
 		os << "Isajet: ";
		os.width(10);
		os << id << " " << hname;
		os << "  HepPID: " ;
		os.width(12);
		os << pid << " " << pn << std::endl;
	    } else if ( strncmp( cname, "ERR", 3 ) != 0 ) {
 		os << "Isajet: ";
		os.width(10);
		os << id << " with name \"" << hname;
		os << "\" has no HepPID translation " << std::endl;
	    }
	}
        // print antiparticle
	if( aid != 0 ) {
	    hname = std::string( acname );
	    pid = HepPID::translateIsajettoPDT( aid );
	    if ( pid != 0 ) {
        	pn = HepPID::particleName( pid );
 		os << "Isajet: ";
		os.width(10);
		os << aid << " " << hname;
		os << "  HepPID: " ;
		os.width(12);
		os << pid << " " << pn << std::endl;
	    } else if ( strncmp( acname, "ERR", 3 ) != 0 ) {
 		os << "Isajet: ";
		os.width(10);
		os << aid << " with name \"" << hname;
		os << "\" has no HepPID translation " << std::endl;
	    }
	}
    }
    
    return 0;
}
