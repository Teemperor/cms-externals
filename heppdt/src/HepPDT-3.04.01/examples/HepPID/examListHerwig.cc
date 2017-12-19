// -------------------------------------------------------------------
// examListHerwig.cc
//
// ----------------------------------------------------------------------
///
/// \author Lynn Garren
///
/// list Herwig particle ID translations
/// Get ID list directly from Herwig
///
/// When mixing C++ and Fortran, the main program must be C++
///
 
#include <fstream>
#include <iostream>
#include <cstdlib>
 
#include "HepPID/ParticleIDTranslations.hh"
#include "HepPID/ParticleName.hh"
#include "HepPID/Version.hh"

extern "C" {
// these functions are defined in examListHerwigInterface.F
    void list_herwig_init_ ( int * nevt );
    void list_herwig_end_ ( );
    void get_list_size_ ( int * );
    void get_herwig_name_( int * ihwg, int * id, char *name );
}

int main()
{
    int nevt=20;
    int i, j, iend, isize;
    int hid, id;
    char cname[10];
    std::string hname;
    std::string pn;
    static char outfile[] = "examListHerwig.out";
    std::string title = "HepPID listing of Herwig translations";

    // initialize herwig 
    list_herwig_init_ ( & nevt );

    // open the output stream
    std::ofstream os( outfile );
    if( !os ) {
        std::cout << "error opening output file"  << std::endl;
        exit(1);
    }
    HepPID::writeVersion(os);

    get_list_size_( & isize );
    os << "      " << title << std::endl;
    os << " number of Herwig particles: " << isize << std::endl;
    
    for( i=1, iend=isize+1; i<iend; ++i ) {
	// get info from herwig
	for( j=0; j<10; ++j) { cname[j] = '\0'; }
	get_herwig_name_( & i, & hid, cname );
	hname = std::string( cname );
	id = HepPID::translateHerwigtoPDT( hid );
        pn = HepPID::particleName( id );
 	os << "Herwig: ";
	os.width(7);
	os << i ;
	os.width(12);
	os << hid << " " << hname;
	os << "  HepPID: " ;
	os.width(12);
	os << id << " " << pn << std::endl;
    }

    list_herwig_end_();
    
    return 0;
}

