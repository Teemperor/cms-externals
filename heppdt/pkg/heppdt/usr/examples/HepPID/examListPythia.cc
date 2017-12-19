// ----------------------------------------------------------------------
// examListPythia.cc
//
// ----------------------------------------------------------------------
///
/// \author Lynn Garren
///
/// list Pythia particle ID translations
/// Get ID list directly from Pythia
///
/// When mixing C++ and Fortran, the main program must be C++
///

#include <fstream>
#include <string>
#include <cstdlib>

#include "HepPID/ParticleIDTranslations.hh"
#include "HepPID/ParticleIDMethods.hh"
#include "HepPID/ParticleName.hh"
#include "HepPID/Version.hh"

extern "C" {
    void list_pythia_ ( );
    void getkf_( int *, int * );
    void getpyname_( int * , int *, char *name );
}

void writeLine( int & i, int & kf,  int & id, 
                std::string & name, std::string & pn, std::ofstream & os );

int main()
{
    int kf,akf,pid,apid,iok;
    char cname[17],caname[17];
    std::string name, aname, pn, apn;
    const char outfile[] = "examListPythia.out";
    std::string title = "HepPID listing of Pythia translations";
    // open the output file
    std::ofstream os( outfile );
    if( !os ) { 
      std::cerr << "cannot open " << outfile << std::endl;
      exit(-1);
    }
    HepPID::writeVersion(os);
    // get Pythia listing using the fortran pylist function
    // write the output of pylist to examListPythia.lpt
    list_pythia_();

    os << "      " << title << std::endl;
    
    for(int i=1; i<501; ++i) {
       getkf_( &i, &kf);
       if( kf != 0 ) {
           getpyname_(&kf, &iok, cname);
	   if( iok == 1 ) {
	       cname[16]='\0';
	       name = std::string( cname );
	       pid = HepPID::translatePythiatoPDT( kf );
               pn = HepPID::particleName( pid );
	       writeLine( i, kf, pid, name, pn, os );
	       if( HepPID::isValid( -pid ) ) {
		 akf=-kf;
        	 getpyname_(&akf, &iok, caname);
	         if( iok == 1 ) {
		     apid = HepPID::translatePythiatoPDT( akf );
        	     apn = HepPID::particleName( apid );
		     caname[16]='\0';
		     aname = std::string( caname );
		     writeLine( i, akf, apid, aname, apn, os );
	         }
	       }
	   }
       }
    }
    return 0;
}

void writeLine( int & i, int & kf, int & id, 
                std::string & name, std::string & pn, std::ofstream & os )
{
    os << "Pythia: ";
    os.width(7);
    os << i ;
    os.width(12);
    os << kf << " " << name;
    os << "  HepPID: " ;
    os.width(12);
    os << id << " " << pn << std::endl;
}
