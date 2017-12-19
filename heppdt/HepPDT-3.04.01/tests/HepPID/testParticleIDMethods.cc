// ----------------------------------------------------------------------
// testParticleIDMethods.cc
// Author: Lynn Garren
//
// test ParticleID
//
// ----------------------------------------------------------------------

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "HepPID/ParticleIDMethods.hh"
#include "HepPID/ParticleName.hh"

using std::setw;

int main()
{
    // open an output file
    //
    const char outfile[] = "testParticleIDMethods.out";
    std::ofstream os( outfile );
    if( !os ) {
      std::cerr << "cannot open " << outfile << std::endl;
      exit(-1);
    }
    // define numbers to check
    int id[40] = { 5, 25, 15, 213, -3214, 10213, 9050225, -200543, 129050225,
                   2000025, 3101, 3301, -2212, 1000020040, -1000060120, 555,
		   5000040, 5100005, 24, 5100024, 5100025, 9221132, 
                   4111370, -4120240, 4110050, 10013730,
		   1000993, 1000612, 1000622, 1000632, 1006213, 1000652, 
		   1009113, 1009213, 1009323,
		   1093114, 1009333, 1006313, 1092214, 1006223 };
    int it;
    int nr, nx;
    int chg, sid, extra;
    int js, ls;
    for( it=0; it < 40; it++ ) {
        int pid = id[it];
	nx = HepPID::digit(HepPID::n,pid);
	nr = HepPID::digit(HepPID::nr,pid);
	extra = HepPID::extraBits(pid);
        os << std::endl;
	os << setw(21) << pid << ": " << nx << " " << nr 
	             << " " << HepPID::digit(HepPID::nl,pid)
                     << " " << HepPID::digit(HepPID::nq1,pid) 
		     << " " << HepPID::digit(HepPID::nq2,pid) 
		     << " " << HepPID::digit(HepPID::nq3,pid)
                     << " " << HepPID::digit(HepPID::nj,pid)
    	             << " extra bits " << extra << std::endl;
	js = HepPID::jSpin(pid);
	ls = HepPID::lSpin(pid);
	sid = HepPID::fundamentalID(pid);
	chg = HepPID::threeCharge(pid);
        if( !HepPID::isValid(pid) ) {
	   os << "**** Invalid PID: " << pid << " ****" << std::endl;
	} else  {
	  if( HepPID::isHadron(pid) ) {
	    if( HepPID::isMeson(pid) ) {
	       os << "meson      " << setw(10) << pid << ": " << nx
                         << " " << nr << " "  << ls << " " 
                         << js << " " << sid << " " << chg
		         << " extra bits " << extra << std::endl;
            } else if( HepPID::isBaryon(pid) ) {
	       os << "baryon     " << setw(10) << pid << ": " << nx
                         << " " << nr << " "  << ls << " " 
                         << js << " " << sid << " " << chg
		         << " extra bits " << extra << std::endl;
	    } else if( HepPID::isRhadron(pid) ) {
	     os << "Rhadron    " << setw(10) << pid << ": " << nx 
                       << " " << nr << " "  << ls << " " 
                       << js << " " << sid << " " << chg
		       << " extra bits " << extra << std::endl;
	    } else if( HepPID::isPentaquark(pid) ) {
	     os << "pentaquark " << setw(10) << pid << ": " << nx 
                       << " " << nr << " "  << ls << " " 
                       << js << " " << sid << " " << chg
		       << " extra bits " << extra << std::endl;
            } else {
	        os << "**** undefined hadron: " << pid
                          << " ****" << std::endl;
	    }
	  } else if( HepPID::isLepton(pid) ) {
	   os << "lepton     " << setw(10) << pid << ": " << nx 
                     << " " << nr << " "  << ls << " " 
                     << js << " " << sid << " " << chg
		     << " extra bits " << extra << std::endl;
	  } else if( HepPID::isDiQuark(pid) ) {
	   os << "diquark    " << setw(10) << pid << ": " << nx 
                     << " " << nr << " "  << ls << " " 
                     << js << " " << sid << " " << chg
		     << " extra bits " << extra << std::endl;
	  } else if( HepPID::isSUSY(pid) ) {
	   os << "SUSY       " << setw(10) << pid << ": " << nx 
                     << " " << nr << " "  << ls << " " 
                     << js << " " << sid << " " << chg
		     << " extra bits " << extra << std::endl;
	  } else if( HepPID::isDyon(pid) ) {
	   os << "Dyon       " << setw(10) << pid << ": " << nx 
                     << " " << nr << " "  << ls << " " 
                     << js << " " << sid << " " << chg
		     << " extra bits " << extra << std::endl;
	  } else if( HepPID::isQBall(pid) ) {
	   os << "QBall      " << setw(10) << pid << ": " << nx 
                     << " " << nr << " "  << ls << " " 
                     << js << " " << sid << " " << chg
		     << " extra bits " << extra << std::endl;
	  } else if( HepPID::isNucleus(pid) ) {
	   os << "ion       " << setw(11) << pid 
                     << ": " << HepPID::digit(HepPID::n10,pid)
                     << " " << setw(3) << HepPID::A(pid) 
                     << " " << setw(3) << HepPID::Z(pid) 
                     << " " << setw(3) << HepPID::lambda(pid) 
                     << " " << js << " " << sid << " " << chg
		     << " extra bits " << extra << std::endl;
	  } else if( sid > 0 && HepPID::abspid(pid) < 100 ) {
	   os << "quark/lepton/boson " << setw(2) << pid << ": " << nx 
                     << " " << nr << " "  << ls << " " 
                     << js << " " << sid << " " << chg
		     << " extra bits " << extra << std::endl;
	  } else {
	   os << "unknown    " << setw(10) << pid << ": " << nx 
                     << " " << nr << " "  << ls << " " 
                     << js << " " << sid << " " << chg
		     << " extra bits " << extra << std::endl;
	  }
	  // special case where a baryon is also a nucleus
	  if( HepPID::isNucleus(pid) && HepPID::isBaryon(pid)) {
	   os << "ion       " << setw(11) << pid 
                     << ": " << HepPID::digit(HepPID::n10,pid)
                     << " " << setw(3) << HepPID::A(pid) 
                     << " " << setw(3) << HepPID::Z(pid) 
                     << " " << setw(3) << HepPID::lambda(pid) 
                     << " " << js << " " << sid << " " << chg
		     << " extra bits " << extra << std::endl;
	  }
        }
        if( HepPID::isValid(pid) ) {
	    os << "total spin: " << HepPID::jSpin(pid) ;
	    os << " orbital angular momentum: " <<  HepPID::lSpin(pid) ;
	    os << " spin: " <<  HepPID::sSpin(pid) ; 
	    os << " charge: " << HepPID::threeCharge(pid) 
	       << " " << HepPID::charge(pid) << std::endl;
	    os << "Particle " <<  HepPID::particleName(pid) << std::endl;
	    // quark content
	    if( HepPID::hasUp(pid) ) {
		os << "Particle " <<  HepPID::particleName(pid) << " " << pid
		   << " has an up quark"  << std::endl;
	    }
	    if( HepPID::hasDown(pid) ) {
		os << "Particle " <<  HepPID::particleName(pid) << " " << pid
		   << " has a down quark"  << std::endl;
	    }
	    if( HepPID::hasStrange(pid) ) {
		os << "Particle " <<  HepPID::particleName(pid) << " " << pid
		   << " has a strange quark"  << std::endl;
	    }
	    if( HepPID::hasCharm(pid) ) {
		os << "Particle " <<  HepPID::particleName(pid) << " " << pid
		   << " has a charm quark"  << std::endl;
	    }
	    if( HepPID::hasBottom(pid) ) {
		os << "Particle " <<  HepPID::particleName(pid) << " " << pid
		   << " has a bottom quark"  << std::endl;
	    }
	    if( HepPID::hasTop(pid) ) {
		os << "Particle " <<  HepPID::particleName(pid) << " " << pid
		   << " has a top quark"  << std::endl;
	    }
        }
	// checking abspid
	int apid = HepPID::abspid(pid);
	int abpid = abs(pid);
	if( apid != abpid ) {
	    os << "Error: Abspid returns " << apid 
        	 << " instead of " << abpid 
        	 << " for " << pid << std::endl;
	}
    }
    return 0;
}
