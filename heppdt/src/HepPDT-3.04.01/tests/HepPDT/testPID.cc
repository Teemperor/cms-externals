// ----------------------------------------------------------------------
// TestPID.cc
// Author: Lynn Garren
//
// test ParticleID
//
// ----------------------------------------------------------------------

#include <iostream>
#include <iomanip>

#include "HepPDT/defs.h"
#include "HepPDT/ParticleID.hh"

void testHadron( HepPDT::ParticleID&, int& );
void testLepton( HepPDT::ParticleID&, int& );
void testDiQuark( HepPDT::ParticleID&, int& );
void testNucleus( HepPDT::ParticleID&, int& );
void testUnknown( HepPDT::ParticleID&, int& );
void testHasMethods( HepPDT::ParticleID&, int& );
void testValid( HepPDT::ParticleID& );

int main()
{
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
        HepPDT::ParticleID pid( id[it] );
	nx = pid.digit(HepPDT::n);
	nr = pid.digit(HepPDT::nr);
	extra = pid.extraBits();
        std::cout << std::endl;
	std::cout << std::setw(18) << id[it] << ": " << nx << " " << nr 
	             << " " << pid.digit(HepPDT::nl)
                     << " " << pid.digit(HepPDT::nq1) << " " 
	             << pid.digit(HepPDT::nq2) << " " << pid.digit(HepPDT::nq3)
                     << " " <<  pid.digit(HepPDT::nj)
    	             << " extra bits " << extra << std::endl;
	js = pid.jSpin();
	ls = pid.lSpin();
	sid = pid.fundamentalID();
	chg = pid.threeCharge();
        if( !pid.isValid() ) {
	   std::cout << "**** Invalid PID: " << pid.pid() 
	             << " ****" << std::endl;
	} else  {
          testHadron( pid, id[it] );
          testLepton( pid, id[it] );
          testDiQuark( pid, id[it] );
          testNucleus( pid, id[it] );
          testUnknown( pid, id[it] );
        }
        testValid( pid );
	testHasMethods( pid, id[it] );
    }
}

void testHadron( HepPDT::ParticleID& pid, int& idn )
{
	int nx = pid.digit(HepPDT::n);
	int nr = pid.digit(HepPDT::nr);
	int extra = pid.extraBits();
	int js = pid.jSpin();
	int ls = pid.lSpin();
	int sid = pid.fundamentalID();
	int chg = pid.threeCharge();
        HepPDT::Quarks cqks = pid.quarks( );
	  if( pid.isHadron() ) {
	    if( pid.isMeson() ) {
	       std::cout << "meson   " << std::setw(10) << idn << ": " << nx
                         << " " << nr << " "  << ls << " " 
	                 << cqks.nq1 << " " << cqks.nq2 << " " << cqks.nq3
                         << " " << js << " " << sid << " " << chg
		         << " extra bits " << extra << std::endl;
            } else if( pid.isRhadron() ) {
	       std::cout << "Rhadron " << std::setw(10) << idn << ": " << nx
                         << " " << nr << " "  << ls << " " 
	                 << cqks.nq1 << " " << cqks.nq2 << " " << cqks.nq3 
                         << " " << js << " " << sid << " " << chg
		         << " extra bits " << extra << std::endl;
            } else if( pid.isBaryon() ) {
	       std::cout << "baryon  " << std::setw(10) << idn << ": " << nx
                         << " " << nr << " "  << ls << " " 
	                 << cqks.nq1 << " " << cqks.nq2 << " " << cqks.nq3 
                         << " " << js << " " << sid << " " << chg
		         << " extra bits " << extra << std::endl;
            } else {
	        std::cout << "**** undefined hadron: " << pid.pid() 
                          << " ****" << std::endl;
	    }
	  }
}

void testLepton( HepPDT::ParticleID& pid, int& idn )
{
        HepPDT::Quarks cqks = pid.quarks( );
	  if( pid.isLepton() ) {
	   std::cout << "lepton  " << std::setw(10) << idn 
	             << ": " << pid.digit(HepPDT::n)
                     << " "  << pid.digit(HepPDT::nr)
		     << " "  << pid.lSpin() << " " 
	             << cqks.nq1 << " " << cqks.nq2 << " " << cqks.nq3
                     << " " <<  pid.jSpin()
		     << " " <<  pid.fundamentalID()
		     << " " <<  pid.threeCharge()
		     << " extra bits " <<  pid.extraBits() << std::endl;
	  }
}

void testDiQuark( HepPDT::ParticleID& pid, int& idn )
{
        HepPDT::Quarks cqks = pid.quarks( );
	  if( pid.isDiQuark() ) {
	   std::cout << "diquark " << std::setw(10) << idn 
	             << ": " << pid.digit(HepPDT::n)
                     << " "  << pid.digit(HepPDT::nr)
		     << " "  << pid.lSpin() << " " 
	             << cqks.nq1 << " " << cqks.nq2 << " " << cqks.nq3 
                     << " " <<  pid.jSpin()
		     << " " <<  pid.fundamentalID()
		     << " " <<  pid.threeCharge()
		     << " extra bits " <<  pid.extraBits() << std::endl;
	  }
}

void testNucleus( HepPDT::ParticleID& pid, int& idn )
{
	  if( pid.isNucleus() ) {
	   std::cout << "ion    " << std::setw(11) << idn 
                     << ": " << pid.digit(HepPDT::n10)
                     << " " << std::setw(3) << pid.A() 
                     << " " << std::setw(3) << pid.Z() 
                     << " " << std::setw(3) << pid.lambda() 
                     << " " <<  pid.jSpin()
		     << " " <<  pid.fundamentalID()
		     << " " <<  pid.threeCharge()
		     << " extra bits " <<  pid.extraBits() << std::endl;
	  }
}

void testUnknown( HepPDT::ParticleID& pid, int& idn )
{
        HepPDT::Quarks cqks = pid.quarks( );
	  if( pid.isHadron() || pid.isLepton() || 
	      pid.isDiQuark() || pid.isNucleus() ) {
	  } else if( pid.isSUSY()  ) {
	   std::cout << "SUSY    " << std::setw(10) << idn 
	             << ": " << pid.digit(HepPDT::n)
                     << " "  << pid.digit(HepPDT::nr)
		     << " "  << pid.lSpin() << " " 
	             << cqks.nq1 << " " << cqks.nq2 << " " << cqks.nq3 
                     << " " <<  pid.jSpin()
		     << " " <<  pid.fundamentalID()
		     << " " <<  pid.threeCharge()
		     << " extra bits " <<  pid.extraBits() << std::endl;
	  } else if( pid.isDyon()  ) {
	   std::cout << "Dyon    " << std::setw(10) << idn 
	             << ": " << pid.digit(HepPDT::n)
                     << " "  << pid.digit(HepPDT::nr)
		     << " "  << pid.lSpin() << " " 
	             << cqks.nq1 << " " << cqks.nq2 << " " << cqks.nq3 
                     << " " <<  pid.jSpin()
		     << " " <<  pid.fundamentalID()
		     << " " <<  pid.threeCharge()
		     << " extra bits " <<  pid.extraBits() << std::endl;
	  } else if( pid.abspid() < 100 ) {
	   std::cout << "quark/lepton/boson " << std::setw(2) << idn 
	             << ": " << pid.digit(HepPDT::n)
                     << " "  << pid.digit(HepPDT::nr)
		     << " "  << pid.lSpin() << " " 
	             << cqks.nq1 << " " << cqks.nq2 << " " << cqks.nq3 
                     << " " <<  pid.jSpin()
		     << " " <<  pid.fundamentalID()
		     << " " <<  pid.threeCharge()
		     << " extra bits " <<  pid.extraBits() << std::endl;
	  } else {
	   std::cout << "unknown " << std::setw(10) << idn 
	             << ": " << pid.digit(HepPDT::n)
                     << " "  << pid.digit(HepPDT::nr)
		     << " "  << pid.lSpin() << " " 
	             << cqks.nq1 << " " << cqks.nq2 << " " << cqks.nq3 
                     << " " <<  pid.jSpin()
		     << " " <<  pid.fundamentalID()
		     << " " <<  pid.threeCharge()
		     << " extra bits " <<  pid.extraBits() << std::endl;
	  }
}

void testValid( HepPDT::ParticleID& pid )
{
        if( pid.isValid() ) {
	    int js = pid.jSpin();
	    int chg = pid.threeCharge();
	    std::cout << "total spin: " << js << " " << HepPDT::spinitod(js) ;
	    std::cout << " orbital angular momentum: " <<  pid.lSpin() ;
	    std::cout << " spin: " <<  pid.sSpin() ; 
	    std::cout << " charge: " << chg 
	              << " " << pid.charge() << std::endl;
            HepPDT::Quarks qlist = pid.quarks( );
	    std::cout <<  "quarks: " << qlist.nq1 
	              << " " << qlist.nq2 << " " << qlist.nq3 << std::endl;
        }
}

void testHasMethods( HepPDT::ParticleID& pid, int& idn )
{
  if( pid.hasUp() ) {
      std::cout << "Particle " << pid.PDTname()  << " " << pid.pid() 
                << " has an up quark" << std::endl;
  }
  if( pid.hasDown() ) {
      std::cout << "Particle " << pid.PDTname()  << " " << pid.pid() 
                << " has a down quark" << std::endl;
  }
  if( pid.hasStrange() ) {
      std::cout << "Particle " << pid.PDTname()  << " " << pid.pid() 
                << " has a strange quark" << std::endl;
  }
  if( pid.hasCharm() ) {
      std::cout << "Particle " << pid.PDTname() << " " << pid.pid() 
                << " has a charmed quark" << std::endl;
  }
  if( pid.hasBottom() ) {
      std::cout << "Particle " << pid.PDTname() << " " << pid.pid() 
                << " has a bottom quark" << std::endl;
  }
  if( pid.hasTop() ) {
      std::cout << "Particle " << pid.PDTname() << " " << pid.pid() 
                << " has a top quark" << std::endl;
  }
} 

