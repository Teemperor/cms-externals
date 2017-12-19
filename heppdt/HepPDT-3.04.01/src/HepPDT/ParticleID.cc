// ----------------------------------------------------------------------
//
// ParticleID.cc
//
// ----------------------------------------------------------------------

#include <stdlib.h>
#include <cmath>	// for pow()

#include "HepPDT/defs.h"
#include "HepPDT/ParticleID.hh"
#include "HepPID/ParticleIDMethods.hh"

namespace HepPDT {


ParticleID::ParticleID( int pid )
: itsPID            ( pid )
{ ; }

ParticleID::ParticleID( const ParticleID & orig )
: itsPID            ( orig.itsPID )
{ ; }

ParticleID& ParticleID::operator=( const ParticleID & rhs )
{ 
  ParticleID temp( rhs );
  swap( temp );
  return *this;
}

void ParticleID::swap( ParticleID & other ) 
{
  std::swap( itsPID, other.itsPID );
}

bool  ParticleID::operator <  ( ParticleID const & other ) const
{
  return itsPID < other.itsPID;
}

bool  ParticleID::operator == ( ParticleID const & other ) const
{
  return itsPID == other.itsPID;
}

int ParticleID::abspid( ) const
{
  /// return a value greater than or equal to zero
  return (itsPID < 0) ? -itsPID : itsPID;
}

//  split the PID into constituent integers
unsigned short ParticleID::digit( location loc ) const
{
    //  PID digits (base 10) are: n nr nl nq1 nq2 nq3 nj
    //  the location enum provides a convenient index into the PID
    int numerator = (int) std::pow(10.0,(loc-1));
    return (abspid()/numerator)%10;
}

// returns everything beyond the 7th digit (e.g. outside the numbering scheme)
int ParticleID::extraBits( ) const
{
    return abspid()/10000000;
}

//  return the first two digits if this is a "fundamental" particle
//  ID = 100 is a special case (internal generator ID's are 81-100)
//  also, 101 and 102 are now used (by HepPID) for geantinos
int ParticleID::fundamentalID( ) const
{
    if( ( digit(n10) == 1 ) && ( digit(n9) == 0 ) ) { return 0; }
    if( digit(nq2) == 0 && digit(nq1) == 0) {
        return abspid()%10000;
    } else if( abspid() <= 102 ) {
        return abspid();
    } else {
        return 0;
    }
}

bool ParticleID::hasUp( ) const
{
    return HepPID::hasUp(itsPID);
}

bool ParticleID::hasDown( ) const
{
    return HepPID::hasDown(itsPID);
}

bool ParticleID::hasStrange( ) const
{
    return HepPID::hasStrange(itsPID);
}

bool ParticleID::hasCharm( ) const
{
    return HepPID::hasCharm(itsPID);
}

bool ParticleID::hasBottom( ) const
{
    return HepPID::hasBottom(itsPID);
}

bool ParticleID::hasTop( ) const
{
    return HepPID::hasTop(itsPID);
}

//  check to see if this is a valid PID
bool ParticleID::isValid( ) const
{
    if( extraBits() > 0 ) {
        if( isNucleus() )   { return true; }
        if( isQBall() )   { return true; }
        return false; 
    }
    // SUSY signature
    if( isSUSY() ) { return true; }
    // R-hadron signature
    if( isRhadron() ) { return true; }
    // Dyon (magnetic monopole) signature
    if( isDyon() ) { return true; }
    // Meson signature
    if( isMeson() )   { return true; }
    // Baryon signature
    if( isBaryon() )  { return true; }
    // DiQuark signature
    if( isDiQuark() ) { return true; }
    // fundamental particle
    if( fundamentalID() > 0 ) { return true; }
    // pentaquark
    if( isPentaquark() ) { return true; }
    // don't recognize this number
    return false;
}

//  check to see if this is a valid pentaquark
bool ParticleID::isPentaquark( ) const
{
    // a pentaquark is of the form 9abcdej,
    // where j is the spin and a, b, c, d, and e are quarks
    if( extraBits() > 0 ) { return false; }
    if( digit(n) != 9 )  { return false; }
    if( digit(nr) == 9 || digit(nr) == 0 )  { return false; }
    if( digit(nj) == 9 || digit(nl) == 0 )  { return false; }
    if( digit(nq1) == 0 )  { return false; }
    if( digit(nq2) == 0 )  { return false; }
    if( digit(nq3) == 0 )  { return false; }
    if( digit(nj) == 0 )  { return false; }
    // check ordering
    if( digit(nq2) > digit(nq1) )  { return false; }
    if( digit(nq1) > digit(nl) )  { return false; }
    if( digit(nl) > digit(nr) )  { return false; }
    return true;
}

// is this an R-hadron?
bool ParticleID::isRhadron( ) const
{
    // an R-hadron is of the form 10abcdj,
    // where j is the spin and a, b, c, and d are quarks or gluons
    if( extraBits() > 0 ) { return false; }
    if( digit(n) != 1 )  { return false; }
    if( digit(nr) != 0 )  { return false; }
    // make sure this isn't a SUSY particle
    if( isSUSY() ) { return false; }
    // All R-hadrons have at least 3 core digits
    if( digit(nq2) == 0 )  { return false; }
    if( digit(nq3) == 0 )  { return false; }
    if( digit(nj) == 0 )  { return false; }
    return true;
}

// is this a valid Dyon (magnetic monopole) ID?
bool ParticleID::isDyon( ) const
{
    ///Magnetic monopoles and Dyons are assumed to have one unit of 
    ///Dirac monopole charge and a variable integer number xyz units 
    ///of electric charge.  
    ///
    ///Codes 411xyz0 are then used when the magnetic and electrical 
    ///charge sign agree and 412xyz0 when they disagree, 
    ///with the overall sign of the particle set by the magnetic charge.  
    ///For now no spin information is provided.
    ///
    if( extraBits() > 0 ) { return false; }
    if( digit(n) != 4 )  { return false; }
    if( digit(nr) != 1 )  { return false; }
    if( (digit(nl) != 1) && (digit(nl) != 2) )  { return false; }
    // All Dyons have at least 1 core digit
    if( digit(nq3) == 0 )  { return false; }
    // Dyons have spin zero for now
    if( digit(nj) != 0 )  { return false; }
    return true;
}

// Check for QBall or any exotic particle with electric charge beyond the qqq scheme
bool ParticleID::isQBall( ) const
{
    // Ad-hoc numbering for such particles is 100xxxx0, 
    // where xxxx is the charge in tenths. 
    if( extraBits() != 1 ) { return false; }
    if( digit(n) != 0 )  { return false; }
    if( digit(nr) != 0 )  { return false; }
    // check the core number
    if( (abspid()/10)%10000 == 0 )  { return false; }
    // these particles have spin zero for now
    if( digit(nj) != 0 )  { return false; }
    return true;
}

// is this a SUSY?
bool ParticleID::isSUSY( ) const
{
    // fundamental SUSY particles have n = 1 or 2
    if( extraBits() > 0 ) { return false; }
    if( digit(n) != 1 && digit(n) != 2 )  { return false; }
    if( digit(nr) != 0 )  { return false; }
    // check fundamental part
    if( fundamentalID() == 0 )  { return false; }
    return true;
}

// This implements the 2006 Monte Carlo nuclear code scheme.
// Ion numbers are +/- 10LZZZAAAI. 
// AAA is A - total baryon number
// ZZZ is Z - total charge
// L is the total number of strange quarks.
// I is the isomer number, with I=0 corresponding to the ground state.
bool ParticleID::isNucleus( ) const
{
     // a proton can also be a Hydrogen nucleus
     if( abspid() == 2212 ) { return true; }
     // new standard: +/- 10LZZZAAAI
     if( ( digit(n10) == 1 ) && ( digit(n9) == 0 ) ) {
        // charge should always be less than or equal to baryon number
        if( A() >= Z() ) { return true; }
     }
     return false;
}

// return A if this is a nucleus
int ParticleID::A( ) const
{
    // a proton can also be a Hydrogen nucleus
    if( abspid() == 2212 ) { return 1; }
    if( ( digit(n10) != 1 ) || ( digit(n9) != 0 ) ) { return 0; }
    return (abspid()/10)%1000;
}

// return Z if this is a nucleus
int ParticleID::Z( ) const
{
    // a proton can also be a Hydrogen nucleus
    if( abspid() == 2212 ) { return 1; }
    if( ( digit(n10) != 1 ) || ( digit(n9) != 0 ) ) { return 0; }
    return (abspid()/10000)%1000;
}

// return nLambda if this is a nucleus
int ParticleID::lambda( ) const
{
    if( ! isNucleus() ) { return 0; }
    // a proton can also be a Hydrogen nucleus
    if( abspid() == 2212 ) { return 0; }
    return digit(n8);
}

//  check to see if this is a valid meson
bool ParticleID::isMeson( ) const
{
    if( extraBits() > 0 ) { return false; }
    if( abspid() <= 100 ) { return false; }
    if( fundamentalID() <= 100 && fundamentalID() > 0 ) { return false; }
    int aid = abspid();
    int id = pid();
    if( aid == 130 || aid == 310 || aid == 210 ) { return true; }
    // EvtGen uses some odd numbers
    if( aid == 150 || aid == 350 || aid == 510 || aid == 530 ) { return true; }
    // pomeron, etc.
    if( id == 110 || id == 990 || id == 9990 ) { return true; }
    if( digit(nj) > 0 && digit(nq3) > 0 && digit(nq2) > 0 && digit(nq1) == 0 ) {
        // check for illegal antiparticles
        if( digit(nq3) == digit(nq2) && pid() < 0 ) {
            return false;
        } else {
            return true;
        }
    }
    return false;
}

//  check to see if this is a valid meson
bool ParticleID::isLepton( ) const
{
    if( extraBits() > 0 ) { return false; }
    if( fundamentalID() >= 11 && fundamentalID() <= 18 ) { return true; }
    return false;
}

//  check to see if this is a valid hadron
bool ParticleID::isHadron( ) const
{
    if( extraBits() > 0 ) { return false; }
    if( isMeson() )   { return true; }
    if( isBaryon() )  { return true; }
    if( isPentaquark() ) { return true; }
    return false;
}

//  check to see if this is a valid diquark
bool ParticleID::isDiQuark( ) const
{
    if( extraBits() > 0 ) { return false; }
    if( abspid() <= 100 ) { return false; }
    if( fundamentalID() <= 100 && fundamentalID() > 0 ) { return false; }
    if( digit(nj) > 0 && digit(nq3) == 0 && digit(nq2) > 0 && digit(nq1) > 0 ) {  // diquark signature
       // EvtGen uses the diquarks for quark pairs, so, for instance, 
       //   5501 is a valid "diquark" for EvtGen
       //if( digit(nj) == 1 && digit(nq2) == digit(nq1) ) { 	// illegal
       //   return false; 
       //} else {
          return true;
       //}
    }
    return false;
}

//  check to see if this is a valid baryon
bool ParticleID::isBaryon( ) const
{
    if( extraBits() > 0 ) { return false; }
    if( abspid() <= 100 ) { return false; }
    if( fundamentalID() <= 100 && fundamentalID() > 0 ) { return false; }
    if( abspid() == 2110 || abspid() == 2210 ) { return true; }
    if( digit(nj) > 0 && digit(nq3) > 0 && digit(nq2) > 0 && digit(nq1) > 0 ) { return true; }
    return false;
}

//  return the total spin as 2J+1
int ParticleID::jSpin( ) const
{
    if( fundamentalID() > 0 && fundamentalID() <= 100 ) { 
	// some of these are known
	int fund = fundamentalID();
	if( fund > 0 && fund < 7 ) return 2;
	if( fund == 9 ) return 3; 
	if( fund > 10 && fund < 17 ) return 2;
	if( fund > 20 && fund < 25 ) return 3;
        return 0; 
    } else if( extraBits() > 0 ) { 
        return 0; 
    }
    return abspid()%10;
}

//  return the orbital angular momentum - valid for mesons only
int ParticleID::lSpin( ) const
{
    int nl, tent, js;
    if( !isMeson() ) { return 0; }
    nl = (abspid()/10000)%10;
    tent = (abspid()/1000000)%10;
    js = abspid()%10;
    if( tent == 9 ) { return 0; }	// tentative assignment
    if( nl == 0 && js == 3 ) { 
        return 0;
    } else if( nl == 0 && js == 5 ) {
        return 1;
    } else if( nl == 0 && js == 7 ) {
        return 2;
    } else if( nl == 0 && js == 9 ) {
        return 3;
    } else if( nl == 0  && js == 1 ) {
        return 0;
    } else if( nl == 1  && js == 3 ) {
        return 1;
    } else if( nl == 1  && js == 5 ) {
        return 2;
    } else if( nl == 1  && js == 7 ) {
        return 3;
    } else if( nl == 1  && js == 9 ) {
        return 4;
    } else if( nl == 2  && js == 3 ) {
        return 1;
    } else if( nl == 2  && js == 5 ) {
        return 2;
    } else if( nl == 2  && js == 7 ) {
        return 3;
    } else if( nl == 2  && js == 9 ) {
        return 4;
    } else if( nl == 1  && js == 1 ) {
        return 1;
    } else if( nl == 3  && js == 3 ) {
        return 2;
    } else if( nl == 3  && js == 5 ) {
        return 3;
    } else if( nl == 3  && js == 7 ) {
        return 4;
    } else if( nl == 3  && js == 9 ) {
        return 5;
    }
    // default to zero
    return 0;
}

//  return the spin - valid for mesons only
int ParticleID::sSpin( ) const
{
    int nl, tent, js;
    if( !isMeson() ) { return 0; }
    nl = (abspid()/10000)%10;
    tent = (abspid()/1000000)%10;
    js = abspid()%10;
    if( tent == 9 ) { return 0; }	// tentative assignment
    if( nl == 0 && js >= 3 ) { 
        return 1;
    } else if( nl == 0  && js == 1 ) {
        return 0;
    } else if( nl == 1  && js >= 3 ) {
        return 0;
    } else if( nl == 2  && js >= 3 ) {
        return 1;
    } else if( nl == 1  && js == 1 ) {
        return 1;
    } else if( nl == 3  && js >= 3 ) {
        return 1;
    }
    // default to zero
    return 0;
}

int ParticleID::threeCharge( ) const
{
    int charge=0;
    int ida, sid;
    unsigned short q1, q2, q3;
    static int ch100[100] = { -1, 2,-1, 2,-1, 2,-1, 2, 0, 0,
                       -3, 0,-3, 0,-3, 0,-3, 0, 0, 0,
                        0, 0, 0, 3, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 3, 0, 0, 3, 0, 0, 0,
                        0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 6, 3, 6, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    q1 = digit(nq1);
    q2 = digit(nq2);
    q3 = digit(nq3);
    ida = abspid();
    sid = fundamentalID();
    if( ida == 0 ) {      // illegal
        return 0;
    } else if( extraBits() > 0 ) {
        if( isNucleus() ) {      // ion
            return 3*Z();
	} else if( isQBall() ) {		// QBall
	   charge = 3*((ida/10)%10000);
	} else { // not an ion
            return 0;
	}
    } else if( isDyon() ) {		// Dyon
        charge = 3*( (ida/10)%1000 );
	// this is half right
	// the charge sign will be changed below if pid < 0
	if( digit(nl) == 2 ) {
            charge = -charge; 
	}
    } else if( sid > 0 && sid <= 100 ) {	// use table
        charge = ch100[sid-1];
        if(ida==1000017 || ida==1000018) { charge = 0; }
        if(ida==1000034 || ida==1000052) { charge = 0; }
        if(ida==1000053 || ida==1000054) { charge = 0; }
        if(ida==5100061 || ida==5100062) { charge = 6; }
    } else if( digit(nj) == 0 ) { 		// KL, Ks, or undefined
        return 0;
    } else if( ( q1 == 0 ) || (isRhadron() && (q1 == 9) ) ) {	// meson			// mesons
        if( q2 == 3 || q2 == 5 ) {
            charge = ch100[q3-1] - ch100[q2-1];
        } else {
            charge = ch100[q2-1] - ch100[q3-1];
        }
    } else if( q3 == 0 ) {			// diquarks
        charge = ch100[q2-1] + ch100[q1-1];
    } else if( isBaryon() || (isRhadron() && (digit(nl) == 9) ) ) {  // baryon 			// baryons
        charge = ch100[q3-1] + ch100[q2-1] + ch100[q1-1];
    }
    if( charge == 0 ) {
        return 0;
    } else if( pid() < 0 ) {
        charge = -charge; 
    }
    return charge;
}

double ParticleID::charge( ) const
{
   int tc = threeCharge();
   if( isQBall() ) {
       return double(tc)/30.;
   } else {
       return double(tc)/3.;
   }
}

} // HepPDT
