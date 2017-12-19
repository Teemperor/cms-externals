// ----------------------------------------------------------------------
//
// ParticleID.hh
// Author:  Lynn Garren
//
//  various utilities to extract information from the particle ID
//
//  In the standard numbering scheme, the PID digits (base 10) are:
//            +/- n nr nl nq1 nq2 nq3 nj
//  It is expected that any 7 digit number used as a PID will adhere to 
//  the Monte Carlo numbering scheme documented by the PDG.
//  Note that many "new" particles not explicitly defined already 
//  can be expressed within this numbering scheme. 
//
// ----------------------------------------------------------------------
#ifndef PARTICLEID_HH
#define PARTICLEID_HH

#include <string>
#include <algorithm>	// swap()

// Particle names
#include "HepPID/ParticleName.hh"
// Translation free functions
#include "HepPID/ParticleIDTranslations.hh"

namespace HepPDT {

/// convert from 2J+1 to the actual spin value 
double spinitod( int js );
/// convert an actual spin to 2J+1
int spindtoi( double spin );

///  PID digits (base 10) are: n nr nl nq1 nq2 nq3 nj
///  The location enum provides a convenient index into the PID.
enum location { nj=1, nq3, nq2, nq1, nl, nr, n, n8, n9, n10 };

/// constituent quarks
struct Quarks {
    // -- constructors
    //
    Quarks( ) : nq1(0), nq2(0), nq3(0)  {}
    Quarks( short q1, short q2, short q3) : nq1(q1), nq2(q2), nq3(q3)  {}
    // -- data members
    //
    short nq1;
    short nq2;
    short nq3;
};

//! The ParticleID has various utilities to extract information from the particle ID

///
/// \class ParticleID
/// \author Lynn Garren
///
///  In the standard numbering scheme, the PID digits (base 10) are:
///            +/- n nr nl nq1 nq2 nq3 nj
///  It is expected that any 7 digit number used as a PID will adhere to 
///  the Monte Carlo numbering scheme documented by the PDG.
///  Note that particles not already explicitly defined 
///  can be expressed within this numbering scheme. 
///
class ParticleID {

public:

  // ---  birth/death:
  //
  /// create from an integer ID
  ParticleID( int pid = 0 );

  // ---  copying:
  //
  ParticleID( const ParticleID & orig );
  ParticleID& operator=( const ParticleID & );
  void swap( ParticleID & other ); 

  bool  operator <  ( ParticleID const & other ) const;
  bool  operator == ( ParticleID const & other ) const;

  // ---  accessors:
  //
  /// get the integer ID
  int    pid( )        const { return itsPID; }
  /// get the absolute value
  int abspid( )        const;

  // ---  boolean methods:
  //
  /// is this a valid ID?
  bool isValid( )   const;
  /// is this a valid meson ID?
  bool isMeson( )   const;
  /// is this a valid baryon ID?
  bool isBaryon( )  const;
  /// is this a valid diquark ID?
  bool isDiQuark( ) const;
  /// is this a valid hadron ID?
  bool isHadron( )  const;
  /// is this a valid lepton ID?
  bool isLepton( )  const;
  /// is this a valid ion ID?
  bool isNucleus( ) const;
  /// is this a valid pentaquark ID?
  bool isPentaquark( ) const;
  /// is this a valid SUSY ID?
  bool isSUSY( ) const;
  /// is this a valid R-hadron ID?
  bool isRhadron( ) const;
  /// is this a valid Dyon (magnetic monopole) ID?
  bool isDyon( ) const;
  /// Check for QBall or any exotic particle with electric charge beyond the qqq scheme
  /// Ad-hoc numbering for such particles is 100xxxx0, where xxxx is the charge in tenths. 
  bool isQBall( ) const;

  /// does this particle contain an up quark?
  bool hasUp( )      const;
  /// does this particle contain a down quark?
  bool hasDown( )    const;
  /// does this particle contain a strange quark?
  bool hasStrange( ) const;
  /// does this particle contain a charm quark?
  bool hasCharm( )   const;
  /// does this particle contain a bottom quark?
  bool hasBottom( )  const;
  /// does this particle contain a top quark?
  bool hasTop( )     const;

  /// jSpin returns 2J+1, where J is the total spin
  int  jSpin( )        const;
  /// sSpin returns 2S+1, where S is the spin
  int  sSpin( )        const;
  /// lSpin returns 2L+1, where L is the orbital angular momentum
  int  lSpin( )        const;
  /// return the first two digits if this is a "fundamental" particle
  int fundamentalID( ) const;
  /// returns everything beyond the 7th digit 
  /// (e.g. outside the standard numbering scheme)
  int extraBits( ) const;
  /// returns a list of 3 constituent quarks
  Quarks quarks( ) const;
  /// this is mostly for use by functions like addPDGParticles that have to 
  /// figure out the charge from the PID
  int threeCharge( ) const;
  /// get the actual charge, which might be fractional
  double charge() const;
  /// if this is a nucleus (ion), get A
  int A( ) const;
  /// if this is a nucleus (ion), get Z
  int Z( ) const;
  /// if this is a nucleus (ion), get nLambda
  int lambda( ) const;
  /// return the digit at a named location in the PID
  unsigned short digit(location) const;
  /// standard particle name
  const std::string PDTname() const { return HepPID::particleName( itsPID ); }

private:

  int itsPID;

};  // ParticleID

inline
void swap( ParticleID & first, ParticleID & second ) { first.swap( second ); }

} // HepPDT

//#include "HepPDT/ParticleID.icc"

#endif // PARTICLEID_HH
