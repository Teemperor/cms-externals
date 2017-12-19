// ----------------------------------------------------------------------
//
// TempParticleData.hh
// Author:  Lynn Garren
//
// temporary structure to hold particles while getting all the info
// ----------------------------------------------------------------------
#ifndef TEMPPARTICLEDATA_HH
#define TEMPPARTICLEDATA_HH

#include <string>
#include <vector>

#include "HepPDT/SpinState.hh"
#include "HepPDT/ParticleID.hh"
#include "HepPDT/Measurement.hh"

namespace HepPDT {

/// free function
/// Given the lifetime, calculate the width.
double calculateWidthFromLifetime( double );

/// Temporary constituent (e.g., quark) information
struct TempConstituent {
    TempConstituent( ParticleID p = ParticleID(0), int m = -1 )
      : tempConstituentPID(p), tempMultiplicity(m)  {}
    ParticleID tempConstituentPID;
    int                tempMultiplicity;
};

/// temporary holder for decay data
struct TempDecayData  {
  /// Type of decay.
  /// Use string to hold int if necessary.
  std::string              tempDecayName;
  double                   tempBranchingFraction;
  /// list of decay particles
  std::vector<std::string> tempDaughterList;
  /// other decay parameters
  std::vector<double>      tempDecayParameters;
};

/// useful typedef
  typedef std::vector<TempDecayData>  TDDlist;

/// Hold Alias information from EvtGen
struct TempAliasData {

  // ---  birth/death:
  //
  TempAliasData();	//!< used in a map<>

  std::string    tempAlias;		//!< the alias
  std::string    tempAliasedParticle;	//!< the "real" particle
  std::string    tempChargeConj;	//!< set if there is a charge conjugate alias
  TDDlist        tempAliasDecayList;	//!< decay list for the alias
};

/// temporary holder for Particle Data information
struct TempParticleData  {

  // ---  birth/death:
  //
  TempParticleData();	//!< used in a map<>
  /// construct a basic TempParticleData from just the particle ID
  explicit TempParticleData( int id );
  /// construct a basic TempParticleData from just the ParticleID
  explicit TempParticleData( ParticleID pid );
  /// given all the information, construct a TempParticleData
  TempParticleData( int id, std::string const & name, 
                    std::string const & source, int oid, double charge,
                    SpinState const & Spin, 
		    Measurement const & mass, Measurement const & wid );

  // ---  copying:
  //
  TempParticleData( TempParticleData const & orig );
  TempParticleData & operator = ( TempParticleData const & rhs );
  void swap( TempParticleData & other );
  /// given a particle definition, create an antiparticle
  TempParticleData  antiparticle( std::string const & name );
 
  /// -- mutator - get spin state and constituent list from PID
  bool processPID();

  // -- data members (public by design)
  //
  ParticleID                   tempID;
  std::string                  tempParticleName;
  std::string                  tempSource;
  int                          tempOriginalID;
  double                       tempCharge;
  double                       tempColorCharge;
  SpinState                    tempSpin;
  Measurement                  tempMass;
  Measurement                  tempWidth;
  double                       tempLowCutoff;
  double                       tempHighCutoff;
  std::vector<TempConstituent> tempQuarks;
  TDDlist                      tempDecayList; 

};  // TempParticleData

inline
void swap( TempParticleData & first, TempParticleData & second ) { first.swap( second ); }

}	// HepPDT

#endif // TEMPPARTICLEDATA_HH
