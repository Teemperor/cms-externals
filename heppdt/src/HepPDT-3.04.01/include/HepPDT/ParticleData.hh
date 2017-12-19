// ----------------------------------------------------------------------
//
// ParticleData.hh
// Author: Lynn Garren
//
// This is a base class for particle data information
// ----------------------------------------------------------------------
#ifndef PARTICLEDATA_HH
#define PARTICLEDATA_HH

#include <string>
#include <vector>

#include "HepPDT/ParticleID.hh"
#include "HepPDT/SpinState.hh"
#include "HepPDT/Constituent.hh"
#include "HepPDT/ResonanceStructure.hh"
#include "HepPDT/TempParticleData.hh"

namespace HepPDT {

//! The ParticleData class holds data for a single particle in the table.

///
/// \class ParticleData
/// \author Lynn Garren
///
/// This class holds the basic paritcle data - name, ID, spin, etc.
/// Methods are provided to get various information about the particle.
/// All the particle information is accessed via this class.
///
class ParticleData {

public:

  // ---  birth/death:
  //
  /// create ParticleData from the temporary information
  ParticleData( const TempParticleData & tpd );
  ~ParticleData();

  // ---  copying:
  //
  void  swap ( ParticleData & rhs );
  ParticleData( const ParticleData & orig );
  ParticleData & operator = ( const ParticleData & rhs );

  // ---  accessors:
  //
  /// return particle name as defined by user input
  const std::string &         name()        const { return itsParticleName; }
  /// return PDG particle name
  const std::string           PDTname()     const { return itsID.PDTname(); }
  /// the name of the input source
  const std::string &         source()      const { return itsSource; }
  /// return the ParticleID
  ParticleID                  ID()          const { return itsID; }
  /// return the integer ID
  int                         pid( )        const { return itsID.pid(); }
  /// return untranslated integer ID
  int                         originalID( ) const { return itsOriginalID; }
  /// return charge
  double                      charge()      const { return itsCharge; }
  /// color information
  double                      color()       const { return itsColorCharge; }
  /// spin information
  SpinState                   spin()        const { return itsSpin; }
  /// mass
  Measurement                 mass()        const { return itsResonance.mass(); }
  /// return the total width
  Measurement                 totalWidth()  const { return itsResonance.totalWidth(); }
  /// lower cutoff of allowed width values
  double                      lowerCutoff() const { return itsResonance.lowerCutoff(); }
  /// upper cutoff of allowed width values
  double                      upperCutoff() const { return itsResonance.upperCutoff(); }
  /// calculate the lifetime
  Measurement          lifetime()    const { return itsResonance.lifetime(); }
  /// number of constituent particles (e.g., quarks)
  int                         numConstituents() const { return itsQuarks.size(); }
  /// constituent information
  Constituent          constituent( unsigned int i ) const;
  /// ParticleID for a constituent particle
  ParticleID           constituentParticle( unsigned int i ) const;
  /// resonance (width) information
  ResonanceStructure const    resonance()   const { return itsResonance; }
  
  /// output information about this particle
  void write( std::ostream & os ) const;
  /// output all information about a particle EXCEPT its decays
  /// This method is redundant with write(), but retained for backwards compatibility.
  void writeParticleInfo( std::ostream & os ) const;
  /// output the translation information for this particle
  void writeParticleTranslation( std::ostream & os ) const;
 
  // --- booleans:
  //
  /// is this a valid meson?
  bool isMeson( )   const { return itsID.isMeson(); }
  /// is this a valid baryon?
  bool isBaryon( )  const { return itsID.isBaryon(); }
  /// is this a valid diquark?
  bool isDiQuark( ) const { return itsID.isDiQuark(); }
  /// is this a valid hadron?
  bool isHadron( )  const { return itsID.isHadron(); }
  /// is this a valid lepton?
  bool isLepton( )  const { return itsID.isLepton(); }
  /// is this a valid ion?
  bool isNucleus( ) const { return itsID.isNucleus(); }
  /// is this a valid pentaquark ID?
  bool isPentaquark( ) const { return itsID.isPentaquark(); }
  /// is this a valid SUSY ID?
  bool isSUSY( ) const { return itsID.isSUSY(); }
  /// is this a valid R-hadron ID?
  bool isRhadron( ) const { return itsID.isRhadron(); }
  /// is this a valid Dyon (magnetic monopole) ID?
  bool isDyon( ) const { return itsID.isDyon(); }
  /// Check for QBall or any exotic particle with electric charge beyond the qqq scheme
  /// Ad-hoc numbering for such particles is 100xxxx0, where xxxx is the charge in tenths. 
  bool isQBall( ) const { return itsID.isQBall(); }
  // has methods look at Constituents, not PID
  /// does this particle contain an up quark?
  bool hasUp()      const;
  /// does this particle contain a down quark?
  bool hasDown()    const;
  /// does this particle contain a strange quark?
  bool hasStrange() const;
  /// does this particle contain a charm quark?
  bool hasCharm()   const;
  /// does this particle contain a bottom quark?
  bool hasBottom()  const;
  /// does this particle contain a top quark?
  bool hasTop()     const;
  /// is this particle allowed to decay?
  bool isStable() const;

  /// compare masses
  bool operator <  ( const ParticleData & other ) const;
  /// use PID and ignore everything else
  bool operator == ( const ParticleData & other ) const;

  // ---  mutators:
  //
  /// change the charge
  void  setCharge( double chg )            { itsCharge = chg; }
  /// change color information
  void  setColor( double col )             { itsColorCharge = col; }
  /// change spin information
  void  setSpin( const SpinState & spin )  { itsSpin = spin; }
  /// add a constituent particle
  void  addConstituent( Constituent c )    { itsQuarks.push_back( c ); }
  /// change the mass
  void  setMass( Measurement const & mass )        {  itsResonance.setMass(mass); }
  /// change the total width
  void  setTotalWidth( Measurement const & width ) {  itsResonance.setTotalWidth(width); }
  /// change the total width using a lifetime
  void  setTotalWidthFromLifetime( Measurement const & lt ) {  itsResonance.setTotalWidthFromLifetime(lt); }
  /// change the total width lower cutoff
  void  setLowerCutoff( double cut )               {  itsResonance.setLowerCutoff(cut); }
  /// change the total width upper cutoff
  void  setUpperCutoff( double cut )               {  itsResonance.setUpperCutoff(cut); }

private:

  // default constructor not allowed
  ParticleData();
 
  // ---  class-specific data:
  //
  ParticleID                itsID;
  std::string               itsParticleName;
  std::string               itsSource;
  int                       itsOriginalID;
  double                    itsCharge;
  double                    itsColorCharge;
  SpinState                 itsSpin;
  std::vector<Constituent>  itsQuarks;
  ResonanceStructure        itsResonance;

  /// helper method
  void getConstituentsFromPID( );

};  // ParticleData

inline void swap( ParticleData & first, ParticleData & second ) {
  first.swap( second );
}

}	// HepPDT

#include "HepPDT/ParticleData.icc"

#endif // PARTICLEDATA_HH
