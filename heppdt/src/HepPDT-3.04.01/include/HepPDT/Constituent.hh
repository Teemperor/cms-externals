// ----------------------------------------------------------------------
//
// Constituent.hh
// Author:  Lynn Garren
//
// Holds a particle constituent (e.g. ParticleID of the quark type)
// ----------------------------------------------------------------------
#ifndef CONSTITUENT_HH
#define CONSTITUENT_HH

#include <algorithm>	// swap()

#include "HepPDT/ParticleID.hh"

namespace HepPDT {

//! The Constituent class has information about constituent particles.

///
/// \class Constituent
/// \author Lynn Garren
///
/// Holds a particle constituent 
/// (e.g. quark type and number of quarks of this type)
///
class Constituent  {

public:
  // ---  birth/death:
  //
  Constituent( ParticleID p = ParticleID(0), int m = -1 )
    : itsPid(p), itsMultiplicity(m)  {}
    
  // ---  copying:
  //
  Constituent( Constituent const & orig )
    :  itsPid(orig.itsPid), itsMultiplicity(orig.itsMultiplicity)  {}
  Constituent &  operator = ( Constituent const & rhs ) {
    Constituent temp( rhs );
    swap( temp );
    return *this;
  }
  void swap( Constituent & other ) {
    std::swap(itsPid, other.itsPid );
    std::swap( itsMultiplicity, other.itsMultiplicity );
  }

  // ---  accessors:
  /// how many of this constituent are there?
  int         multiplicity() const  { return itsMultiplicity; }
  /// ParticleID of this constituent
  ParticleID  pid()          const  { return itsPid; }

  // ---  booleans:
  /// is this an up quark?
  bool  isUp() const;
  /// is this a down quark?
  bool  isDown() const;
  /// is this a strange quark?
  bool  isStrange() const;
  /// is this a charm quark?
  bool  isCharm() const;
  /// is this a bottom quark?
  bool  isBottom() const;
  /// is this a top quark?
  bool  isTop() const;

private:
  // -- data members
  ParticleID  itsPid;
  int         itsMultiplicity;

};  // Constituent


inline
void swap( Constituent & first, Constituent & second ) { first.swap( second ); }

}	// HepPDT

#endif // CONSTITUENT_HH
