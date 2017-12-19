// ----------------------------------------------------------------------
//
// SpinState.hh
// Author:  Lynn Garren
//
// ----------------------------------------------------------------------
#ifndef SPINSTATE_HH
#define SPINSTATE_HH

namespace HepPDT {

//! The SpinState class holds spin information

///
/// \class SpinState
/// \author Lynn Garren
///
/// SpinState contains total spin, spin, and orbital angular momentum.
///
class SpinState  {

public:

  // ---  birth/death:
  //
  /// (default) constructor
  inline SpinState( double ts = 0., double spin = 0., double oam =0. );

  // ---  copying:
  //
  inline SpinState( const SpinState & orig );
  inline SpinState& operator = ( const SpinState & rhs );
  inline void swap( SpinState & other );

  // ---  booleans:
  //
  /// all three spins must match
  inline bool  operator == ( const SpinState & rhs ) const;

  // ---  accessors:
  //
  /// return the total spin
  double totalSpin()     const {  return itsTotalSpin; }
  /// return the spin
  double spin()          const {  return itsSpin; }
  /// return the orbital angular momentum
  double orbAngMom()     const {  return itsOrbAngMom; }

  // ---  mutators:
  //
  /// change the total spin
  void  setTotalSpin( double spin )  {  itsTotalSpin = spin; }
  /// change the spin
  void  setSpin( double spin )       {  itsSpin = spin; }
  /// change the orbital angular momentum
  void  setOrbAngMom( double ang )   {  itsOrbAngMom = ang; }
  
private:

  // ---  class-specific data:
  //
  double  itsTotalSpin;
  double  itsSpin;
  double  itsOrbAngMom;

};  // SpinState

inline
void swap( SpinState & first,SpinState & second )  { first.swap( second ); }

}	// HepPDT

#include "HepPDT/SpinState.icc"

#endif // SPINSTATE_HH
