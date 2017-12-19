// ----------------------------------------------------------------------
//
// ResonanceStructure.hh
// Author:  Lynn Garren
//
// ResonanceStructure contains the minimum information for a Breit-Wigner
// distribution about a given mass.
// ----------------------------------------------------------------------
#ifndef RESONANCESTRUCTURE_HH
#define RESONANCESTRUCTURE_HH

#include <algorithm>	// swap()

#include "HepPDT/Measurement.hh"

namespace HepPDT {

//! The ResonanceStructure class is holds mass and width information.

///
/// \class ResonanceStructure
/// \author Lynn Garren
///
/// ResonanceStructure contains the minimum information for a Breit-Wigner
/// distribution about a given mass.
///
class ResonanceStructure  {

public:

  // ---  birth/death:
  //
  /// construct from mass and total width
  ResonanceStructure( Measurement mass = Measurement(),
                      Measurement width = Measurement(),
                      double min = 0., double max = 0. );
  virtual ~ResonanceStructure();

  // ---  copying:
  //
  ResonanceStructure( const ResonanceStructure & orig );
  ResonanceStructure& operator= ( const ResonanceStructure & rhs );
  void swap ( ResonanceStructure & other );

  // ---  accessors:
  //
  /// get the mass
  Measurement   const & mass()         const { return itsMass; }
  /// get the total width
  Measurement   const & totalWidth()   const { return itsTotalWidth; }
  /// calculate lifetime from total width
  Measurement   lifetime()             const;
  /// lower cutoff of allowed width values
  double        lowerCutoff()          const { return itsLowerCutoff; }
  /// upper cutoff of allowed width values
  double        upperCutoff()          const { return itsUpperCutoff; }

  // ---  mutators:
  //
  /// change the mass
  void  setMass( Measurement const & mass )        {  itsMass = mass; }
  /// change the total width
  void  setTotalWidth( Measurement const & width ) {  itsTotalWidth = width; }
  /// change the total width using a lifetime
  void  setTotalWidthFromLifetime( Measurement const & lt );
  /// change the lower cutoff of allowed width values
  void  setLowerCutoff( double cut )               {  itsLowerCutoff = cut; }
  /// change the upper cutoff of allowed width values
  void  setUpperCutoff( double cut )               {  itsUpperCutoff = cut; }

private:

  // ---  class-specific data:
  //
  Measurement   itsMass;
  Measurement   itsTotalWidth;
  double        itsLowerCutoff;
  double        itsUpperCutoff;

};  // ResonanceStructure

inline
void swap( ResonanceStructure & first, ResonanceStructure & second )  {
  first.swap( second );
}

}	// HepPDT

#endif // RESONANCESTRUCTURE_HH
