// ----------------------------------------------------------------------
//
// ParticleDataTableComparison.hh
// Author: Lynn Garren
//
// ----------------------------------------------------------------------
#ifndef ParticleDataTableComparison_HH
#define ParticleDataTableComparison_HH

#include "HepPDT/ParticleID.hh"

namespace HepPDT {
  
//! The ParticleDataTableComparison class provides a utility for sorting the PDT.

///
/// \class ParticleDataTableComparison
/// \author Lynn Garren
///
///  Sort the PDT by absolute value of the particle ID.
///  If the absolute values are equal, the positive number comes first.
///
class ParticleDataTableComparison  {

public:
  /// stateless class
  ParticleDataTableComparison( ) {;}
  
  /// Comparison method sorts by absolute value.
  /// If the absolute values are the same, the positive entry comes first.
  inline bool operator() ( const ParticleID&, const ParticleID & ) const;

};  // ParticleDataTableComparison


bool ParticleDataTableComparison::operator() ( const ParticleID& p1, const ParticleID & p2 ) const
{
    if( p1.abspid() == p2.abspid() ) {
      return  p2.pid() <  p1.pid();
    }
    return p1.abspid() < p2.abspid();
}


}	// HepPDT

#endif // ParticleDataTableComparison_HH
