// ----------------------------------------------------------------------
//
// Measurement.hh
// author:  Mark Fischler
// ----------------------------------------------------------------------
#ifndef MEASUREMENT_HH
#define MEASUREMENT_HH

namespace HepPDT {

//! The Measurement class defines a value with its error

///
/// \class Measurement
/// \author Mark Fischler
///
/// Contain a value with associated errors.
/// Provide simple comparisons.
///
class Measurement {
public:
  inline Measurement ();  
  inline Measurement ( double value, double sigma );

  inline Measurement ( const Measurement &m );
  inline void swap( Measurement & other );
  inline Measurement& operator = ( Measurement const & rhs );

  /// compare the value, ignore the error
  inline bool  operator <  ( Measurement const & other ) const;
  /// compare the value, ignore the error
  inline bool  operator == ( Measurement const & other ) const;

  inline double   value()  const;
  inline double   sigma()  const;
  inline operator double() const;

private:  
  double val;
  double sig;

};	// Measurement

inline
void swap( Measurement & first, Measurement & second ) { first.swap( second ); }

}      // end of namespace HepPDT

#include "HepPDT/Measurement.icc"

#endif // end of codeguard MEASUREMENT_HH
