// ----------------------------------------------------------------------
//
// DefTable.hh
// Author: Lynn Garren
//
// ----------------------------------------------------------------------
#ifndef DEFTABLE_HH
#define DEFTABLE_HH

#include <string>
#include <map>

namespace HepPDT {

//! The DefTable class holds EvtGen definitions

///
/// \class DefTable
/// \author Lynn Garren
///
/// This is temporary information storage used when reading EvtGen input.
///
class DefTable  {
  typedef  std::map<std::string,double>            TempDefMap;

public:
  typedef  TempDefMap::const_iterator  const_iterator;
  typedef  TempDefMap::iterator        iterator;

  DefTable( );
  ~DefTable()  { ; }

  // --- mutator
  //
  /// add a definition to the map
  void addDefinition( std::string const & def, double val ) { tempDefines[def] = val; }

  // --- boolean
  //
  /// is this definition already defined?
  bool hasDefinition( std::string const & def ) const;

  // --- accessors
  //
  /// get the size of the definition map
  int             size()  const { return tempDefines.size(); }
  /// return the definition of this parameter
  double          definition( std::string const & def );
  /// use for diagnostics
  void            writeDefinitions() const;	// intended for diagnostic use

  /// begin iterating over the definition map
  iterator        begin()       { return tempDefines.begin(); }
  /// begin iterating over the definition map
  const_iterator  begin() const { return tempDefines.begin(); }

  /// end iterating over the definition map
  iterator        end()       { return tempDefines.end(); }
  /// end iterating over the definition map
  const_iterator  end() const { return tempDefines.end(); }

private:
  TempDefMap           tempDefines;	// hold EvtGen "Define" info

  // forbidden copy:
  DefTable &  operator = ( DefTable const & );
  DefTable( DefTable const & );

};  // DefTable<>

}  // namespace HepPDT

#endif // DEFTABLE_HH
