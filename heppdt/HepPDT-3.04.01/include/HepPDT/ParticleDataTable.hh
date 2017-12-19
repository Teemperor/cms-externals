// ----------------------------------------------------------------------
//
// ParticleDataTable.hh
// Author: Lynn Garren, Walter Brown
//
// ----------------------------------------------------------------------
#ifndef ParticleDataTable_HH
#define ParticleDataTable_HH

#include <iostream>
#include <string>
#include <map>

#include "HepPDT/ParticleID.hh"
#include "HepPDT/ParticleData.hh"
#include "HepPDT/ParticleDataTableComparison.hh"
#include "HepPDT/ProcessUnknownID.hh"
#include "HepPDT/SimpleProcessUnknownID.hh"
#include "HepPDT/Version.hh"

//! HepPDT is a Particle Data Table namespace 

///
/// \namespace HepPDT
/// The HepPDT classes are used to create a Particle Data Table.
///
namespace HepPDT {
  
//! The ParticleDataTable class is the core of HepPDT.

///
/// \class ParticleDataTable
/// \author Lynn Garren, Walter Brown
///
/// This is the table of particle data information.
/// This table is created once at the beginning of a job and referenced
/// as needed.  
/// We expect that the table will be saved as part of the standard job output.
/// Methods exist to fill the table from a variety of input formats.
/// The user may fill the table from as many differernt input streams 
/// as desired.
/// See the examples.
///
class ParticleDataTable  {
  
public:
  typedef  ParticleData                     CPD;

  typedef  std::map<ParticleID,TempParticleData>  TempMap;
  typedef  std::map<ParticleID,ParticleData,ParticleDataTableComparison> PDTMap;
  typedef  std::map<std::string,ParticleID>       PDTNameMap;

  typedef PDTMap::const_iterator                  const_iterator;
  typedef PDTNameMap::const_iterator              const_iteratorByName;

  // ---  birth/death:
  //
  /// Construct a particle data table with an identifying name.
  /// Require a method to deal with unknown PID's.
  ParticleDataTable( std::string name = " ", ProcessUnknownID* = new SimpleProcessUnknownID );
  ~ParticleDataTable();

  // ---  accessors:
  //
  /// size of the particle data table
  int             size()  const { return itsMap.size(); }
  /// begin iterating over the particle data table
  const_iterator  begin() const { return itsMap.begin(); }
  /// end iterating over the particle data table
  const_iterator  end()   const { return itsMap.end(); }
  /// size of the map of particle names
  int                   sizeNameMap()  const { return itsNameMap.size(); }
  /// begin iterating over the map of particle names
  const_iteratorByName  beginNameMap() const { return itsNameMap.begin(); }
  /// end iterating over the map of particle names
  const_iteratorByName  endNameMap()   const { return itsNameMap.end(); }
  /// return the name of this particle data table
  std::string     tableName() const { return itsTableName; }

  /// access particle information via ParticleID
  ParticleData const * particle( ParticleID ) const;
  /// access particle information via ParticleID
  ParticleData       * particle( ParticleID );
  /// access particle information via a particle name
  ParticleData const * particle( std::string ) const;
  /// access particle information via a particle name
  ParticleData       * particle( std::string );

  /// access particle information via ParticleID
  ParticleData       * operator [] ( ParticleID );
  /// access particle information via ParticleID
  ParticleData const * operator [] ( ParticleID ) const;
  /// access particle information via a particle name
  ParticleData       * operator [] ( std::string );
  /// access particle information via a particle name
  ParticleData const * operator [] ( std::string ) const;

  /// output all information in the PDT
  void writeParticleData( std::ostream & outstr );
  /// output all information about a particle EXCEPT its decays
  void writeParticleInfo( std::ostream & outstr );
  /// output a list of original IDs and their translations
  void writeParticleTranslation( std::ostream & outstr );
  /// primarily useful for testing
  void writeParticleStatus( std::ostream & );
  
  /// used by the TableBuilder destructor to fill the PDT 
  void convertTemporaryMap( TempMap & tempPDT, std::ostream & err );

private:
  // --- local data:
  mutable PDTMap      itsMap;
  mutable PDTNameMap  itsNameMap;
  std::string itsTableName;
  ProcessUnknownID * itsProcessUnknownID;

  // ---  copying; forbidden:
  ParticleDataTable( const ParticleDataTable & orig );
  ParticleDataTable& operator=( const ParticleDataTable & );
  
  // --- private methods
  void addParticle( ParticleData const & p ) const;
  // private interface to ProcessUnknownID
  ParticleData * DealWithUnknownID( ParticleID & ) const;
  // private non const iterators
  typedef PDTMap::iterator        iterator;
  typedef PDTNameMap::iterator    nameIterator;

};  // ParticleDataTable

// --- free functions
//

bool writePDGStream (std::ostream & os, const ParticleDataTable & table);

bool writePythiaStream (std::ostream & os, const ParticleDataTable & table);

bool writeHerwigStream (std::ostream & os, const ParticleDataTable & table);

bool writeIsajetStream (std::ostream & os, const ParticleDataTable & table);

bool writeQQStream (std::ostream & os, const ParticleDataTable & table);

bool writeEvtGenStream (std::ostream & os, const ParticleDataTable & table);


}	// HepPDT

#include "HepPDT/ParticleDataTable.icc"

#endif // ParticleDataTable_HH
