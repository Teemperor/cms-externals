// ----------------------------------------------------------------------
//
// hasMethods.cc
// Author: Lynn Garren
//
//  check to see if this ParticleData has a particular constituent
//  look at Constituents, not PID
//
// ----------------------------------------------------------------------

#include "HepPDT/defs.h"
#include "HepPDT/ParticleData.hh"

namespace HepPDT {

bool ParticleData::hasUp( ) const
{
    unsigned int i;
    if( itsQuarks.size() == 0 ) { return false; }
    for( i=0; i<itsQuarks.size(); ++i ) {
       if( itsQuarks[i].isUp() ) { return true; }
    }
    return false;
}

bool ParticleData::hasDown( ) const
{
    unsigned int i;
    if( itsQuarks.size() == 0 ) { return false; }
    for( i=0; i<itsQuarks.size(); ++i ) {
       if( itsQuarks[i].isDown() ) { return true; }
    }
    return false;
}

bool ParticleData::hasStrange( ) const
{
    unsigned int i;
    if( itsQuarks.size() == 0 ) { return false; }
    for( i=0; i<itsQuarks.size(); ++i ) {
       if( itsQuarks[i].isStrange() ) { return true; }
    }
    return false;
}

bool ParticleData::hasCharm( ) const
{
    unsigned int i;
    if( itsQuarks.size() == 0 ) { return false; }
    for( i=0; i<itsQuarks.size(); ++i ) {
       if( itsQuarks[i].isCharm() ) { return true; }
    }
    return false;
}

bool ParticleData::hasBottom( ) const
{
    unsigned int i;
    if( itsQuarks.size() == 0 ) { return false; }
    for( i=0; i<itsQuarks.size(); ++i ) {
       if( itsQuarks[i].isBottom() ) { return true; }
    }
    return false;
}

bool ParticleData::hasTop( ) const
{
    unsigned int i;
    if( itsQuarks.size() == 0 ) { return false; }
    for( i=0; i<itsQuarks.size(); ++i ) {
       if( itsQuarks[i].isTop() ) { return true; }
    }
    return false;
}

} // HepPDT
