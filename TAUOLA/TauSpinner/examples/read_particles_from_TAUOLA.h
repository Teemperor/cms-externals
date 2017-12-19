#ifndef _READ_PARTICLES_FROM_TAUOLA_H_
#define _READ_PARTICLES_FROM_TAUOLA_H_

// HepMC headers
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/IO_GenEvent.h"

#include "TauSpinner/SimpleParticle.h"
#include <vector>
using namespace std;
using namespace TauSpinner;

/** Get HepMC::Event

  Return last event parsed by readParticlesFromTAUOLA_HepMC */
HepMC::GenEvent* readParticlesFromTAUOLA_HepMC_getEvent();

/** Read HepMC::IO_GenEvent.

  Read HepMC event from data file
  and return particles needed for tau spin weight calculation.
  
  This routine is prepared for use with files generated by Pythia8.
  Fills:
  
  'X'              - Heavy particle (W+/-, H+/-, H, Z)
  'tau'            - first tau
  'tau2'           - second tau or nu_tau, if 'X' decays to one tau only
  'tau_daughters'  - daughters of 'tau'
  'tau2_daughters' - daughters of 'tau2' or empty list, if 'tau2' is nu_tau.
  
  Returns:
  0 - no more events to read               (finished processing the file)
  1 - no decay found in the event          (finished processing current event)
  2 - decay found and processed correctly.
      Event will continue to be processed
      with next function call. */
int readParticlesFromTAUOLA_HepMC(HepMC::IO_GenEvent &input_file, SimpleParticle &X, SimpleParticle &tau, SimpleParticle &tau2, vector<SimpleParticle> &tau_daughters, vector<SimpleParticle> &tau2_daughters);

/** Get daughters of HepMC::GenParticle

  Recursively searches for final-state daughters of 'x' */
vector<SimpleParticle> *getDaughters(HepMC::GenParticle *x);

#endif
