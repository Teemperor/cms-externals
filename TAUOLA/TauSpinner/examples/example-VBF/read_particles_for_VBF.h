#ifndef READ_PARTICLES_FOR_VBF_H
#define READ_PARTICLES_FOR_VBF_H

#include <vector>


#include "Tauola/Log.h"
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"

#include "HepMC/IO_GenEvent.h"
#include "TauSpinner/SimpleParticle.h"
using TauSpinner::SimpleParticle;
using std::vector;

/** @brief Process HepMC event with events: ff -> H ff; H -> tau tau
 *
 *  Finds appropriate HepMC particles and returns them as SimpleParticle
 *  for TauSpinner
 */
int read_particles_for_VBF(HepMC::IO_GenEvent     &file,            // HepMC input file
                           SimpleParticle         &p1,              // incoming quark 1
                           SimpleParticle         &p2,              // incoming quark 2
                           SimpleParticle         &X,               // intermediate boson
                           SimpleParticle         &p3,              // jet 1
                           SimpleParticle         &p4,              // jet 2
                           SimpleParticle         &tau1,            // tau 1
                           SimpleParticle         &tau2,            // tau 2
                           vector<SimpleParticle> &tau1_daughters,  // tau 1 daughters
                           vector<SimpleParticle> &tau2_daughters); // tau 2 daughters

#endif
