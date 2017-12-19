#ifndef _TAU_REWEIGHT_LIB_H_
#define _TAU_REWEIGHT_LIB_H_

// Debug mode
#ifdef DEBUG_MODE
#define DEBUG(arg) arg;
#else
#define DEBUG(arg)
#endif

// TAUOLA header
#include "Tauola/Tauola.h"

#include "TauSpinner/Tauola_wrapper.h"
#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/Particle.h"

// LHAPDF header
#include "LHAPDF/LHAPDF.h"

#include <vector>
#include <iostream>
using std::vector;
using std::cout;
using std::endl;

namespace TauSpinner {

/** Definition of REAL*8 FUNCTION DISTH(S,T,H1,H2) from disth.f calculating 
    SM glue glue-->Higgs --> tau tau*/
extern "C" double disth_(double *SVAR, double *COSTHE, int *TA, int *TB);

/** Initialize TauSpinner

    Print info and set global variables */
void initialize_spinner(bool _Ipp, int _Ipol, int _nonSM2, int _nonSMN, double _CMSENE);

/** Set flag for calculating relative (NONSM/SM) or absolute that is matrix element^2 (spin averaged), 
  weight for X-section calculated as by product in longitudinal polarization method. */
void setRelWTnonSM(int _relWTnonSM);

/** set flag for simple formula (just Breit Wigner) of Higgs cross section, 
    set its parameters: Higgs mass, width and normalization for Higgs born (BW) function */
 void setHiggsParameters(int jak, double mass, double width, double normalization);

/** Get Higgs mass, width and normalization of Higgs born function */
void getHiggsParameters(double *mass, double *width, double *normalization);

/**  Set flag defining what type of spin treatment was used in analyzed  sample */
void setSpinOfSample(bool _Ipol);

/**  Turn  non Standard Model calculation  on/off */
void setNonSMkey(int key);

/** set transverse components for Higgs spin density matrix */
void setHiggsParametersTR(double Rxx, double Ryy, double Rxy, double Ryx);

/** set multipliers for transverse components  for Drell-Yan spin density matrix */
void setZgamMultipliersTR(double Rxx, double Ryy, double Rxy, double Ryx);

/** get transverse components  for Drell-Yan spin density matrix */
void getZgamParametersTR(double &Rxx, double &Ryy, double &Rxy, double &Ryx);

/** Get nonSM weight of production matrix element */
double getWtNonSM();

/** Get tau+ weight of its decay matrix element */
double getWtamplitP();

/** Get tau- weight of its decay matrix element */
double getWtamplitM();

/** Get tau Helicity
    Used after event  is analyzed to obtain  attributed tau helicity.
    Returns -1 or 1. Note that 0 is returned if attribution fails. Method
    assumes that sample has spin effects taken into account. Effects can 
    be taken into account with the help of spin weights. */
double getTauSpin();

/** Calculate weights.

  Determines decay channel, calculates all necessary components for 
  calculation of all weights, calculates weights.
  Function for W+/- and H+/-  */
double calculateWeightFromParticlesWorHpn(SimpleParticle &W, SimpleParticle &tau, SimpleParticle &nu_tau, vector<SimpleParticle> &tau_daughters);

/** Calculate weights.

  Determines decay channel, calculates all necessary components for 
  calculation of all weights, calculates weights.
  Function for H/Z */
double calculateWeightFromParticlesH(SimpleParticle &sp_X, SimpleParticle &sp_tau1, SimpleParticle &sp_tau2, vector<SimpleParticle> &sp_tau1_daughters, vector<SimpleParticle> &sp_tau2_daughters);

/** Prepare kinematics for HH calculation
  
  Boost particles to effective pair rest frame, and rotate them so that tau is on Z axis.
  Then rotate again with theta2 phi2 so neutrino (first tau decay product) from tau decay is along Z. */
void prepareKinematicForHH(Particle &tau, Particle &nu_tau, vector<Particle> &tau_daughters, double *phi2, double *theta2);

/** Calculate polarization vector.

  We use FORTRAN metdods to calculate HH.
  First decide what is the channel. After that, 4-vectors
  are moved to tau rest frame of tau.
  Polarimetric vector HH is then rotated using angles phi and theta.
  
  Order of the particles does not matter. */
double* calculateHH(int tau_pdgid, vector<Particle> &tau_daughters, double phi, double theta);

/**
 Get Longitudinal polarization
 
 Returns longitudinal polarization in Z/gamma* -> tau+ tau-
 S: invariant mass^2 of the bozon */
double getLongitudinalPolarization(double S, SimpleParticle &sp_tau, SimpleParticle &sp_nu_tau);

/** Check if the vector of particles match the listed order of pdgid's:

  1) Returns true if 'particles' contain all of the listed pdgid-s.
  2) If it does - 'particles' will be reordered so that they have
     the same order as listed pdgid-s.

  It is done so the order of particles is the same as the order mandatory
  for TAUOLA Fortran routines. */
bool channelMatch(vector<Particle> &particles, int p1, int p2=0, int p3=0, int p4=0, int p5=0, int p6=0);

//------------------------------------------------------------------------------
//-- Useful debug methods ------------------------------------------------------
//------------------------------------------------------------------------------

/** Prints out two vertices:
      W   -> tau, nu_tau
      tau -> tau_daughters */
void print(Particle &W, Particle &nu_tau, Particle &tau, vector<Particle> &tau_daughters);

/** Sums all 4-vectors of the particles on the list */
Particle *vector_sum(vector<Particle> &x);

} // namespace TauSpinner
#endif
