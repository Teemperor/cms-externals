#ifndef _TAU_SPINNER_NONSM_H_
#define _TAU_SPINNER_NONSM_H_
#include "TauSpinner/tau_reweight_lib.h"

namespace TauSpinner {

/** Default nonSM function

    Always returns 1.0. Prints warnings if used with key!=0 */
double default_nonSM_born(int ID, double S, double cost, int H1, int H2, int key);

/** Default nonSMH function

    Always returns 1.0. Prints warnings if used with key!=0 */
double default_nonSM_bornH(int ID, double S, double cost, int H1, int H2, int key);

/** nonSM born

    Switches between nonSM born for Z and H depending on 'IfHiggs' flag
    (global variable of TauSpinner namespace) */
double nonSM_born(int ID, double S, double cost, int H1, int H2, int key);

/** Set nonSM function

    Set function for user-defined  born, including new physics
    This function must be of format:

    double fun(int ID, double S, double cost, int H1, int H2, int key)

    Where: ID    -- 1 for down-type qqbar, 2 for up-type qqbar
           S     -- cm qqbar energy^2
           cost  -- cosTheta (scattering angle of tau- with respect to quark)
           H1,H2 -- spin state of first, second tau 
           key   -- 1 when new effect is added, 0  Standard Model 

    If set to NULL, default function will be used.      */
void set_nonSM_born( double (*fun)(int, double, double, int, int, int) );


/** Set nonSMH function

    Set function for user-defined  born, including new physics
    This function must be of format:

    double fun(int ID, double S, double cost, int H1, int H2, int key)

    Where: ID    -- 0 for glu-glu, 1 for down-type qqbar, 2 for up-type qqbar (may be dummy)
           S     -- cm qqbar energy^2
           cost  -- cosTheta (scattering angle of tau- with respect to gluon/quark)
           H1,H2 -- spin state of first, second tau 
           key   -- 1 when new effect is added, 0  Standard Model 

    If set to NULL, default function will be used.      */
void set_nonSM_bornH( double (*fun)(int, double, double, int, int, int) );


/** plzap2

    This function calculates probability for the helicity +1 +1 configuration
    of taus for given Z/gamma transfer and COSTH0 cosine of scattering angle.
  
    This is the copy of plzap0 using 'nonSM_born' instead
    of 't_born' when switch for nonSM is turned on. */
double plzap2(int ide, int idf, double svar, double costhe);


/** plweight

    This function calculates ratio of cross section if switch  nonSM2
    contribution is turned on. There are two options possible. One which
    change costhe distribution only and the one changing x-section as well */
double plweight(int ide, double svar, double costhe);


/** plnorm

    This function calculates ratio of integrated cross section if nonSMN
    is turned on. We integrate distr over cos-theta, assuming that distr is
    polynomial of the 4-th order in cos-theta (at most) */
double plnorm(int ide, double svar);


/** Calculate nonSM sig and polarization for Higgs

   Used only if nonSM2=1

   Input:   S - virtuality of H (note that S may not be ivariant of tau pair
                                 only, extra photons could contribute)
            tau1,tau2 - tau pair
   Returns: corrX2, polX2
            TauSpinner::WTnonSM */
void nonSMHcorrPol(double S, SimpleParticle &tau1, SimpleParticle &tau2,
                   double *corrX2, double *polX2);

} // namespace TauSpinner
#endif
