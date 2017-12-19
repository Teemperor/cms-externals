#ifndef _TAUOLA_WRAPPER_
#define _TAUOLA_WRAPPER_
/** 
  Wrapper for tauola routines calculating tau polarimetric vector HH.
  Note that order of input 4-vectors matters and must be prepared
  accordingly to choice in FORTRAN version of TAUOLA.
  For details see documentation of old TAUOLA fortran.

  Methods available for some decay modes only.
*/

#include "Tauola/f_Variables.h"

namespace TauSpinner {

extern "C" {
  void dam2pi_(int *MNUM, float *PT, float *PN, float *PIM1, float *PIM2, float *AMPLIT, float *HH);

  void dam4pi_(int *MNUM, float *PT, float *PN, float *PIM1, float *PIM2, float *PIZ, float *PIPL,
               float *AMPLIT, float *HH);

  void damppk_(int *MNUM, float *PT, float *PN, float *PIM1, float *PIM2, float *PIPL,
               float *AMPLIT, float *HH);

  void dampry_(int *ITDKRC, double *XK0DEC, double *XK, double *XA, double *QP, double *XN,
               double *AMPLIT, double *HV);

  void initwk_(int *IDE, int *IDF, double *SVAR);

  double t_born_(int *MODE, double *SVAR, double *COSTHE, double *TA, double *TB);
  
  // COMMON /CHANOPT/ used to switch between pi- pi- pi+ and pi0 pi0 pi-
  //                  in RChL currents
  extern "C" struct CHANOPT {
    int JJ;
  } chanopt_;

}

} // namespace TauSpinner
#endif
