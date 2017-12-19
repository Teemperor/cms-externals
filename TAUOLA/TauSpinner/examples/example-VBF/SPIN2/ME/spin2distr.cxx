#include "../include/spin2distr.h"
#include <cstdio>
#include <cstring>

namespace SPIN2 {

/** Wrapper to VBDISTR and place for interface to user provided modification*/
double spin2distr(int I1, int I2, int I3, int I4, int H1, int H2, double P[6][4], int KEY, double vbfdistr_result)
{
  double P_copy[6][4];
  double original_result = 0.;
  
  memcpy(P_copy,P,sizeof(double)*6*4);

  if(  KEY<2) {
    // SM mode should not be used, SPIN2 is not prepared for that 
    return  spin2distr_(&I1, &I2, &I3, &I4, &H1, &H2, P_copy, &KEY); 
  }
  else  {
    // modification mode
    int KEY_BUF=KEY-2;
    return  spin2distr_(&I1, &I2, &I3, &I4, &H1, &H2, P_copy, &KEY);
  }

}

} // namespace SPIN2
