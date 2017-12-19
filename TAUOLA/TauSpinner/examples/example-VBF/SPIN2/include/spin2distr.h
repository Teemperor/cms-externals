namespace SPIN2 {

/** Wrapper to VBDISTR and place for interface to user provided modification*/
double spin2distr(int I1, int I2, int I3, int I4, int H1, int H2, double P[6][4], int KEY, double vbfdistr_result);

}

// extern "C" { extern void spin2_reinit_(int *key);}
/** Choses variant for spin2 amplitudes initializations reference and modified*/
extern "C" { extern void spin2init_(int *ref, int *variant);}

/** Definition of REAL*8 FUNCTION SPIN2DISTR(I1,I2,I3,I4,H1,H2,P,KEY) from SPIN2_UD.f */
extern "C" double spin2distr_(int *I1, int *I2, int *I3, int *I4, int *H1, int *H2, double P[6][4], int *KEY);
