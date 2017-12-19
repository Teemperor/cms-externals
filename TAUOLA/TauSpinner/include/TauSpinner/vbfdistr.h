#include "TauSpinner/SimpleParticle.h"
#include <vector>
using std::vector;

namespace TauSpinner {
extern "C" {
  extern struct {
    double sqrt__as;
    double g__exp__2;
    double conjg__ckm3x3;
    double ckm3x3;
    double   lamws__exp__2;
    double lamws__exp__3;
    double mz__exp__2;
    double mz__exp__4;
    double sqrt__2;
    double   mh__exp__2;
    double aew;
    double sqrt__aew;
    double ee;
    double mw__exp__2;
    double sw2;
    double cw;
    double sqrt__sw2;
    double sw;
    double g1;
    double  gw;
    double vev;
    double vev__exp__2;
    double lam;
    double yb;
    double yt;
    double ytau;
    double muh;
    double ee__exp__2;
    double sw__exp__2;
    double  cw__exp__2;
    double aewm1;
    double gf;
    double as;
    double lamws;
    double aws;
    double rhows;
    double etaws;
    double ymb;
    double ymt;
    double ymtau;
  } params_r_;

  extern struct {
    int icp; //tau pdg id
  } cpstatus_;
}
/**  vbf amplitudes re-initialization (after chanche of alpha_s */
 extern "C" { extern void vbf_reinit_(int *key);}

/** Choses variant for vbf amplitudes initializations reference and modified*/
 extern "C" { extern void vbfinit_(int *ref, int *variant);}

/** Definition of REAL*8 FUNCTION VBDISTR(I1,I2,I3,I4,H1,H2,P,KEY) from VBF_UD.f */
extern "C" double vbfdistr_(int *I1, int *I2, int *I3, int *I4, int *H1, int *H2, double P[6][4], int *KEY);

/** Wrapper to VBDISTR and entry point for vbfdistrModif*/
double vbfdistr(int I1, int I2, int I3, int I4, int H1, int H2, SimpleParticle &p1, SimpleParticle &p2, SimpleParticle &tau1, SimpleParticle &tau2, int KEY);


/** set option for QCD initialization */
void setPDFOpt(int QCDdefault, int QCDvariant);

/** Set vbfdistrModif function. Function arguments as of vbfdistr(), last one would be  its result*/
void set_vbfdistrModif(double (*function)(int, int, int, int, int, int, double[6][4], int, double) );

/** Set alphasModif function. Function arguments as of alphas() */
 void set_alphasModif(void (*function)(double, int, int) );

/** Get VBF ME2
    Returns array W[2][2] */
void getME2VBF(SimpleParticle &p3, SimpleParticle &p4, SimpleParticle &sp_X,SimpleParticle &tau1, SimpleParticle &tau2, double (&W)[2][2], int KEY);

double calculateWeightFromParticlesVBF(SimpleParticle &p3, SimpleParticle &p4,SimpleParticle &sp_X, SimpleParticle &sp_tau1, SimpleParticle &sp_tau2, vector<SimpleParticle> &sp_tau1_daughters, vector<SimpleParticle> &sp_tau2_daughters);

} // namespace TauSpinner
