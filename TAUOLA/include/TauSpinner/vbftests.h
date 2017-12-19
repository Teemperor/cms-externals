#include "TauSpinner/SimpleParticle.h"
#include <vector>
using std::vector;

namespace TauSpinner {

/** Definition of REAL*8 FUNCTION VBDISTR(I1,I2,I3,I4,H1,H2,P,KEY) from VBF_UD.f */
extern "C" double vbfdistr_(int *I1, int *I2, int *I3, int *I4, int *H1, int *H2, double P[6][4], int *KEY);

/** Simple test
    only pritnout */
void makeSimpleTestME2();
 void calcTestME2(int iter, double P[6][4]);

/** Simple test for PDFs
    only printout */
void makeSimpleTestPDF();

/** calcXsect
    Returns array W[2][2] */
 void calcXsect(int IDPROD, SimpleParticle &p3, SimpleParticle &p4, SimpleParticle &sp_X,SimpleParticle &tau1, SimpleParticle &tau2, double (&W)[2][2], int KEY);

/** calcProdMatrix
    Returns array W[2][2] */
 void calcProdMatrix(SimpleParticle &p3, SimpleParticle &p4, SimpleParticle &sp_X,SimpleParticle &tau1, SimpleParticle &tau2, double (&W)[2][2], int KEY, int ID1, int ID2, int ID3, int ID4, int pdfOpt);

/** calcSumME2
    Returns array W[2][2] */
 void calcSumME2(SimpleParticle &p3, SimpleParticle &p4, SimpleParticle &sp_X,SimpleParticle &tau1, SimpleParticle &tau2, double (&W)[2][2], int KEY, int ID1, int ID2, int ID3, int ID4);

/** calcPDFs
    Returns array W[2][2] */
 void calcPDFs(SimpleParticle &p3, SimpleParticle &p4, SimpleParticle &sp_X,SimpleParticle &tau1, SimpleParticle &tau2, double (&W)[2][2], int KEY, int ID1, int ID2, int ID3, int ID4, int pdfOpt);


 double calculateWeightFromParticlesVBFPROD(int IDPROD, SimpleParticle &p3, SimpleParticle &p4,SimpleParticle &sp_X, SimpleParticle &sp_tau1, SimpleParticle &sp_tau2, vector<SimpleParticle> &sp_tau1_daughters, vector<SimpleParticle> &sp_tau2_daughters, int KEY);


} // namespace TauSpinner
