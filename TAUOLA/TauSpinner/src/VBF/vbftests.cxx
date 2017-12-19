#include "TauSpinner/tau_reweight_lib.h"
#include "TauSpinner/Particle.h"
#include "TauSpinner/vbftests.h"
#include <cstdlib>
#include <string.h>
using namespace Tauolapp;

/*This is a group of methods prepared for the testing purpose. The tests include those where fixed kinematics is 
forced. Then external calculation can be used for matrix element calculation and/or PDFs etc. 
Gradually better documentation will be introduced. It will be documented with comment lines and in the draft of the paper.
To use one has to modify the Makefile the appropriate line of 
objects has to be activated (# removed). Then invoking lines in code have to be commented out as well.
At a certain moment such commented  out lines will be removed from the code or will remain. We will have to make up the mind.
-*/

namespace TauSpinner {

extern double CMSENE;
extern int    relWTnonSM;
extern double WTnonSM;   
extern double Polari;
extern double WTamplit;
extern double WTamplitP;
extern double WTamplitM;
extern double CMSENE;
extern double f(double x, int ID, double SS, double cmsene);

const bool  DEBUG = 0;

/** Wrapper to VBDISTR */
double vbfdistrWRAP(int I1, int I2, int I3, int I4, int H1, int H2, double P[6][4], int KEY)
{
  double P_copy[6][4];
  memcpy(P_copy,P,sizeof(double)*6*4);

  return vbfdistr_(&I1, &I2, &I3, &I4, &H1, &H2, P_copy, &KEY);

}

//---------------------------------------------------------------------------------------------------------------------------------------------------
  void calcXsect(int IDPROD, SimpleParticle &p3, SimpleParticle &p4, SimpleParticle &sp_X,SimpleParticle &tau1, SimpleParticle &tau2, double (&W)[2][2], int KEY)
{
  /*
  // this may be necessary because of matrix element calculation may require absolute energy-momentum conservation!
  // FSR photons may need to be treated explicitely or with interpolation procedures.
  Particle P_QQ( p3.px()+p4.px()+tau1.px()+tau2.px(),
                 p3.py()+p4.py()+tau1.py()+tau2.py(),
                 p3.pz()+p4.pz()+tau1.pz()+tau2.pz(),
                 p3.e() +p4.e() +tau1.e() +tau2.e(), 0 );
  */
  Particle P_QQ( p3.px()+p4.px()+sp_X.px(),
                 p3.py()+p4.py()+sp_X.py(),
                 p3.pz()+p4.pz()+sp_X.pz(),
                 p3.e() +p4.e() +sp_X.e() , 0 );


  double SS = P_QQ.recalculated_mass()*P_QQ.recalculated_mass(); 
  
  double x1x2  = SS/CMSENE/CMSENE;
  double x1Mx2 = P_QQ.pz()/CMSENE*2;
  
  double x1 = (  x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;
  double x2 = ( -x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;

  //---------------------------------------------------------------------------
  // Construct the matrix for FORTRAN function
  // NOTE: different order of indices than in FORTRAN!
  double P[6][4] = { { CMSENE/2*x1, 0.0,       0.0,          CMSENE/2*x1  },
                     { CMSENE/2*x2, 0.0,       0.0,         -CMSENE/2*x2  }, 
                     { p3.e(),      p3.px(),   p3.py(),      p3.pz()      }, 
                     { p4.e(),      p4.px(),   p4.py(),      p4.pz()      },
                     { tau1.e(),    tau1.px(), tau1.py(),    tau1.pz()    },
                     { tau2.e(),    tau2.px(), tau2.py(),    tau2.pz()    } };
 

  W[0][0]=0.0;
  W[0][1]=0.0;
  W[1][0]=0.0;
  W[1][1]=0.0;
  
  //
  // Calculate 'f' function for all x1 and all ID1, ID2
  //
  double  f_x1_ARRAY[11] = { 0.0 };
  double  f_x2_ARRAY[11] = { 0.0 };
  double *f_x1 = f_x1_ARRAY+5;     // This way f_x1[i],f_x2[i] can be used with 
  double *f_x2 = f_x2_ARRAY+5;     // i going from -5 to 5

  double SSfix = 91.188*91.188;
  double CMSENEfix = 91.188;
  for(int i=-5;i<=5;++i) {
    f_x1[i] = f(x1,i,SSfix,CMSENEfix);
    f_x2[i] = f(x2,i,SSfix,CMSENEfix);
    //    std::cout << "ERW: x, Q2, f(x,Q2), f(x,Q)=" << x1 << "  " << SS << "  " << f(x1,i,SS,CMSENE) << "  " << f(x1,i,sqrt(SS),CMSENE) << std::endl;
  }

  W[0][0]=0.0;
  W[0][1]=0.0;
  W[1][0]=0.0;
  W[1][1]=0.0;
  
  // these loops need to be cleaned from zero contributions! 
  for(int I1=-4;I1<=4;I1++){  // for test of single production process fix flavour
    for(int I2=-4;I2<=4;I2++){  // for test of single production process fix flavour
      for(int I3=-4;I3<=4;I3++){
	for(int I4=-4;I4<=4;I4++){
	  
          int ID1 = I1; if( ID1 == 0 ) ID1 = 21;
          int ID2 = I2; if( ID2 == 0 ) ID2 = 21;
          int ID3 = I3; if( ID3 == 0 ) ID3 = 21;
          int ID4 = I4; if( ID4 == 0 ) ID4 = 21;

	  //preselect production group
	  bool accept = false;
	  // any process 
          if( IDPROD == 0)  accept = true;
          // only GG processes
          if( IDPROD == 1 &&  ID1 == 21 && ID2 == 21 && abs (ID3) < 5 &&  abs (ID4) < 5 ) accept = true;
          // only GG processes ... details
          if( IDPROD == 110 &&  ID1 == 21 && ID2 == 21 ) {
	    accept = true;
            if( abs (ID3) != 2 || abs (ID4) != 2  ) accept = false;
	  }
          if( IDPROD == 111 &&  ID1 == 21 && ID2 == 21 ) {
	    accept = true;
            if( abs (ID3) != 1 || abs (ID4) != 1  ) accept = false;
	  }
          if( IDPROD == 112 &&  ID1 == 21 && ID2 == 21 ) {
	    accept = true;
            if( abs (ID3) != 3 || abs (ID4) != 3  ) accept = false;
	  }
          if( IDPROD == 113 &&  ID1 == 21 && ID2 == 21 ) {
	    accept = true;
            if( abs (ID3) != 4 || abs (ID4) != 4  ) accept = false;
	  }
          // only GQ processes
          if( IDPROD == 2 &&  ID1 == 21 && abs(ID2) < 5 ) accept = true ;
          if( IDPROD == 2 &&  ID2 == 21 && abs(ID1) < 5 ) accept = true ;
          // only GQ processes ... details
          if( (IDPROD == 210) && ( ID1 == 21) && ( ID2 < 5 ) && (ID2 > 0) ) accept = true;
          if( (IDPROD == 211) && ( ID1 == 21) && ( ID2 < 5 ) && (ID2 < 0) ) accept = true;
	  if( (IDPROD == 212) && ( ID2 == 21) && ( ID1 < 5 ) && (ID1 > 0) ) accept = true;
	  if( (IDPROD == 213) && ( ID2 == 21) && ( ID1 < 5 ) && (ID1 < 0) ) accept = true;
         // only QQ, QXQX processes
	  if( IDPROD == 3 &&  (ID1 * ID2 > 0 ) && abs (ID1) < 5 &&  abs (ID2) < 5  ) accept = true;
         // only QQ processes .. details
 	  if( IDPROD == 31 &&  ID1 > 0  &&  ID2 > 0  && abs (ID1) < 5 &&  abs (ID2) < 5  ) accept = true;
 	  if( IDPROD == 32 &&  ID1 < 0  &&  ID2 < 0  && abs (ID1) < 5 &&  abs (ID2) < 5  ) accept = true;
         // only QQ processes .. details
	  if( IDPROD == 310 &&  ID1 > 0 &&  ID2 > 0  && ID1 < 5 && ID2 < 5 && ID1 == ID2 ) accept = true;
	  if( IDPROD == 311 &&  ID1 > 0 &&  ID2 > 0  && ID1 < 5 && ID2 < 5 && ID1 != ID2 ) accept = true;
	  if( IDPROD == 312 &&  ID1 < 0 &&  ID2 < 0  && ID1 == ID2 ) accept = true;
	  if( IDPROD == 313 &&  ID1 < 0 &&  ID2 < 0  && ID1 != ID2 ) accept = true;
          // only QQ processes .. details
	  if( IDPROD == 320 &&  ID1 *  ID2 > 0  && abs(ID1) < 5 && abs(ID2) < 5 && abs(ID1) == abs(ID2) ) accept = true;
	  if( IDPROD == 321 &&  ID1 *  ID2 > 0  && abs(ID1) < 5 && abs(ID2) < 5 && abs(ID1) != abs(ID2) ) accept = true;
          // only QQ processes .. details
          if( IDPROD == 3101 &&  ID1 ==  2 &&  ID2 ==  2 &&  ID3 ==  2 &&  ID4 ==  2 ) accept = true;
          if( IDPROD == 3102 &&  ID1 == -2 &&  ID2 == -2 &&  ID3 == -2 &&  ID4 == -2  ) accept = true;
          if( IDPROD == 3103 &&  ID1 ==  2 &&  ID2 ==  2 ) accept = true;
          if( IDPROD == 3104 &&  ID1 == -2 &&  ID2 == -2 ) accept = true;
        // only QQX processes
          if( IDPROD == 4 &&  (ID1 * ID2 < 0 ) && abs (ID1) < 5 &&  abs (ID2) < 5  ) accept = true;
          // only QQX processes .. details
	  if( IDPROD == 420 &&  ID1 * ID2 < 0  && abs(ID1) < 5 && abs(ID2) < 5 && abs(ID1) == abs(ID2) ) accept = true;
	  if( IDPROD == 421 &&  ID1 * ID2 < 0  && abs(ID1) < 5 && abs(ID2) < 5 && abs(ID1) != abs(ID2) ) accept = true;
          // only QQX processes, first family in initial/final
          if( IDPROD == 41 &&  (ID1 * ID2 < 0 ) && abs(ID1) < 5 &&  abs(ID2) < 5 && abs(ID3) < 5 &&  abs(ID4) < 5 ) {
	    accept = true;
            if( abs (ID1) > 2 || abs (ID2) > 2 || abs (ID3) >2 || abs (ID4) > 2  ) accept = false;
	  }
          // only QQX processes, second family in initial/final state
          if( IDPROD == 42 &&  (ID1 * ID2 < 0 ) && abs(ID1) < 5 &&  abs(ID2) < 5 && abs(ID3) < 5 &&  abs(ID4) < 5 ) {
	    accept = true;
            if( abs (ID1) < 3 || abs (ID2) < 3 || abs (ID3) < 3 || abs (ID4) < 3  ) accept = false;
	  }
          // only QQX processes, mixed families in initial/final state
          if( IDPROD == 43 &&  ( ID1 * ID2 < 0) && abs(ID1) < 5 &&  abs(ID2) < 5 && abs(ID3) < 5 &&  abs(ID4) < 5 ) {
	    accept = true;
            if(  abs (ID1) < 3  &&  abs (ID2) < 3  &&  abs (ID3) < 3  &&  abs (ID4) < 3  ) accept = false;
            if(  abs (ID1) > 2  &&  abs (ID2) > 2  &&  abs (ID3) > 2  &&  abs (ID4) > 2  ) accept = false; 
	  }
          // only QQX processes, first family in initial state
          if( IDPROD == 44 &&  (ID1 * ID2 < 0 ) && abs(ID1) < 5 &&  abs(ID2) < 5 ) {
	    accept = true;
            if( abs (ID1) > 2 || abs (ID2) > 2 ) accept = false;
	  }
          // only QQX processes ... details
          if( IDPROD == 4101 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> U UX
            if( abs (ID1) != 2 || abs (ID2) != 2 || abs (ID3) != 2 || abs (ID4) != 2  ) accept = false;
            if( ID1 != 2 ) accept = false;
	  }
          if( IDPROD == 4102 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> U UX
            if( abs (ID1) != 2 || abs (ID2) != 2 || abs (ID3) != 2 || abs (ID4) != 2  ) accept = false;
            if( ID1 != -2 ) accept = false;
	  }
          if( IDPROD == 4103&&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> Q QX
            if( abs (ID1) != 2 || abs (ID2) != 2 ) accept = false;
            if( ID1 != 2 ) accept = false;
	  }
          if( IDPROD == 4104&&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only UX U --> Q QX
            if( abs (ID1) != 2 || abs (ID2) != 2 ) accept = false;
            if( ID1 != -2 ) accept = false;
	  }
          if( IDPROD == 4105&&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only D DX --> Q QX
            if( abs (ID1) != 1 || abs (ID2) != 1 ) accept = false;
            if( ID1 != 1 ) accept = false;
	  }
          if( IDPROD == 4106&&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only DX D --> Q QX
            if( abs (ID1) != 1 || abs (ID2) != 1 ) accept = false;
            if( ID1 != -1 ) accept = false;
	  }
          if( IDPROD == 4107&&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only Q QX --> Q QX
            if( ID1 < 0  ||  ID2 > 0) accept = false;
	  }
          if( IDPROD == 4108&&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only QX Q --> Q QX
            if( ID1 > 0  ||  ID2 < 0) accept = false;
	  }
           // only QQX processes ... details
          if( IDPROD == 410 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> U UX
            if( abs (ID1) != 2 || abs (ID2) != 2 || abs (ID3) != 2 || abs (ID4) != 2  ) accept = false;
	  }
          if( IDPROD == 411 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> C CX
            if( abs (ID1) != 2 || abs (ID2) != 2 || abs (ID3) != 4 || abs (ID4) != 4  ) accept = false;
	  }
          if( IDPROD == 412 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> S DX
            if( abs (ID1) != 2 || abs (ID2) != 2 || abs (ID3) != 3 || abs (ID4) != 1  ) accept = false;
	  }
          if( IDPROD == 413 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> D SX
            if( abs (ID1) != 2 || abs (ID2) != 2 || abs (ID3) != 1 || abs (ID4) != 3  ) accept = false;
	  }
          if( IDPROD == 414 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> D DX
            if( abs (ID1) != 2 || abs (ID2) != 2 || abs (ID3) != 1 || abs (ID4) != 1  ) accept = false;
	  }
          if( IDPROD == 415 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> S SX
             if( abs (ID1) != 2 || abs (ID2) != 2 || abs (ID3) != 3 || abs (ID4) != 3  ) accept = false;
	  }
          if( IDPROD == 416 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> Q QX 
            if( abs (ID1) != 2 || abs (ID2) != 2 || abs (ID3) != 21 || abs (ID4) != 21  ) accept = false;
	  }
          if( IDPROD == 417 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> G G 
            if( abs (ID1) != 2 || abs (ID2) != 2 || abs (ID3) == 21 || abs (ID4) == 21  ) accept = false;
	  }
          if( IDPROD == 418 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> J J  
            if( abs (ID1) != 2 || abs (ID2) != 2 ) accept = false;
	  }
          if( IDPROD == 419 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only D DX --> D DX  
            if( abs (ID1) != 1 || abs (ID2) != 1 || abs (ID3) != 1 || abs (ID4) != 1  ) accept = false;
	  }
          if( IDPROD == 511 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only D DX, DX D --> G G 
            if( abs (ID1) != 1 || abs (ID2) != 1 || abs (ID3) != 21 || abs (ID4) != 21  ) accept = false;
	  }
          if( IDPROD == 512 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only D DX --> G G 
            if( ID1 != 1 || ID2 != -1 || abs (ID3) != 21 || abs (ID4) != 21  ) accept = false;
	  }
          if( IDPROD == 513 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only DX D --> G G 
            if( ID1 != -1 || ID2 != 1 || abs (ID3) != 21 || abs (ID4) != 21  ) accept = false;
	  }
          if( IDPROD == 521 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX, UX U --> G G 
            if( abs (ID1) != 2 || abs (ID2) != 2 || abs (ID3) != 21 || abs (ID4) != 21  ) accept = false;
	  }
          if( IDPROD == 522 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only U UX --> G G 
            if( ID1 != 2 || ID2 != -2 || abs (ID3) != 21 || abs (ID4) != 21  ) accept = false;
	  }
          if( IDPROD == 523 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only UX U --> G G 
            if( ID1 != -2 || ID2 != 2 || abs (ID3) != 21 || abs (ID4) != 21  ) accept = false;
	  }
          if( IDPROD == 531 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only Q QX, QX Q --> G G 
            if( abs(ID1) > 4 || abs(ID2) > 4 || abs (ID3) != 21 || abs (ID4) != 21  ) accept = false;
	  }
          if( IDPROD == 532 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only Q QX  --> G G 
            if( abs(ID1) > 4 || abs(ID2) > 4 || ID1 < 0 || abs(ID3) != 21 || abs(ID4) != 21  ) accept = false;
		}
          if( IDPROD == 533 &&  (ID1 * ID2 < 0 ) ) {
	    accept = true;
	    // only QX Q  --> G G 
            if( abs(ID1) > 4 || abs(ID2) > 4 || ID2 < 0 || abs (ID3) != 21 || abs (ID4) != 21  ) accept = false;
	  }
          if( IDPROD == 541 ) {
	    accept = false;
	    // only U SX  --> D CX 
            if( ID1 ==  2 && ID2 == -3 && ID3 == 1 && ID4 == -4 ) accept = true;
            if( ID1 == -3 && ID2 ==  2 && ID3 == 1 && ID4 == -4 ) accept = true;
            if( ID1 ==  2 && ID2 == -3 && ID3 ==-4 && ID4 ==  1 ) accept = true;
            if( ID1 == -3 && ID2 ==  2 && ID3 ==-4 && ID4 ==  1 ) accept = true;
	  }

 	  if(accept ){

	    W[0][0] += f_x1[I1]*f_x2[I2]*vbfdistrWRAP(ID1,ID2,ID3,ID4, 1, -1, P, KEY);    // as in case of nonSM_adopt we change the sign for 
	    W[0][1] += f_x1[I1]*f_x2[I2]*vbfdistrWRAP(ID1,ID2,ID3,ID4, 1,  1, P, KEY);    // second tau helicity assuming without checking that
	    W[1][0] += f_x1[I1]*f_x2[I2]*vbfdistrWRAP(ID1,ID2,ID3,ID4,-1, -1, P, KEY);    // for VBF quantization frame conventions are the same.
	    W[1][1] += f_x1[I1]*f_x2[I2]*vbfdistrWRAP(ID1,ID2,ID3,ID4,-1,  1, P, KEY);

	  }	    
	}
      }
    }
  }
}



void makeSimpleTestME2()
{
 
  //---------------------------------------------------------------------------
  // Construct the matrix for FORTRAN function
  // NOTE: different order of indices than in FORTRAN!

  std::cout << "ERW: ----------------------------- " << std::endl;
  double P1[6][4]  = {{ 0.5000000E+03,  0.0,            0.0,            0.5000000E+03   },
                      { 0.5000000E+03,  0.0,            0.0,           -0.5000000E+03   }, 
                      { 0.8855009E+02, -0.2210038E+02,  0.4007979E+02, -0.7580437E+02   }, 
                      { 0.3283248E+03, -0.1038482E+03, -0.3019295E+03,  0.7649385E+02   },
                      { 0.1523663E+03, -0.1058795E+03, -0.9770827E+02,  0.4954769E+02   },
                      { 0.4307588E+03,  0.2318280E+03,  0.3595580E+03, -0.5023717E+02   } };
  calcTestME2(1, P1);

  Particle p1         (  0.0,            0.0,            0.5000000E+03,  0.5000000E+03,   1  );
  Particle p2         (  0.0,            0.0,           -0.5000000E+03,  0.5000000E+03,  -1  ); 
  Particle p3         ( -0.2210038E+02,  0.4007979E+02, -0.7580437E+02,  0.8855009E+02,  2 ); 
  Particle p4         ( -0.1038482E+03, -0.3019295E+03,  0.7649385E+02,  0.3283248E+03,  -2  );
  Particle tau1       ( -0.1058795E+03, -0.9770827E+02,  0.4954769E+02,  0.1523663E+03,  15 );
  Particle tau2       (  0.2318280E+03,  0.3595580E+03, -0.5023717E+02,  0.4307588E+03,  -15  );

  std::cout << "ERW:  4-momenta before boost ----------------------------- " << std::endl;

  p1.print();
  p2.print();
  p3.print();
  p4.print();
  tau1.print();
  tau2.print();
  double x1=0.03;  // x1 for the first parton
  double xm = 1000.0;
  double cmsene = 13000.0;
  double x2=xm/cmsene*xm/cmsene/x1;    // x2 adjusted from  cms energy (13000), and virtuality of the whole system 1000 
  std::cout << "x1, x2= " << x1 << " " << x2 << std::endl;
  if (x2>=1.0) std::cout << "ERW: wrong x1, x2 cannot be defined (is above 1);  x2=" << x2 << std::endl;

  Particle P_Z0( tau1.px()+tau2.px(), tau1.py()+tau2.py(), tau1.pz()+tau2.pz(), tau1.e()+tau2.e(), 23 );
  P_Z0.print();
  double Q2 = P_Z0.recalculated_mass() *  P_Z0.recalculated_mass();
  std::cout << " pdfs = " << f(x1,0,Q2,cmsene) * f(x2,1,Q2,cmsene) << std::endl;

  Particle P_QQ( 0.0, 1.0e-10, (x1-x2)*cmsene/2.0, (x1+x2)*cmsene/2.0, 10 );
  P_QQ.print();
 
  p1.boostFromRestFrame(P_QQ);
  p2.boostFromRestFrame(P_QQ);
  p3.boostFromRestFrame(P_QQ);
  p4.boostFromRestFrame(P_QQ);
  P_Z0.boostFromRestFrame(P_QQ);
  tau1.boostFromRestFrame(P_QQ);
  tau2.boostFromRestFrame(P_QQ);

  std::cout << "ERW:  4-momenta after boost ----------------------------- " << std::endl;

  p1.print();
  p2.print();
  p3.print();
  p4.print();
  P_Z0.print();
  tau1.print();
  tau2.print();

  std::cout << "ERW: ----------------------------- " << std::endl;
  double P2[6][4]  = {{ 0.5000000E+03,  0.0,            0.0,            0.5000000E+03   },
                      { 0.5000000E+03,  0.0,            0.0,           -0.5000000E+03   }, 
                      { 0.1177462E+03,  -0.6070512E+02,   0.7123011E+02,   0.7145150E+02   }, 
                      { 0.3509495E+03,  -0.3178775E+02,   0.8393832E+02,   0.3392779E+03   },
                      { 0.3493321E+03,   0.1840069E+03,  -0.5152712E+02,  -0.2924315E+03   },
                      { 0.1819722E+03,  -0.9151401E+02,  -0.1036413E+03,  -0.1182978E+03   } };
  calcTestME2(2, P2);
  std::cout << "ERW: ----------------------------- " << std::endl;
  double P3[6][4]  = {{ 0.5000000E+03,  0.0,            0.0,            0.5000000E+03   },
                      { 0.5000000E+03,  0.0,            0.0,           -0.5000000E+03   }, 
                      { 0.2586900E+03,   0.1324670E+03,  -0.1696171E+03,  -0.1435378E+03   }, 
                      { 0.1084567E+03,  -0.5735712E+02,  -0.2162482E+02,  -0.8947281E+02   },
                      { 0.4005742E+03,  -0.1580760E+03,   0.3563160E+03,   0.9223569E+02   },
                      { 0.2322791E+03,   0.8296613E+02,  -0.1650741E+03,   0.1407749E+03   } };
  calcTestME2(3, P3);
  std::cout << "ERW: ----------------------------- " << std::endl;
  double P4[6][4]  = {{ 0.5000000E+03,  0.0,            0.0,            0.5000000E+03   },
                      { 0.5000000E+03,  0.0,            0.0,           -0.5000000E+03   }, 
                      { 0.1595700E+03,  -0.6917808E+02,  -0.1395175E+03,  -0.3481123E+02   }, 
                      { 0.2247758E+03,  -0.1360140E+03,   0.1650340E+03,  -0.6919641E+02   },
                      { 0.2508802E+03,   0.1447863E+01,   0.2499830E+03,  -0.2107335E+02   },
                      { 0.3647740E+03,   0.2037442E+03,  -0.2754995E+03,   0.1250810E+03   } };
  calcTestME2(4, P4);

}

/*******************************************************************************/
  void calcTestME2(int iter, double P1[6][4])
{

  double ME2ref, ME2;
  int ID1, ID2, ID3, ID4;

  int KEY=0;

  // case GD -> GD
  ID1 = 21; ID2 = 1; ID3 = 21; ID4 = 1;

  if(iter == 1) ME2ref =  1.5961579933867344E-010;
  if(iter == 2) ME2ref =  3.8749742050329834E-010;
  if(iter == 3) ME2ref =  5.0434636937545397E-012;
  if(iter == 4) ME2ref =  7.9077184257060427E-012;

  ME2    =   vbfdistrWRAP(ID1,ID2,ID3,ID4,  1, -1, P1,  KEY) 
           + vbfdistrWRAP(ID1,ID2,ID3,ID4,  1,  1, P1,  KEY)  
           + vbfdistrWRAP(ID1,ID2,ID3,ID4, -1,  1, P1,  KEY) 
           + vbfdistrWRAP(ID1,ID2,ID3,ID4, -1, -1, P1,  KEY);  

  std::cout << "ERW: ME2 GD->GD  (ref)     = " << ME2ref << std::endl;
  std::cout << "ERW: ME2 GD->GD            = " << ME2 << "       ratio to ref = " << ME2/ME2ref  << std::endl;

  // case GU -> GU
  ID1 = 21; ID2 = 2; ID3 = 21; ID4 = 2;

  if(iter == 1) ME2ref =  2.9195503763051040E-010;

  ME2    =   vbfdistrWRAP(ID1,ID2,ID3,ID4,  1, -1, P1,  KEY) 
           + vbfdistrWRAP(ID1,ID2,ID3,ID4,  1,  1, P1,  KEY)  
           + vbfdistrWRAP(ID1,ID2,ID3,ID4, -1,  1, P1,  KEY) 
           + vbfdistrWRAP(ID1,ID2,ID3,ID4, -1, -1, P1,  KEY);  

  if(iter == 1)std::cout << "ERW: ME2 GU->GU  (ref)     = " << ME2ref << std::endl;
  if(iter == 1)std::cout << "ERW: ME2 GU->GU            = " << ME2 << "       ratio to ref = " << ME2/ME2ref  << std::endl;

  // case DD -> DD
  ID1 = 1; ID2 = 1; ID3 = 1; ID4 = 1;

  if(iter == 1) ME2ref =  3.3953129762581284E-017;
  if(iter == 2) ME2ref =  6.0054072781075002E-017;
  if(iter == 3) ME2ref =  3.9707833415682912E-018;
  if(iter == 4) ME2ref =  1.5342177192807347E-018;

  ME2    =   vbfdistrWRAP(ID1,ID2,ID3,ID4,  1, -1, P1,  KEY) 
           + vbfdistrWRAP(ID1,ID2,ID3,ID4,  1,  1, P1,  KEY)  
           + vbfdistrWRAP(ID1,ID2,ID3,ID4, -1,  1, P1,  KEY) 
           + vbfdistrWRAP(ID1,ID2,ID3,ID4, -1, -1, P1,  KEY);  

  std::cout << "ERW: ME2 DD->DD  (ref)     = " << ME2ref << std::endl;
  std::cout << "ERW: ME2 DD->DD            = " << ME2 << "       ratio to ref = " << ME2/ME2ref  << std::endl;

  // case UU -> UU
  ID1 = 2; ID2 = 2; ID3 = 2; ID4 = 2;

  if(iter == 1) ME2ref =  1.9412924824120248E-017;
  if(iter == 2) ME2ref =  4.0830132078559964E-017;
  if(iter == 3) ME2ref =  2.1297931556738857E-018;
  if(iter == 4) ME2ref =  7.9215386777281075E-019;

  ME2    =   vbfdistrWRAP(ID1,ID2,ID3,ID4,  1, -1, P1,  KEY) 
           + vbfdistrWRAP(ID1,ID2,ID3,ID4,  1,  1, P1,  KEY)  
           + vbfdistrWRAP(ID1,ID2,ID3,ID4, -1,  1, P1,  KEY) 
           + vbfdistrWRAP(ID1,ID2,ID3,ID4, -1, -1, P1,  KEY);  

  std::cout << "ERW: ME2 UU->UU  (ref)     = " << ME2ref << std::endl;
  std::cout << "ERW: ME2 UU->UU            = " << ME2 << "       ratio to ref = " << ME2/ME2ref  << std::endl;

  // case DDX -> CCX
  ID1 = 1; ID2 = -1; ID3 = 4; ID4 = -4;

  if(iter == 1) ME2ref =  3.6780119739685137E-020;
  if(iter == 2) ME2ref =  5.2280923900274694E-018;
  if(iter == 3) ME2ref =  2.9001589049209953E-019;
  if(iter == 4) ME2ref =  5.8509026904432882E-020;

  ME2    =   vbfdistrWRAP(ID1,ID2,ID3,ID4,  1, -1, P1,  KEY) 
           + vbfdistrWRAP(ID1,ID2,ID3,ID4,  1,  1, P1,  KEY)  
           + vbfdistrWRAP(ID1,ID2,ID3,ID4, -1,  1, P1,  KEY) 
           + vbfdistrWRAP(ID1,ID2,ID3,ID4, -1, -1, P1,  KEY);  

  std::cout << "ERW: ME2 DDX->CCX  (ref)   = " << ME2ref << std::endl;
  std::cout << "ERW: ME2 DDX->CCX          = " << ME2 << "       ratio to ref = " << ME2/ME2ref  << std::endl;

}

  //
  // -------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------
  //

void makeSimpleTestPDF()
{
 
  std::cout << "------------------------------------" << std::endl;
  std::cout << "ERW: TEST PDF only = " << std::endl;
  std::cout << "ERW: PDF:LHAPDFset = cteq6ll.LHpdf " << std::endl;
  std::cout << "ERW: reference from PYTHIA printout " << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "ERW: ID1,ID2 =2 -2 " << std::endl;
  std::cout << "ERW: Q,Q2 = 97.868, 9578.118 " << std::endl;
  double pdf1= f(2.51297e-01,2,9578.118,97.868);
  std::cout << "ERW: x1, xpdf1= "  <<  2.51297e-01  << "  " << 2.51297e-01  *pdf1  << "  ref pdf = " <<  0.411  << std::endl;
  double pdf2= f(1.94463e-04,-2,9578.118,97.868);
  std::cout << "ERW: x2, xpdf2= " <<  1.94463e-04 << "  " << 1.94463e-04  *pdf2  <<  "  ref pdf = " <<  2.989  << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "ERW: ID1,ID2 =2 -2 " << std::endl;
  std::cout << "ERW: Q,Q2 = 38.747, 1501.314 " << std::endl;
  pdf1= f(4.17758e-04,-2,1501.314,38.747);
  std::cout << "ERW: x1, xpdf1 = "  <<  4.17758e-04  << "  " <<  4.17758e-04 *pdf1  << "  ref pdf = " <<  1.768 << std::endl;
  pdf2= f(1.83354e-02, 2,1501.314,38.747);
  std::cout << "ERW: x2, xpdf2 = " <<  1.83354e-02 << "  " << 1.83354e-02  *pdf2  <<  "  ref pdf = " <<  0.610 << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "ERW: ID1,ID2 =1 -1 " << std::endl;
  std::cout << "ERW: Q,Q2 =  91.585, 8387.831 " << std::endl;
  pdf1= f(6.29049e-02, 1,8387.831,91.585);
  std::cout << "ERW: x1, xpdf1 = "  <<  6.29049e-02  << "  " <<  6.29049e-02 *pdf1 << "  ref pdf = " <<  0.386  << std::endl;
  pdf2= f(6.80314e-04, -1,8387.831,91.585);
  std::cout << "ERW: x2, xpdf2 = "  <<  6.80314e-04 << "  " << 6.80314e-04  *pdf2  << "  ref pdf = " <<  1.751 << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "ERW: ID1,ID2 =-2 2 " << std::endl;
  std::cout << "ERW: Q,Q2 =  12.979, 168.453 " << std::endl;
  pdf1= f( 1.29739e-01,  -2,168.453,12.979);
  std::cout << "ERW: x1, xpdf1 = "  <<  1.29739e-01  << "  " << 1.29739e-01*pdf1 << "  ref pdf = " <<  0.066  << std::endl;
  pdf2= f( 6.62449e-06,  2,168.453,12.979);
  std::cout << "ERW: x2, xpdf2 = "  <<  6.62449e-06 << "  " << 6.62449e-06  *pdf2  << "  ref pdf = " << 4.559  << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "ERW: ID1,ID2 = 2 -2 " << std::endl;
  std::cout << "ERW: Q,Q2 =  91.707, 8410.127 " << std::endl;
  pdf1= f(1.73131e-02 ,  2,8410.127,91.707);
  std::cout << "ERW: x1, xpdf1 = "  <<  1.73131e-02  << "  " << 1.73131e-02 *pdf1 << "  ref pdf = " <<   0.643  << std::endl;
  pdf2= f( 2.47840e-03,  -2,8410.127,91.707);
  std::cout << "ERW: x2, xpdf2 = "  <<  2.47840e-03 << "  " << 2.47840e-03  *pdf2  << "  ref pdf = " << 0.986  << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "ERW: ID1,ID2 = 1 -1 " << std::endl;
  std::cout << "ERW: Q,Q2 =  91.203 8317.936 " << std::endl;
  pdf1= f( 3.87999e-04 , 1, 8317.936, 91.203);
  std::cout << "ERW: x1, xpdf1 = "  <<  3.87999e-04  << "  " <<  3.87999e-04 *pdf1 << "  ref pdf = " <<   2.256  << std::endl;
  pdf2= f(  1.09378e-01, -1,8317.936, 91.203);
  std::cout << "ERW: x2, xpdf2 = "  <<  1.09378e-01  << "  " <<   1.09378e-01 *pdf2  << "  ref pdf = " << 0.106  << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "ERW: ID1,ID2 = 1 0 " << std::endl;
  std::cout << "ERW: Q,Q2 =  51.963, 2700.201 " << std::endl;
  pdf1= f( 9.22157e-02 , 1, 2700.201, 51.963);
  std::cout << "ERW: x1, xpdf1 = "  <<  9.22157e-02  << "  " <<  9.22157e-02 *pdf1 << "  ref pdf = " << 0.353   << std::endl;
  pdf2= f( 2.50786e-03 , 0,2700.201, 51.963);
  std::cout << "ERW: x2, xpdf2 = "  <<  2.50786e-03  << "  " <<  2.50786e-03 *pdf2  << "  ref pdf = " <<  20.362  << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "ERW: ID1,ID2 = 1 0 " << std::endl;
  std::cout << "ERW: Q,Q2 =  55.719, 3104.648  " << std::endl;
  pdf1= f( 1.28045e-01 , 1,  3104.648, 55.719);
  std::cout << "ERW: x1, xpdf1 = "  <<   1.28045e-01 << "  " <<  1.28045e-01 *pdf1 << "  ref pdf = " << 0.312   << std::endl;
  pdf2= f( 3.12973e-03 , 0, 3104.648, 55.719);
  std::cout << "ERW: x2, xpdf2 = "  <<   3.12973e-03 << "  " <<  3.12973e-03 *pdf2  << "  ref pdf = " <<  17.938  << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "ERW: ID1,ID2 = 3 0 " << std::endl;
  std::cout << "ERW: Q,Q2 =  106.126, 11262.651  " << std::endl;
  pdf1= f(  8.02832e-02 , 3,  11262.651, 106.126);
  std::cout << "ERW: x1, xpdf1 = "  <<   8.02832e-02 << "  " <<   8.02832e-02  *pdf1 << "  ref pdf = " << 0.077   << std::endl;
  pdf2= f(  3.85011e-03 , 0, 11262.651, 106.126);
  std::cout << "ERW: x2, xpdf2 = "  <<    3.85011e-03 << "  " <<   3.85011e-03 *pdf2  << "  ref pdf = " <<  16.542  << std::endl;
  std::cout << "------------------------------------" << std::endl;

  std::cout << "ERW: Q= 98.069051344553060, x =1.36851248162658170E-002 " << std::endl;
  double x = 1.36851248162658170E-002;
  std::cout << "ERW: f(-5) = " <<  f(  1.36851248162658170E-002 , -5, 9617.53, 98.069 ) << " ref=12.649157008859417" << std::endl;
  std::cout << "ERW: f(-4) = " <<  f(  1.36851248162658170E-002 , -4, 9617.53, 98.069 ) << " ref=19.196982264477899" << std::endl;
  std::cout << "ERW: f(-3) = " <<  f(  1.36851248162658170E-002 , -3, 9617.53, 98.069 ) << " ref=23.971647011421041" << std::endl;
  std::cout << "ERW: f(-2) = " <<  f(  1.36851248162658170E-002 , -2, 9617.53, 98.069 ) << " ref=30.602439987620933" << std::endl;
  std::cout << "ERW: f(-1) = " <<  f(  1.36851248162658170E-002 , -1, 9617.53, 98.069 ) << " ref=31.664848276050574" << std::endl;
  std::cout << "ERW: f(0)  = " <<  f(  1.36851248162658170E-002 ,  0, 9617.53, 98.069 ) << " ref=483.13004227606962" << std::endl;
  std::cout << "ERW: f( 1) = " <<  f(  1.36851248162658170E-002 ,  1, 9617.53, 98.069 ) << " ref=41.868651977841559" << std::endl;
  std::cout << "ERW: f( 2) = " <<  f(  1.36851248162658170E-002 ,  2, 9617.53, 98.069 ) << " ref=49.256024133250129" << std::endl;
  std::cout << "ERW: f( 3) = " <<  f(  1.36851248162658170E-002 ,  3, 9617.53, 98.069 ) << " ref=23.971647011421041" << std::endl;
  std::cout << "ERW: f( 4) = " <<  f(  1.36851248162658170E-002 ,  4, 9617.53, 98.069 ) << " ref=19.196982264477899" << std::endl;
  std::cout << "ERW: f( 5) = " <<  f(  1.36851248162658170E-002 ,  5, 9617.53, 98.069 ) << " ref=12.649157008859417" << std::endl;

  std::cout << "ERW: f(-5) = " <<  f(  1.36851248162658170E-002 , -5, 98.069, 98.069 ) << " ref=12.649157008859417" << std::endl;
  std::cout << "ERW: f(-4) = " <<  f(  1.36851248162658170E-002 , -4, 98.069, 98.069 ) << " ref=19.196982264477899" << std::endl;
  std::cout << "ERW: f(-3) = " <<  f(  1.36851248162658170E-002 , -3, 98.069, 98.069 ) << " ref=23.971647011421041" << std::endl;
  std::cout << "ERW: f(-2) = " <<  f(  1.36851248162658170E-002 , -2, 98.069, 98.069 ) << " ref=30.602439987620933" << std::endl;
  std::cout << "ERW: f(-1) = " <<  f(  1.36851248162658170E-002 , -1, 98.069, 98.069 ) << " ref=31.664848276050574" << std::endl;
  std::cout << "ERW: f(0)  = " <<  f(  1.36851248162658170E-002 ,  0, 98.069, 98.069 ) << " ref=483.13004227606962" << std::endl;
  std::cout << "ERW: f( 1) = " <<  f(  1.36851248162658170E-002 ,  1, 98.069, 98.069 ) << " ref=41.868651977841559" << std::endl;
  std::cout << "ERW: f( 2) = " <<  f(  1.36851248162658170E-002 ,  2, 98.069, 98.069 ) << " ref=49.256024133250129" << std::endl;
  std::cout << "ERW: f( 3) = " <<  f(  1.36851248162658170E-002 ,  3, 98.069, 98.069 ) << " ref=23.971647011421041" << std::endl;
  std::cout << "ERW: f( 4) = " <<  f(  1.36851248162658170E-002 ,  4, 98.069, 98.069 ) << " ref=19.196982264477899" << std::endl;
  std::cout << "ERW: f( 5) = " <<  f(  1.36851248162658170E-002 ,  5, 98.069, 98.069 ) << " ref=12.649157008859417" << std::endl;


}


//---------------------------------------------------------------------------------------------------------------------------------------------------
void calcSumME2(SimpleParticle &p3, SimpleParticle &p4, SimpleParticle &sp_X,SimpleParticle &tau1, SimpleParticle &tau2, double (&W)[2][2], int KEY,
                    int ID1, int ID2, int ID3, int ID4)
{
  /*
  // this may be necessary because of matrix element calculation may require absolute energy-momentum conservation!
  // FSR photons may need to be treated explicitely or with interpolation procedures.
  Particle P_QQ( p3.px()+p4.px()+tau1.px()+tau2.px(),
                 p3.py()+p4.py()+tau1.py()+tau2.py(),
                 p3.pz()+p4.pz()+tau1.pz()+tau2.pz(),
                 p3.e() +p4.e() +tau1.e() +tau2.e(), 0 );
  */
  Particle P_QQ( p3.px()+p4.px()+sp_X.px(),
                 p3.py()+p4.py()+sp_X.py(),
                 p3.pz()+p4.pz()+sp_X.pz(),
                 p3.e() +p4.e() +sp_X.e() , 0 );


  double SS = P_QQ.recalculated_mass()*P_QQ.recalculated_mass(); 
  
  double x1x2  = SS/CMSENE/CMSENE;
  double x1Mx2 = P_QQ.pz()/CMSENE*2;
  
  double x1 = (  x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;
  double x2 = ( -x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;

  //---------------------------------------------------------------------------
  // Construct the matrix for FORTRAN function
  // NOTE: different order of indices than in FORTRAN!
  double P[6][4] = { { CMSENE/2*x1, 0.0,       0.0,          CMSENE/2*x1  },
                     { CMSENE/2*x2, 0.0,       0.0,         -CMSENE/2*x2  }, 
                     { p3.e(),      p3.px(),   p3.py(),      p3.pz()      }, 
                     { p4.e(),      p4.px(),   p4.py(),      p4.pz()      },
                     { tau1.e(),    tau1.px(), tau1.py(),    tau1.pz()    },
                     { tau2.e(),    tau2.px(), tau2.py(),    tau2.pz()    } };
 

  W[0][0] = vbfdistrWRAP(ID1,ID2,ID3,ID4, 1, -1, P, KEY);   
  W[0][1] = vbfdistrWRAP(ID1,ID2,ID3,ID4, 1,  1, P, KEY);  
  W[1][0] = vbfdistrWRAP(ID1,ID2,ID3,ID4,-1, -1, P, KEY); 
  W[1][1] = vbfdistrWRAP(ID1,ID2,ID3,ID4,-1,  1, P, KEY);

  double  sumW =(W[0][0]+W[0][1]+ W[1][0]+W[1][1]); 
  if( sumW < 0 ){
    std::cout << "ERW: ==== " << std::endl; 
    std::cout << "ERW: W[]  = " << W[0][0] << " " << W[0][1] << " " << W[1][0] << " " << W[1][1]	<< std::endl; 
    std::cout << "ERW: ==== " << std::endl; 
    std::cout << "ERW: p1   = " << CMSENE/2*x1 << "  " << 0.0 << "  " << 0.0 << "  " <<   CMSENE/2*x1 	<< std::endl; 
    std::cout << "ERW: p2   = " << CMSENE/2*x2 << "  " << 0.0 << "  " << 0.0 << "  " <<  -CMSENE/2*x2 	<< std::endl; 
    std::cout << "ERW:  X   = " << sp_X.e()  << "  " << sp_X.px()  << "  " << sp_X.py()  << "  " <<  sp_X.pz() 	<< std::endl; 
    std::cout << "ERW: p3   = " << p3.e() << "  " << p3.px() << "  " << p3.py() << "  " <<  p3.pz() 	<< std::endl; 
    std::cout << "ERW: p4   = " << p4.e() << "  " << p4.px() << "  " << p4.py() << "  " <<  p4.pz() 	<< std::endl; 
    std::cout << "ERW: tau1 = " << tau1.e() << "  " << tau1.px() << "  " << tau1.py() << "  " <<  tau1.pz() 	<< std::endl; 
    std::cout << "ERW: tau2 = " << tau2.e() << "  " << tau2.px() << "  " << tau2.py() << "  " <<  tau2.pz() 	<< std::endl; 
    std::cout << "ERW: ==== " << std::endl; 
    std::cout << "ERW:   p1+p2= " << P[0][0]+P[1][0] <<  "  " <<  P[0][1]+P[1][1] << "  "  <<  P[0][2]+P[1][2] <<  "  "  <<  P[0][3]+P[1][3]	<< std::endl;
    std::cout << "ERW: X+p3+p4= " << P[2][0]+P[3][0]+P[4][0]+P[5][0] <<  "  " <<  P[2][1]+P[3][1]+P[4][1]+P[5][1] << "  " 
                                  << P[2][2]+P[3][2]+P[4][2]+P[5][2] <<  "  " <<  P[2][3]+P[3][3]+P[4][3]+P[5][3] << std::endl;

    double p3sqr = p3.e()*p3.e()-p3.px()*p3.px()-p3.py()*p3.py()-p3.pz()*p3.pz();
    double p4sqr = p4.e()*p4.e()-p4.px()*p4.px()-p4.py()*p4.py()-p4.pz()*p4.pz();
    std::cout << "ERW: p3^2, p4^2 =" <<  p3sqr  << "  " <<  p4sqr << std::endl;
  }

}

//---------------------------------------------------------------------------------------------------------------------------------------------------
void calcPDFs(SimpleParticle &p3, SimpleParticle &p4, SimpleParticle &sp_X,SimpleParticle &tau1, SimpleParticle &tau2, double (&W)[2][2], int KEY,
                    int ID1, int ID2, int ID3, int ID4, int pdfOpt)
{
  /*
  // this may be necessary because of matrix element calculation may require absolute energy-momentum conservation!
  // FSR photons may need to be treated explicitely or with interpolation procedures.
  Particle P_QQ( p3.px()+p4.px()+tau1.px()+tau2.px(),
                 p3.py()+p4.py()+tau1.py()+tau2.py(),
                 p3.pz()+p4.pz()+tau1.pz()+tau2.pz(),
                 p3.e() +p4.e() +tau1.e() +tau2.e(), 0 );
  */
  Particle P_QQ( p3.px()+p4.px()+sp_X.px(),
                 p3.py()+p4.py()+sp_X.py(),
                 p3.pz()+p4.pz()+sp_X.pz(),
                 p3.e() +p4.e() +sp_X.e() , 0 );


  double SS = P_QQ.recalculated_mass()*P_QQ.recalculated_mass(); 
  
  double x1x2  = SS/CMSENE/CMSENE;
  double x1Mx2 = P_QQ.pz()/CMSENE*2;
  
  double x1 = (  x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;
  double x2 = ( -x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;
  
  //---------------------------------------------------------------------------
  // Construct the matrix for FORTRAN function
  // NOTE: different order of indices than in FORTRAN!
  double P[6][4] = { { CMSENE/2*x1, 0.0,       0.0,          CMSENE/2*x1  },
                     { CMSENE/2*x2, 0.0,       0.0,         -CMSENE/2*x2  }, 
                     { p3.e(),      p3.px(),   p3.py(),      p3.pz()      }, 
                     { p4.e(),      p4.px(),   p4.py(),      p4.pz()      },
                     { tau1.e(),    tau1.px(), tau1.py(),    tau1.pz()    },
                     { tau2.e(),    tau2.px(), tau2.py(),    tau2.pz()    } };
 

  //
  // Calculate 'f' function for all x1 and all ID1, ID2
  //
  double  f_x1_ARRAY[11] = { 0.0 };
  double  f_x2_ARRAY[11] = { 0.0 };
  double *f_x1 = f_x1_ARRAY+5;     // This way f_x1[i],f_x2[i] can be used with 
  double *f_x2 = f_x2_ARRAY+5;     // i going from -5 to 5


  Particle P_X(  sp_X.px(), sp_X.py(), sp_X.pz(), sp_X.e() , 0 );


  double Q2  =  P_X.recalculated_mass()*P_X.recalculated_mass();
  double PT2 =  sp_X.px()* sp_X.px() +  sp_X.py()*  sp_X.py(); 

  if(pdfOpt == 0){
    for(int i=-5;i<=5;++i) {
      f_x1[i] = f(x1,i,Q2,CMSENE);
      f_x2[i] = f(x2,i,Q2,CMSENE);
    }
  } else if (pdfOpt == 1) {
    for(int i=-5;i<=5;++i) {
      f_x1[i] = f(x1,i,PT2,CMSENE);
      f_x2[i] = f(x2,i,PT2,CMSENE);
    }
  } else if (pdfOpt == 2) {
    double PT24 = PT2/4.;
    for(int i=-5;i<=5;++i) {
      f_x1[i] = f(x1,i,PT24,CMSENE);
      f_x2[i] = f(x2,i,PT24,CMSENE);
    }
  }


  int I1=ID1;
  int I2=ID2;
  if( I1 == 21) I1=0;
  if( I2 == 21) I2=0;

  W[0][0] = f_x1[I1]*f_x2[I2];   
  W[0][1] = f_x1[I1]*f_x2[I2];  
  W[1][0] = f_x1[I1]*f_x2[I2]; 
  W[1][1] = f_x1[I1]*f_x2[I2];


}


//---------------------------------------------------------------------------------------------------------------------------------------------------
void calcProdMatrix(SimpleParticle &p3, SimpleParticle &p4, SimpleParticle &sp_X,SimpleParticle &tau1, SimpleParticle &tau2, double (&W)[2][2], int KEY,
                    int ID1, int ID2, int ID3, int ID4, int pdfOpt)
{
  /*
  // this may be necessary because of matrix element calculation may require absolute energy-momentum conservation!
  // FSR photons may need to be treated explicitely or with interpolation procedures.
  Particle P_QQ( p3.px()+p4.px()+tau1.px()+tau2.px(),
                 p3.py()+p4.py()+tau1.py()+tau2.py(),
                 p3.pz()+p4.pz()+tau1.pz()+tau2.pz(),
                 p3.e() +p4.e() +tau1.e() +tau2.e(), 0 );
  */
  Particle P_QQ( p3.px()+p4.px()+sp_X.px(),
                 p3.py()+p4.py()+sp_X.py(),
                 p3.pz()+p4.pz()+sp_X.pz(),
                 p3.e() +p4.e() +sp_X.e() , 0 );


  double SS = P_QQ.recalculated_mass()*P_QQ.recalculated_mass(); 
  
  double x1x2  = SS/CMSENE/CMSENE;
  double x1Mx2 = P_QQ.pz()/CMSENE*2;
  
  double x1 = (  x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;
  double x2 = ( -x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;
  
  //---------------------------------------------------------------------------
  // Construct the matrix for FORTRAN function
  // NOTE: different order of indices than in FORTRAN!
  double P[6][4] = { { CMSENE/2*x1, 0.0,       0.0,          CMSENE/2*x1  },
                     { CMSENE/2*x2, 0.0,       0.0,         -CMSENE/2*x2  }, 
                     { p3.e(),      p3.px(),   p3.py(),      p3.pz()      }, 
                     { p4.e(),      p4.px(),   p4.py(),      p4.pz()      },
                     { tau1.e(),    tau1.px(), tau1.py(),    tau1.pz()    },
                     { tau2.e(),    tau2.px(), tau2.py(),    tau2.pz()    } };
 

  //
  // Calculate 'f' function for all x1 and all ID1, ID2
  //
  double  f_x1_ARRAY[11] = { 0.0 };
  double  f_x2_ARRAY[11] = { 0.0 };
  double *f_x1 = f_x1_ARRAY+5;     // This way f_x1[i],f_x2[i] can be used with 
  double *f_x2 = f_x2_ARRAY+5;     // i going from -5 to 5


  Particle P_X(  sp_X.px(), sp_X.py(), sp_X.pz(), sp_X.e() , 0 );


  double Q2  =  P_X.recalculated_mass()*P_X.recalculated_mass();
  double PT2 =  sp_X.px()* sp_X.px() +  sp_X.py()*  sp_X.py(); 

  if(pdfOpt == 0){
    for(int i=-5;i<=5;++i) {
      f_x1[i] = f(x1,i,Q2,CMSENE);
      f_x2[i] = f(x2,i,Q2,CMSENE);
    }
  } else if (pdfOpt == 1) {
    for(int i=-5;i<=5;++i) {
      f_x1[i] = f(x1,i,PT2,CMSENE);
      f_x2[i] = f(x2,i,PT2,CMSENE);
    }
  } else if (pdfOpt == 2) {
    double PT24 = PT2/4.;
    for(int i=-5;i<=5;++i) {
      f_x1[i] = f(x1,i,PT24,CMSENE);
      f_x2[i] = f(x2,i,PT24,CMSENE);
    }
  }


  int I1=ID1;
  int I2=ID2;
  if( I1 == 21) I1=0;
  if( I2 == 21) I2=0;

  W[0][0] = f_x1[I1]*f_x2[I2]*vbfdistrWRAP(ID1,ID2,ID3,ID4, 1, -1, P, KEY);   
  W[0][1] = f_x1[I1]*f_x2[I2]*vbfdistrWRAP(ID1,ID2,ID3,ID4, 1,  1, P, KEY);  
  W[1][0] = f_x1[I1]*f_x2[I2]*vbfdistrWRAP(ID1,ID2,ID3,ID4,-1, -1, P, KEY); 
  W[1][1] = f_x1[I1]*f_x2[I2]*vbfdistrWRAP(ID1,ID2,ID3,ID4,-1,  1, P, KEY);

  double  sumW =(W[0][0]+W[0][1]+ W[1][0]+W[1][1]); 
  if( sumW < 0 ){
    std::cout << "ERW: f1*f2= " << f_x1[I1]*f_x2[I2]  << "  x1=" << x1 << " x2=" << x2 << " f1=" << f_x1[I1] << " f2=" << f_x2[I1]	<< std::endl; 
    std::cout << "ERW: ==== " << std::endl; 
    std::cout << "ERW: W[]  = " << W[0][0] << " " << W[0][1] << " " << W[1][0] << " " << W[1][1]	<< std::endl; 
    std::cout << "ERW: ==== " << std::endl; 
    std::cout << "ERW: p1   = " << CMSENE/2*x1 << "  " << 0.0 << "  " << 0.0 << "  " <<   CMSENE/2*x1 	<< std::endl; 
    std::cout << "ERW: p2   = " << CMSENE/2*x2 << "  " << 0.0 << "  " << 0.0 << "  " <<  -CMSENE/2*x2 	<< std::endl; 
    std::cout << "ERW:  X   = " << sp_X.e()  << "  " << sp_X.px()  << "  " << sp_X.py()  << "  " <<  sp_X.pz() 	<< std::endl; 
    std::cout << "ERW: p3   = " << p3.e() << "  " << p3.px() << "  " << p3.py() << "  " <<  p3.pz() 	<< std::endl; 
    std::cout << "ERW: p4   = " << p4.e() << "  " << p4.px() << "  " << p4.py() << "  " <<  p4.pz() 	<< std::endl; 
    std::cout << "ERW: tau1 = " << tau1.e() << "  " << tau1.px() << "  " << tau1.py() << "  " <<  tau1.pz() 	<< std::endl; 
    std::cout << "ERW: tau2 = " << tau2.e() << "  " << tau2.px() << "  " << tau2.py() << "  " <<  tau2.pz() 	<< std::endl; 
    std::cout << "ERW: ==== " << std::endl; 
    std::cout << "ERW:   p1+p2= " << P[0][0]+P[1][0] <<  "  " <<  P[0][1]+P[1][1] << "  "  <<  P[0][2]+P[1][2] <<  "  "  <<  P[0][3]+P[1][3]	<< std::endl;
    std::cout << "ERW: X+p3+p4= " << P[2][0]+P[3][0]+P[4][0]+P[5][0] <<  "  " <<  P[2][1]+P[3][1]+P[4][1]+P[5][1] << "  " 
                                  << P[2][2]+P[3][2]+P[4][2]+P[5][2] <<  "  " <<  P[2][3]+P[3][3]+P[4][3]+P[5][3] << std::endl;

    double p3sqr = p3.e()*p3.e()-p3.px()*p3.px()-p3.py()*p3.py()-p3.pz()*p3.pz();
    double p4sqr = p4.e()*p4.e()-p4.px()*p4.px()-p4.py()*p4.py()-p4.pz()*p4.pz();
    std::cout << "ERW: p3^2, p4^2 =" <<  p3sqr  << "  " <<  p4sqr << std::endl;
  }

}


/*******************************************************************************
    Calculate weights, case of event record vertex like Z/gamma/H ... -> tau tau decay plus 2 jets.

  Determine taus decay channel, calculates all necessary components from decay HHp HHm and production W[][] for 
  calculation of all weights, later calculates weights.

  Input:        X four momentum of Z/Higgs may be larger than sum of tau1 tau2, missing component
                is assumed to be QED brem four momenta of outgoing partons and taus; vectors of decay products for first and second tau

  Hidden input:  relWTnonSM,  switch of what is  WTnonSM
  Hidden output: WTamplitM or WTamplitP matrix elements of tau1 tau2 decays


                    Polari - helicity attributed to taus, 100% correlations between tau+ and tau-
                    accesible with getTauSpin()
                   WTnonSM weight for introduction of matrix element due to nonSM in cross section, or just ME^2*PDF's for  relWTnonS==false 
  Explicit output: WT spin correlation weight 
*******************************************************************************/
double calculateWeightFromParticlesVBFPROD(int IDPROD, SimpleParticle &p3, SimpleParticle &p4,SimpleParticle &sp_X, SimpleParticle &sp_tau1, SimpleParticle &sp_tau2, vector<SimpleParticle> &sp_tau1_daughters, vector<SimpleParticle> &sp_tau2_daughters, int KEY)
{
  SimpleParticle         sp_tau;
  SimpleParticle         sp_nu_tau;
  vector<SimpleParticle> sp_tau_daughters;
  
  // First iteration is for tau plus, so the 'nu_tau' is tau minus
  if (sp_tau1.pdgid() == -15 )
  {
    sp_tau           = sp_tau1;
    sp_nu_tau        = sp_tau2;
    sp_tau_daughters = sp_tau1_daughters;
  }
  else
  {
    sp_tau           = sp_tau2;
    sp_nu_tau        = sp_tau1;
    sp_tau_daughters = sp_tau2_daughters;
  }

  double *HHp, *HHm;
  
  // We use this to separate namespace for tau+ and tau-
  if(true)
  {
    // Create Particles from SimpleParticles
    Particle X     (      sp_X.px(),      sp_X.py(),      sp_X.pz(),      sp_X.e(),      sp_X.pdgid() );
    Particle tau   (    sp_tau.px(),    sp_tau.py(),    sp_tau.pz(),    sp_tau.e(),    sp_tau.pdgid() );
    Particle nu_tau( sp_nu_tau.px(), sp_nu_tau.py(), sp_nu_tau.pz(), sp_nu_tau.e(), sp_nu_tau.pdgid() );

    vector<Particle> tau_daughters;

    // tau pdgid
    int tau_pdgid = sp_tau.pdgid();

    // Create list of tau daughters
    for(unsigned int i=0; i<sp_tau_daughters.size(); i++)
    {
      Particle pp(sp_tau_daughters[i].px(),
                  sp_tau_daughters[i].py(),
                  sp_tau_daughters[i].pz(),
                  sp_tau_daughters[i].e(),
                  sp_tau_daughters[i].pdgid() );

      tau_daughters.push_back(pp);
    }

    double phi2 = 0.0, theta2 = 0.0;


    //  Move decay kinematics first to tau rest frame  with z axis pointing along nu_tau direction
    //  later rotate again to have neutrino from tau decay along z axis: angles phi2, theta2
    prepareKinematicForHH   (tau, nu_tau, tau_daughters, &phi2, &theta2);


    //  Determine decay channel and then calculate polarimetric vector HH
    HHp = calculateHH(tau_pdgid, tau_daughters, phi2, theta2);

    DEBUG
    (
      cout<<tau_pdgid<<" -> ";
      for(unsigned int i=0;i<tau_daughters.size();i++) cout<<tau_daughters[i].pdgid()<<" ";
      cout<<" (HHp: "<<HHp[0]<<" "<<HHp[1]<<" "<<HHp[2]<<" "<<HHp[3]<<") ";
      cout<<endl;
    )

    WTamplitP = WTamplit;
  } // end of tau+

  // Second iteration is for tau minus, so the 'nu_tau' is tau minus
  if(sp_tau1.pdgid() == 15 )
  {
    sp_tau           = sp_tau1;
    sp_nu_tau        = sp_tau2;
    sp_tau_daughters = sp_tau1_daughters;
  }
  else
  {
    sp_tau           = sp_tau2;
    sp_nu_tau        = sp_tau1;
    sp_tau_daughters = sp_tau2_daughters;
  }
  
  // We use this to separate namespace for tau+ and tau-
  if(true)
  {
    // Create Particles from SimpleParticles
    Particle X     (      sp_X.px(),      sp_X.py(),      sp_X.pz(),      sp_X.e(),      sp_X.pdgid() );
    Particle tau   (    sp_tau.px(),    sp_tau.py(),    sp_tau.pz(),    sp_tau.e(),    sp_tau.pdgid() );
    Particle nu_tau( sp_nu_tau.px(), sp_nu_tau.py(), sp_nu_tau.pz(), sp_nu_tau.e(), sp_nu_tau.pdgid() );

    vector<Particle> tau_daughters;

    // tau pdgid
    int tau_pdgid = sp_tau.pdgid();

    // Create list of tau daughters
    for(unsigned int i=0; i<sp_tau_daughters.size(); i++)
    {
      Particle pp(sp_tau_daughters[i].px(),
                  sp_tau_daughters[i].py(),
                  sp_tau_daughters[i].pz(),
                  sp_tau_daughters[i].e(),
                  sp_tau_daughters[i].pdgid() );

      tau_daughters.push_back(pp);
    }

    double phi2 = 0.0, theta2 = 0.0;


    //  Move decay kinematics first to tau rest frame  with z axis pointing along nu_tau direction
    //  later rotate again to have neutrino from tau decay along z axis: angles phi2, theta2
    prepareKinematicForHH   (tau, nu_tau, tau_daughters, &phi2, &theta2);


    //  Determine decay channel and then calculate polarimetric vector HH
    HHm = calculateHH(tau_pdgid, tau_daughters, phi2, theta2);

    DEBUG
    (
      cout<<tau_pdgid<<" -> ";
      for(unsigned int i=0;i<tau_daughters.size();i++) cout<<tau_daughters[i].pdgid()<<" ";
      cout<<" (HHm: "<<HHm[0]<<" "<<HHm[1]<<" "<<HHm[2]<<" "<<HHm[3]<<") ";
      cout<<endl;
    )

    WTamplitM = WTamplit; 
  } // end of tau-


  double W[2][2] = { { 0.25, 0.25 },
                     { 0.25, 0.25 } };  // this  is trivial W spin (helicity only) density matrix
                                        // consist  of  ME^2*PDF's for production.
                                        // pre-initialization all equal i.e. no spin effects
  
 
  // calculate ME^2*PDF's for SM
  if(KEY==0 || KEY==1 )calcXsect(IDPROD, p3, p4, sp_X, sp_tau1, sp_tau2, W, KEY); 
  else                 calcXsect(IDPROD, p3, p4, sp_X, sp_tau1, sp_tau2, W, 0);


  double sum=(W[0][0]+W[0][1]+ W[1][0]+W[1][1]); // getVBFspins calculated PDF's times matrix elements squared for each of four helicity configurations  
                                                 // tau+ and tau- using SM process KEY=0 DY, KEY=1 Higgs. 

  if(KEY==0 || KEY==1 ) { 
    WTnonSM=1.0;
    if(relWTnonSM==0)  WTnonSM=sum;  // The variable accessible for user stores production weight ME^2 * PDF's summed over flavours and spin states
  }
  if(KEY>=2) { 
    // now we re-calculate W  using anomalous variant of matrix elements.
    calcXsect(IDPROD, p3, p4, sp_X, sp_tau1, sp_tau2, W, KEY);

    double sum2=(W[0][0]+W[0][1]+ W[1][0]+W[1][1]);

    WTnonSM=sum2/sum;
    if(relWTnonSM==0)  WTnonSM=sum2;
  }

  double WT = W[0][0]*(1+HHp[2])*(1+HHm[2])+W[0][1]*(1+HHp[2])*(1-HHm[2])+ W[1][0]*(1-HHp[2])*(1+HHm[2])+W[1][1]*(1-HHp[2])*(1-HHm[2]);
  WT = WT/(W[0][0]+W[0][1]+ W[1][0]+W[1][1]);

  // we separate cross section into first tau helicity + and - parts. Second tau follow.
  // This must be tested  especially for  KEY >=2
  // case for scalar (H) has flipped sign of second tau spin? Tests needed
   double RRR = Tauola::randomDouble();  
  Polari=1.0;
  if (RRR<(W[0][0]*(1+HHp[2])*(1+HHm[2])+W[0][1]*(1+HHp[2])*(1-HHm[2]))/(W[0][0]*(1+HHp[2])*(1+HHm[2])+W[0][1]*(1+HHp[2])*(1-HHm[2])+W[1][0]*(1-HHp[2])*(1+HHm[2])+W[1][1]*(1-HHp[2])*(1-HHm[2]))) Polari=-1.0;



  // Print out some info about the channel
  DEBUG( cout<<" WT: "<<WT<<endl; )

  if (WT<0.0) {
    printf("WT is: %13.10f. Setting WT = 0.0\n",WT);
    WT = 0.0;
   }

  if (WT>4.0) {
    printf("WT is: %13.10f. Setting WT = 4.0\n",WT);
    WT = 4.0;
  }

  delete[] HHp;
  delete[] HHm;
  
  return WT; 
}

} // namespace TauSpinner

