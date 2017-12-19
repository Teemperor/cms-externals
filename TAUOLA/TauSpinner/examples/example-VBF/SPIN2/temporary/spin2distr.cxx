#include "TauSpinner/tau_reweight_lib.h"
#include "TauSpinner/Particle.h"
#include "TauSpinner/spin2distr.h"
#include <cstdlib>
#include <string.h>
#include <iomanip>
#include <fstream>
using namespace Tauolapp;

namespace TauSpinner {

extern double CMSENE;
extern int    nonSM2;
extern int    relWTnonSM;
extern double WTnonSM;   
extern double Polari;
extern double WTamplit;
extern double WTamplitP;
extern double WTamplitM;
extern double CMSENE;
extern double f(double x, int ID, double SS, double cmsene);

const bool  DEBUG = 1;
  //  int _QCDdefault=1;
  // int _QCDvariant=1;

/** @brief Pointer to a function that modify rezult of alpha_ calc. in spin2dist
 *  
 *  By default this function is not used (pointer is NULL).
 *  It can be changed by the user through TauSpinner::set_spin2distrModif
 */
//void (*alphasModif)(double, int, int) = NULL;

/** @brief Set spin2distrModif function */
//  void set_alphasModif(void (*function)(double, int, int) )
//{
  //   alphasModif = function;
  // }

/** @brief Pointer to a function that modify rezult of Matrix Element calc. in spin2dist
 *  
 *  By default this function is not used (pointer is NULL).
 *  It can be changed by the user through TauSpinner::set_spin2distrModif
 */
double (*spin2distrModif)(int, int, int, int, int, int, double[6][4], int, double) = NULL;

/** @brief Set spin2distrModif function */
void set_spin2distrModif(double (*function)(int, int, int, int, int, int, double[6][4], int, double) )
{
    spin2distrModif = function;
}


/** Wrapper to VBDISTR and place for interface to user provided modification*/
double spin2distr(int I1, int I2, int I3, int I4, int H1, int H2, double P[6][4], int KEY)
{
  double P_copy[6][4];
  double original_result = 0.;
  
  memcpy(P_copy,P,sizeof(double)*6*4);

  if(  KEY<2) {
    // SM mode
    return  spin2distr_(&I1, &I2, &I3, &I4, &H1, &H2, P_copy, &KEY); 
  }
  else if(  !spin2distrModif) {
    printf("TauSpinner::spin2distr: User function  spin2distrModif not declared. Setting WT_contrib = 0.0 Failed attempt with KEY = %i5. \n",KEY);
    return original_result;
  }
  else if(  KEY<4) {
    // modification mode
    int KEY_BUF=KEY-2;
    original_result=  spin2distr_(&I1, &I2, &I3, &I4, &H1, &H2, P_copy, &KEY_BUF);
    return spin2distrModif(I1,I2,I3,I4,H1,H2,P_copy,KEY,original_result);
  }
  else {
    // replacement mode
    return spin2distrModif(I1,I2,I3,I4,H1,H2,P_copy,KEY,original_result);
  }
}


/** Get SPIN2 ME2
    Returns array W[2][2] */
//---------------------------------------------------------------------------------------------------------------------------------------------------
void getME2SPIN2(SimpleParticle &p3i, SimpleParticle &p4i, SimpleParticle &sp_X,SimpleParticle &tau1i, SimpleParticle &tau2i, double (&W)[2][2], int KEY)
{
  
  // this may be necessary because of matrix element calculation may require absolute energy-momentum conservation!
  // FSR photons may need to be treated explicitely or with interpolation procedures.
  Particle p3(p3i.px(),p3i.py(),p3i.pz(),p3i.e(),0);
  Particle p4(p4i.px(),p4i.py(),p4i.pz(),p4i.e(),0);
  double   mtaul=sqrt(tau1i.e()*tau1i.e()-tau1i.px()*tau1i.px()-tau1i.py()*tau1i.py()-tau1i.pz()*tau1i.pz());
  Particle tau1(tau1i.px(),tau1i.py(),tau1i.pz(),tau1i.e(),0);
  Particle tau2(tau2i.px(),tau2i.py(),tau2i.pz(),tau2i.e(),0);




  // we may want to force p3,p4 to be masless.
  Particle v1(0.,0., 1.,1.,0.);
  Particle v2(0.,0.,-1.,1.,0.);

  Particle P_QQ( p3.px()+p4.px()+tau1.px()+tau2.px()+1e-8*p3.pz(),
                 p3.py()+p4.py()+tau1.py()+tau2.py()-2e-8*p3.pz(),
                 p3.pz()+p4.pz()+tau1.pz()+tau2.pz(),
                 p3.e() +p4.e() +tau1.e() +tau2.e(), 0 );
  tau1.boostToRestFrame(P_QQ);
  tau2.boostToRestFrame(P_QQ);
  p3.boostToRestFrame(P_QQ);
  p4.boostToRestFrame(P_QQ);
  v1.boostToRestFrame(P_QQ);
  v2.boostToRestFrame(P_QQ);

  // Now we can define vers1 vers2
  double xn1=sqrt(v1.px()*v1.px()+v1.py()*v1.py()+v1.pz()*v1.pz());
  double xn2=sqrt(v2.px()*v2.px()+v2.py()*v2.py()+v2.pz()*v2.pz());
  double xn12=sqrt( (v1.px()-v2.px())*(v1.px()-v2.px()) +(v1.py()-v2.py())*(v1.py()-v2.py())+(v1.pz()-v2.pz())*(v1.pz()-v2.pz()));

  double SS = P_QQ.recalculated_mass()*P_QQ.recalculated_mass(); 
  
  double x1x2  = SS/CMSENE/CMSENE;
  double x1Mx2 = P_QQ.pz()/CMSENE*2;
  
  double x1 = (  x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;
  double x2 = ( -x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;
 
  //---------------------------------------------------------------------------
  // Construct the matrix for FORTRAN function
  // NOTE: different order of indices than in FORTRAN!
  // four options for partonic beams.
  double P[6][4] = { { sqrt(SS)/2,   sqrt(SS)/2/xn12*(v1.px()-v2.px()),    sqrt(SS)/2/xn12*(v1.py()-v2.py()),   sqrt(SS)/2/xn12*(v1.pz()-v2.pz())  },
		     { sqrt(SS)/2,  -sqrt(SS)/2/xn12*(v1.px()-v2.px()),   -sqrt(SS)/2/xn12*(v1.py()-v2.py()),  -sqrt(SS)/2/xn12*(v1.pz()-v2.pz())  },
		     // double P[6][4] = { { sqrt(SS)/2,  -sqrt(SS)/2/xn1*v2.px(),   -sqrt(SS)/2/xn1*v2.py(),  -sqrt(SS)/2/xn1*v2.pz()  },
		     //                    { sqrt(SS)/2,   sqrt(SS)/2/xn1*v2.px(),    sqrt(SS)/2/xn1*v2.py(),   sqrt(SS)/2/xn1*v2.pz()  },
		     // double P[6][4] = { { sqrt(SS)/2,   sqrt(SS)/2/xn1*v1.px(),    sqrt(SS)/2/xn1*v1.py(),   sqrt(SS)/2/xn1*v1.pz()  },
		     //                    { sqrt(SS)/2,  -sqrt(SS)/2/xn1*v1.px(),   -sqrt(SS)/2/xn1*v1.py(),  -sqrt(SS)/2/xn1*v1.pz()  },
		     // double P[6][4] = { { sqrt(SS)/2,                      0.0,                       0.0,              sqrt(SS)/2   },
		     //                    { sqrt(SS)/2,                      0.0,                       0.0,             -sqrt(SS)/2   }, 
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

  Particle P_tautau( tau1.px()+tau2.px(),  tau1.py()+tau2.py(), tau1.pz()+tau2.pz(), tau1.e()+tau2.e(),0  );

  double Q2  =  P_tautau.recalculated_mass()*P_tautau.recalculated_mass();
  double QQ2 =  P_QQ.recalculated_mass()*P_QQ.recalculated_mass();
  double PT2 =  P_QQ.px() * P_QQ.px() +  P_QQ.py()* P_QQ.py(); 

  double sumET =   p3.e()*sqrt( p3.px()*p3.px()+p3.py()*p3.py())/sqrt( p3.px()*p3.px()+p3.py()*p3.py()+p3.pz()*p3.pz())
                 + p4.e()*sqrt( p4.px()*p4.px()+p4.py()*p4.py())/sqrt( p4.px()*p4.px()+p4.py()*p4.py()+p4.pz()*p4.pz())
                 + tau1.e()*sqrt( tau1.px()*tau1.px()+tau1.py()*tau1.py())/sqrt( tau1.px()*tau1.px()+tau1.py()*tau1.py()+tau1.pz()*tau1.pz() )
                 + tau2.e()*sqrt( tau2.px()*tau2.px()+tau2.py()*tau2.py())/sqrt( tau2.px()*tau2.px()+tau2.py()*tau2.py()+tau2.pz()*tau2.pz() ) ;

  double sumMT =    sqrt( p3.recalculated_mass()*p3.recalculated_mass() + p3.px()*p3.px()+p3.py()*p3.py() )
                 +  sqrt( p4.recalculated_mass()*p4.recalculated_mass() + p4.px()*p4.px()+p4.py()*p4.py() )
                 +  sqrt( tau1.recalculated_mass()*tau1.recalculated_mass() + tau1.px()*tau1.px()+tau1.py()*tau1.py() )
                 +  sqrt( tau2.recalculated_mass()*tau2.recalculated_mass() + tau2.px()*tau2.px()+tau2.py()*tau2.py() );

  double fixed_scale = 200.;
 
  // reset to zero
  W[0][0]=0.0;
  W[0][1]=0.0;
  W[1][0]=0.0;
  W[1][1]=0.0;
  //add by marzieh  
  /*  if(DEBUG )
    {
      std::cout << " " <<  std::endl;
      std::cout << "---podstawiamy 4-pedy na jakies ---" <<  std::endl;
      std::cout << "ERW: ----------------------------- " << std::endl;

      
      //check point for gg > X2 > u u~
      double P1[6][4] = {{0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03}, 
			 {0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03},  
			 {0.5000000E+03,  0.1109243E+03,  0.4448308E+03, -0.1995529E+03},
			 {0.5000000E+03, -0.1109243E+03, -0.4448308E+03,  0.1995529E+03}};   
for (int I1=0;I1<=3;I1++){for (int I2=0;I2<=3;I2++){
      P[I1][I2]=P1[I1][I2];
    }}

    printf("  our event     :            P[0,i]=   %16.10f  %16.10f  %16.10f  %16.10f  \n",P[0][0],P[0][1],P[0][2],P[0][3]);
    int II=1;
    printf("                             P[1,i]=   %16.10f  %16.10f  %16.10f  %16.10f  \n",P[II][0],P[II][1],P[II][2],P[II][3]);
    II=2;
   
    printf("                             P[2,i]=   %16.10f  %16.10f  %16.10f  %16.10f  \n",P[II][0],P[II][1],P[II][2],P[II][3]);
    II=3;
    printf("                             P[3,i]=   %16.10f  %16.10f  %16.10f  %16.10f  \n",P[II][0],P[II][1],P[II][2],P[II][3]);
    printf(" ===========================\n");
//#################################################################
    KEY=0; //#################################################################
//#################################################################
    printf(" ============== KEY= %2i   ==\n",KEY);
    printf(" ===========================\n");
    printf(" \n");
    printf(" ===========================================================\n");
    printf(" ============== non-zero contributions to <|ME|^2>_spin   ==\n");
    printf(" ===========================================================\n");
    printf(" \n");

 
  
 for(int I1=-5;I1<=5;I1++){  // for test of single production process fix flavour
    for(int I2=-5;I2<=5;I2++){  // for test of single production process fix flavour
       for(int I3=-5;I3<=5;I3++){
	for(int I4=-5;I4<=5;I4++){
	  
	  
          int ID1 = I1; if( ID1 == 0 ) ID1 = 21;
          int ID2 = I2; if( ID2 == 0 ) ID2 = 21;
          int ID3 = I3; 
          int ID4 = I4; 
          W[0][0] += f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4, 1, -1, P, KEY);    // as in case of nonSM_adopt we change the sign for 
          W[0][1] += f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4, 1,  1, P, KEY);    // second tau helicity assuming without checking that
          W[1][0] += f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4,-1, -1, P, KEY);    // for SPIN2 quantization frame conventions are the same.
          W[1][1] += f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4,-1,  1, P, KEY);

	  if( ( W[0][0] > 10) || (W[0][1] > 10) || (W[1][0] > 10) || ( W[1][1] > 10)) { 
	      std::cout << "ERWxxx: ID= " <<ID1 << " " << ID2 << " " << ID3 << " " << ID4 << " W[]= " 
			<< f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4, 1, -1, P, KEY) << " " 
			<< f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4, 1,  1, P, KEY) << " "
			<< f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4,-1, -1, P, KEY) << " "
			<< f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4,-1,  1, P, KEY)
			<< std::endl; 
	    }
	  if( spin2distr(ID1,ID2,ID3,ID4, 1, -1, P, KEY)+spin2distr(ID1,ID2,ID3,ID4, -1, 1, P, KEY)+spin2distr(ID1,ID2,ID3,ID4, 1, 1, P, KEY) +spin2distr(ID1,ID2,ID3,ID4, -1, -1, P, KEY)> 0) {
	    float check = spin2distr(ID1,ID2,ID4,ID3, 1, -1, P, KEY)+spin2distr(ID1,ID2,ID3,ID4, -1, 1, P, KEY);
	      if( check != 0 ||  check == 0){ 
                int CPSTATUSICP=0;
		printf(" ID-s= %2i %2i %2i %2i  CP used= %1i ## VALUE:  %16.10e   ##  Spin contr.:  (+-)=  %9.3e  (-+)=  %9.3e  (--)=  %9.3e  (++)=  %9.3e  %9.3e \n", ID1, ID2, ID3,ID4, cpstatus_.icp,
		       spin2distr(ID1,ID2,ID3,ID4, 1, -1, P, KEY)+ spin2distr(ID1,ID2,ID3,ID4, -1, 1, P, KEY)+spin2distr(ID1,ID2,ID3,ID4, -1, -1, P, KEY)+ spin2distr(ID1,ID2,ID3,ID4, 1, 1, P, KEY),
		       spin2distr(ID1,ID2,ID3,ID4,  1, -1, P, KEY),
		       spin2distr(ID1,ID2,ID3,ID4, -1,  1, P, KEY),
		       spin2distr(ID1,ID2,ID3,ID4, -1, -1, P, KEY), 
		       spin2distr(ID1,ID2,ID3,ID4,  1,  1, P, KEY));
	    }
	  }

	}
      }
    }
 }
    }*/ 
  
  //end add by marzieh

    
  if( DEBUG ){
    // here you can overwrite kinematical configurations for tests
    std::cout << " " <<  std::endl;
    std::cout << "---podstawiamy 4-pedy na jakies ---" <<  std::endl;
    std::cout << "ERW: ----------------------------- " << std::endl;
    /* double P1[6][4]={{48.6174492572,  0.0000000220,  -0.0000000440,  48.6174492572},  
		    {48.6174492572,  -0.0000000220,  0.0000000440, -48.6174492572},  
		    {3.5929914991,   -3.3365984737,  0.5550238929,  -1.2118774680},  
		    {1.5522322367,   -1.3604813110,  0.5148003252,  -0.5417528648},  
		    {46.4880385810,  -8.4517422131,  -26.1544601820, 37.4498633797},  
		    {45.6016361975,  13.1488225544,  25.0846348505,  -35.6962330469}};  */

     double P1[6][4]  = {{0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03},
			{0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03},
			{0.8855009E+02, -0.2210038E+02,  0.4007979E+02, -0.7580437E+02},  
			{0.3283248E+03, -0.1038482E+03, -0.3019295E+03,  0.7649385E+02}, 
			{0.1523663E+03, -0.1058795E+03, -0.9770827E+02,  0.4954769E+02},  
			{0.4307588E+03,  0.2318280E+03,  0.3595580E+03, -0.5023717E+02}}; 
     /*  std::cout << "ERW: ----------------------------- " << std::endl;
  double P2[6][4]  = {{ 0.5000000E+03,  0.0,            0.0,            0.5000000E+03   },
                      { 0.5000000E+03,  0.0,            0.0,           -0.5000000E+03   }, 
                      { 0.1177462E+03,  -0.6070512E+02,   0.7123011E+02,   0.7145150E+02   }, 
                      { 0.3509495E+03,  -0.3178775E+02,   0.8393832E+02,   0.3392779E+03   },
                      { 0.3493321E+03,   0.1840069E+03,  -0.5152712E+02,  -0.2924315E+03   },
                      { 0.1819722E+03,  -0.9151401E+02,  -0.1036413E+03,  -0.1182978E+03   } };

  std::cout << "ERW: ----------------------------- " << std::endl;
  double P3[6][4]  = {{ 0.5000000E+03,  0.0,            0.0,            0.5000000E+03   },
                      { 0.5000000E+03,  0.0,            0.0,           -0.5000000E+03   }, 
                      { 0.2586900E+03,   0.1324670E+03,  -0.1696171E+03,  -0.1435378E+03   }, 
                      { 0.1084567E+03,  -0.5735712E+02,  -0.2162482E+02,  -0.8947281E+02   },
                      { 0.4005742E+03,  -0.1580760E+03,   0.3563160E+03,   0.9223569E+02   },
                      { 0.2322791E+03,   0.8296613E+02,  -0.1650741E+03,   0.1407749E+03   } };
   std::cout << "ERW: ----------------------------- " << std::endl;
  double P4[6][4]  = {{ 0.5000000E+03,  0.0,            0.0,            0.5000000E+03   },
                      { 0.5000000E+03,  0.0,            0.0,           -0.5000000E+03   }, 
                      { 0.1595700E+03,  -0.6917808E+02,  -0.1395175E+03,  -0.3481123E+02   }, 
                      { 0.2247758E+03,  -0.1360140E+03,   0.1650340E+03,  -0.6919641E+02   },
                      { 0.2508802E+03,   0.1447863E+01,   0.2499830E+03,  -0.2107335E+02   },
                      { 0.3647740E+03,   0.2037442E+03,  -0.2754995E+03,   0.1250810E+03   } };*/

  for (int I1=0;I1<=5;I1++){for (int I2=0;I2<=3;I2++){
      P[I1][I2]=P1[I1][I2];
    }}

    printf("  our event     :            P[0,i]=   %16.10f  %16.10f  %16.10f  %16.10f  \n",P[0][0],P[0][1],P[0][2],P[0][3]);
    int II=1;
    printf("                             P[1,i]=   %16.10f  %16.10f  %16.10f  %16.10f  \n",P[II][0],P[II][1],P[II][2],P[II][3]);
    II=2;
    printf("                             P[2,i]=   %16.10f  %16.10f  %16.10f  %16.10f  \n",P[II][0],P[II][1],P[II][2],P[II][3]);
    II=3;
    printf("                             P[3,i]=   %16.10f  %16.10f  %16.10f  %16.10f  \n",P[II][0],P[II][1],P[II][2],P[II][3]);
    II=4;
    printf("                             P[4,i]=   %16.10f  %16.10f  %16.10f  %16.10f  \n",P[II][0],P[II][1],P[II][2],P[II][3]);
    II=5;
    printf("                             P[5,i]=   %16.10f  %16.10f  %16.10f  %16.10f  \n",P[II][0],P[II][1],P[II][2],P[II][3]);
    printf(" ===========================\n");
//#################################################################
    KEY=0; //#################################################################
//#################################################################
    printf(" ============== KEY= %2i   ==\n",KEY);
    printf(" ===========================\n");
    printf(" \n");
    printf(" ===========================================================\n");
    printf(" ============== non-zero contributions to <|ME|^2>_spin   ==\n");
    printf(" ===========================================================\n");
    printf(" \n");

  }

     // these loops need to be cleaned from zero contributions! 
  for(int I1=-5;I1<=5;I1++){  // for test of single production process fix flavour
    for(int I2=-5;I2<=5;I2++){  // for test of single production process fix flavour
      for(int I3=-5;I3<=5;I3++){
	for(int I4=-5;I4<=5;I4++){
	  
          int ID1 = I1; if( ID1 == 0 ) ID1 = 21;
          int ID2 = I2; if( ID2 == 0 ) ID2 = 21;
          int ID3 = I3; if( ID3 == 0 ) ID3 = 21;
          int ID4 = I4; if( ID4 == 0 ) ID4 = 21;

	  W[0][0] += f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4, 1, -1, P, KEY);  // as in case of nonSM_adopt we change the sign for 
	  W[0][1] += f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4, 1,  1, P, KEY);  // second tau helicity assuming without checking that
	  W[1][0] += f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4,-1, -1, P, KEY); //for SPIN2 quantization frame conventions are the same
	  W[1][1] += f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4,-1,  1, P, KEY);

	  if( ( W[0][0] > 10) || (W[0][1] > 10) || (W[1][0] > 10) || ( W[1][1] > 10)) { 
	    std::cout << "ERWxxx: ID= " <<ID1 << " " << ID2 << " " << ID3 << " " << ID4 << " W[]= " 
		      << f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4, 1, -1, P, KEY) << " " 
		      << f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4, 1,  1, P, KEY) << " "
		      << f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4,-1, -1, P, KEY) << " "
		      << f_x1[I1]*f_x2[I2]*spin2distr(ID1,ID2,ID3,ID4,-1,  1, P, KEY)
		      << std::endl; 
	  }
	  if( spin2distr(ID1,ID2,ID3,ID4, 1, -1, P, KEY)+spin2distr(ID1,ID2,ID3,ID4, -1, 1, P, KEY)+spin2distr(ID1,ID2,ID3,ID4, 1, 1, P, KEY) +spin2distr(ID1,ID2,ID3,ID4, -1, -1, P, KEY)> 0) {
	    float check = spin2distr(ID1,ID2,ID4,ID3, 1, -1, P, KEY)+spin2distr(ID1,ID2,ID3,ID4, -1, 1, P, KEY);
	    if( check != 0 ||  check == 0){ 
	      int CPSTATUSICP=0;
	      /*    printf(" ID-s= %2i %2i %2i %2i  CP used= %1i ## VALUE:  %16.10e   ##  Spin contr.:  (+-)=  %9.3e  (-+)=  %9.3e  (--)=  %9.3e  (++)=  %9.3e  \n", ID1, ID2, ID3,ID4, cpstatus2_.icp,
		     spin2distr(ID1,ID2,ID3,ID4, 1, -1, P, KEY)+ spin2distr(ID1,ID2,ID3,ID4, -1, 1, P, KEY)+spin2distr(ID1,ID2,ID3,ID4, -1, -1, P, KEY)+ spin2distr(ID1,ID2,ID3,ID4, 1, 1, P, KEY),
		     spin2distr(ID1,ID2,ID3,ID4,  1, -1, P, KEY),
		     spin2distr(ID1,ID2,ID3,ID4, -1,  1, P, KEY),
		     spin2distr(ID1,ID2,ID3,ID4, -1, -1, P, KEY), 
		     spin2distr(ID1,ID2,ID3,ID4,  1,  1, P, KEY) );*/
	    
	        std::cout << ID1 << " "<< ID2 << " "<<ID3 << " "<<ID4 << " "<< std::setprecision(14) << spin2distr(ID1,ID2,ID3,ID4, 1, 1, P, KEY)<<'\n';
	      
	      
	    }
	  }
	}
	
      }
     
    }
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
double calculateWeightFromParticlesSPIN2(SimpleParticle &p3, SimpleParticle &p4,SimpleParticle &sp_X, SimpleParticle &sp_tau1, SimpleParticle &sp_tau2, vector<SimpleParticle> &sp_tau1_daughters, vector<SimpleParticle> &sp_tau2_daughters)
{
  SimpleParticle         sp_tau;
  SimpleParticle         sp_nu_tau;
  vector<SimpleParticle> sp_tau_daughters;
  int KEY=0;
  // here we impose that it is for SM Higgs only
  if(sp_X.pdgid()==25) KEY=1;

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
                                        // pre-initialization all are equal i.e. no spin effects
  
 
  // calculate ME^2*PDF's for SM
  getME2SPIN2(p3, p4, sp_X, sp_nu_tau, sp_tau, W, KEY);
  //  if(KEY==0 || KEY==1 )getME2SPIN2(p3, p4, sp_X, sp_tau1, sp_tau2, W, KEY); 
  //else                 getME2SPIN2(p3, p4, sp_X, sp_tau1, sp_tau2, W, 0);


  double sum=(W[0][0]+W[0][1]+ W[1][0]+W[1][1]); // getME2SPIN2 calculated PDF's times matrix elements squared for each of four helicity configurations  
                                                 // tau+ and tau- using SM process KEY=0 DY, KEY=1 Higgs. 

  if(nonSM2==0 ) { 
    WTnonSM=1.0;
    if(relWTnonSM==0)  WTnonSM=sum;  // The variable accessible for user stores production weight ME^2 * PDF's summed over flavours and spin states
  }
  if(nonSM2==1) { 
    // now we re-calculate W  using anomalous variant of matrix elements.
    getME2SPIN2(p3, p4, sp_X, sp_nu_tau, sp_tau, W, KEY+2);

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
