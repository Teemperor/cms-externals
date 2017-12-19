#include <iostream>
#include "TauSpinner/nonSM.h"
using std::cout;
using std::endl;

namespace TauSpinner {

/*****   GLOBAL VARIABLES   *****/

double (*nonSM_bornZ)(int, double, double, int, int, int) = default_nonSM_born;
double (*nonSM_bornH)(int, double, double, int, int, int) = default_nonSM_bornH;

extern int nonSM2;
extern int nonSMN;
extern double WTnonSM;
extern bool   IfHiggs;
/********************************/

double nonSM_born(int ID, double S, double cost, int H1, int H2, int key)
{
  if(IfHiggs) return nonSM_bornH(ID,S,cost,H1,H2,key);
  else        return nonSM_bornZ(ID,S,cost,H1,H2,key);
}


double default_nonSM_born(int ID, double S, double cost, int H1, int H2, int key)
{
  cout<<"TauSpinner::default_nonSM_born: this function is dummy\n"<<endl;
  cout<<"             user must provide his own nonSM_born"<<endl;
  cout<<"             see: nonSM_adopt() and set_nonSM_born( NULL )"<<endl;
  cout<<"             in TauSpinner/examples/tau-reweight-test.cxx for details"<<endl;
  exit(-1);
  return 1.0;
}


double default_nonSM_bornH(int ID, double S, double cost, int H1, int H2, int key)
{
  cout<<"TauSpinner::default_nonSM_born: this function is dummy\n"<<endl;
  cout<<"             user must provide his own nonSM_born"<<endl;
  cout<<"             see: nonSM_adopt() and set_nonSM_born( NULL )"<<endl;
  cout<<"             in TauSpinner/examples/tau-reweight-test.cxx for details"<<endl;
  exit(-1);
  return 1.0;
}

void set_nonSM_born( double (*fun)(int, double, double, int, int, int) )
{
  if(fun==NULL) nonSM_bornZ = default_nonSM_born;
  else          nonSM_bornZ = fun;
}

void set_nonSM_bornH( double (*fun)(int, double, double, int, int, int) )
{
  if(fun==NULL) nonSM_bornH = default_nonSM_bornH;
  else          nonSM_bornH = fun;
}

/*******************************************************************************
 Probability of the helicity state 
*******************************************************************************/
double plzap2(int ide, int idf, double svar, double costhe)
{
  if(ide!=0)
  {
    if(idf > 0) initwk_(&ide,&idf,&svar);
    else
    {
      int mide=-ide;
      int midf=-idf;

      initwk_(&mide,&midf,&svar);
    }
  }

  int    zero =  0;
  double one  =  1.0;
  double mone = -1.0;
  double ret  =  0.0;

  if(nonSM2==0)
  {
    ret =  t_born_(&zero,&svar,&costhe, &one, &one)
         /(t_born_(&zero,&svar,&costhe, &one, &one)
          +t_born_(&zero,&svar,&costhe,&mone,&mone));
  }
  else if(nonSM2==1)
  {
    ret =  nonSM_born(ide,svar,costhe, 1, 1, nonSM2)
         /(nonSM_born(ide,svar,costhe, 1, 1, nonSM2)
          +nonSM_born(ide,svar,costhe,-1,-1, nonSM2));

  // test of user prepared born cross section can be prepared here. 
  // convention between t_born and nonSM_born in choice of flavours
  // sign of costhe helicity signs etc have to be performed using SM version
  // of nonSM_born. Matching (up to may be overall s-dependent factor) between
  // t_born and nonSM_born must be achieved, see Section 4 for details
  // on technical tests. 
  DEBUG(
    double sm=  t_born_(&zero,&svar,&costhe, &one, &one)
              /(t_born_(&zero,&svar,&costhe, &one, &one)
               +t_born_(&zero,&svar,&costhe,&mone,&mone));
    double nsm= nonSM_born(ide,svar,costhe, 1, 1, nonSM2)
              /(nonSM_born(ide,svar,costhe, 1, 1, nonSM2)
               +nonSM_born(ide,svar,costhe,-1,-1, nonSM2));
    double smn= nonSM_born(ide,svar,costhe, 1, 1, 0     )
              /(nonSM_born(ide,svar,costhe, 1, 1, 0     )
               +nonSM_born(ide,svar,costhe,-1,-1, 0     ));

    cout<<"test of nonSM Born nonsm2="<<nonSM2 << endl;
    cout<<"ide,svar,costhe="<<ide <<" " << svar <<" "  << costhe << endl;
    cout<<"sm="<<sm <<" sm (new)="<<smn <<" nsm="<<nsm << endl;
    cout<<"sm and sm (new) should be essentially equal" << endl << endl;
    if (IfHiggs) cout << "(for Higgs we need to improve algorithm)" << endl;

  )

  }
  return ret;
}

double plweight(int ide, double svar, double costhe)
{
  if(nonSM2==0) return 1.0;
  if (ide==0 && !IfHiggs) return 1.0;

  double ret      =( nonSM_born(ide,svar,costhe, 1, 1, nonSM2)/svar
                    +nonSM_born(ide,svar,costhe,-1,-1, nonSM2)/svar)
                  /( nonSM_born(ide,svar,costhe, 1, 1, 0)/svar
		    +nonSM_born(ide,svar,costhe,-1,-1, 0)/svar);   // svar is introduced for future use
  // another option angular dependence only. Effect on x-section removed
  if(nonSMN==1) ret = ret * plnorm(ide,svar);

  return ret;
}

double plnorm(int ide, double svar)
{
  if(nonSMN==0) return 1.0;

  double c1 = 1.0/sqrt(3.0);
  double c2 = sqrt(2.0/3.0);

  double alpha  = 2*nonSM_born(ide,svar,0.0, 1, 1,0)+
                  2*nonSM_born(ide,svar,0.0,-1,-1,0);
  double beta   =   nonSM_born(ide,svar, c1, 1, 1,0) + nonSM_born(ide,svar,-c1, 1, 1,0)+
                    nonSM_born(ide,svar, c1,-1,-1,0) + nonSM_born(ide,svar,-c1,-1,-1,0);
  double gamma  =   nonSM_born(ide,svar, c2, 1, 1,0) + nonSM_born(ide,svar,-c2, 1, 1,0)+
                    nonSM_born(ide,svar, c2,-1,-1,0) + nonSM_born(ide,svar,-c2,-1,-1,0);
  double ret = ( alpha + 0.9*(gamma+alpha-2*beta) + 0.5*(4*beta-3*alpha-gamma) );

  alpha  = 2*nonSM_born(ide,svar,0.0, 1, 1,nonSM2)+
           2*nonSM_born(ide,svar,0.0,-1,-1,nonSM2);
  beta   =   nonSM_born(ide,svar, c1, 1, 1,nonSM2) + nonSM_born(ide,svar,-c1, 1, 1,nonSM2)+
             nonSM_born(ide,svar, c1,-1,-1,nonSM2) + nonSM_born(ide,svar,-c1,-1,-1,nonSM2);
  gamma  =   nonSM_born(ide,svar, c2, 1, 1,nonSM2) + nonSM_born(ide,svar,-c2, 1, 1,nonSM2)+
             nonSM_born(ide,svar, c2,-1,-1,nonSM2) + nonSM_born(ide,svar,-c2,-1,-1,nonSM2);

  ret = ret / ( alpha + 0.9*(gamma+alpha-2*beta) + 0.5*(4*beta-3*alpha-gamma) );

  return ret;
}

void nonSMHcorrPol(double S, SimpleParticle &tau1, SimpleParticle &tau2,
                   double *corrX2, double *polX2)
{  // tau+ and tau- in lab frame
  Particle tau_plus (    tau1.px(),    tau1.py(),    tau1.pz(),    tau1.e(),    tau1.pdgid() );
  Particle tau_minus(    tau2.px(),    tau2.py(),    tau2.pz(),    tau2.e(),    tau2.pdgid() );

  // P_QQ = sum of tau+ and tau- in lab frame
  Particle P_QQ( tau_plus.px()+tau_minus.px(), tau_plus.py()+tau_minus.py(), tau_plus.pz()+tau_minus.pz(), tau_plus.e()+tau_minus.e(), 0 );
  
  Particle P_B1(0, 0, 1, 1, 0);
  Particle P_B2(0, 0,-1, 1, 0);

  tau_plus. boostToRestFrame(P_QQ);
  tau_minus.boostToRestFrame(P_QQ);
  P_B1.     boostToRestFrame(P_QQ);
  P_B2.     boostToRestFrame(P_QQ);
  
  double costheta1 = (tau_plus.px()*P_B1.px()    +tau_plus.py()*P_B1.py()    +tau_plus.pz()*P_B1.pz()    ) /
                 sqrt(tau_plus.px()*tau_plus.px()+tau_plus.py()*tau_plus.py()+tau_plus.pz()*tau_plus.pz()) /
                 sqrt(P_B1.px()    *P_B1.px()    +P_B1.py()    *P_B1.py()    +P_B1.pz()    *P_B1.pz()    );

  double costheta2 = (tau_minus.px()*P_B2.px()    +tau_minus.py()*P_B2.py()    +tau_minus.pz()*P_B2.pz()    ) /
                 sqrt(tau_minus.px()*tau_minus.px()+tau_minus.py()*tau_minus.py()+tau_minus.pz()*tau_minus.pz()) /
                 sqrt(P_B2.px()    *P_B2.px()    +P_B2.py()    *P_B2.py()    +P_B2.pz()    *P_B2.pz()    );
               
  double sintheta1 = sqrt(1-costheta1*costheta1);
  double sintheta2 = sqrt(1-costheta2*costheta2);
  
  // Cosine of hard scattering
  double costhe = (costheta1*sintheta2 + costheta2*sintheta1) / (sintheta1 + sintheta2);
  
  // Invariant mass of tau+tau- pair
  double SS     = S; // other option is for tests: P_QQ.recalculated_mass()*P_QQ.recalculated_mass();
  

  // corrX2 calculation
  // polX2 calculation
  // WTnonSM calculation
  *corrX2 = -1.0;  // default
  *polX2  =  0.0;  // default
  WTnonSM =  1.0;  // default
}

} // namespace TauSpinner
