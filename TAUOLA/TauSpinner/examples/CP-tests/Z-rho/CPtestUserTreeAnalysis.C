//#include "UserTreeAnalysis.H"   // remove if copied to user working directory
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "MC4Vector.H"
#include "HEPParticle.H"
#include "TH1.h"
#include "Setup.H"
#include "TObjArray.h"
#include "TMath.h"

using namespace std;


// very similar to  MC_FillUserHistogram from MC-TESTER/src/Generate.cxx
// it fills or defines (if missing) new TH1D histogram (root)
inline void fillUserHisto(char *name,double val, double weight=1.0, 
                          double min=Setup::bin_min[0][0], 
                          double max=Setup::bin_max[0][0]){

    TH1D *h=(TH1D*)(Setup::user_histograms->FindObject(name));
    if(!h){
      // define new histogram (if missing)
      h=new TH1D(name,name,Setup::nbins[0][0],min,max);
      if(!h) return;
      Setup::user_histograms->Add(h);
      //   printf("user histogram created %s\n", name);
    }
    // fill histogram
    h->Fill(val,weight);

}

// print 3-vector. Note that 'v' is in fact 4-vector with energy= v[3] and momenta= v[0] - v[2] because
// of loading with fortran libraries
void print(double * v){
  cout << "("<<v[0]<<","<<v[1]<<","<<v[2]<<",E: "<<v[3]<<")"<<endl;
}

// calculates vector product of 3-vectors:  result = (v1 x v2), normalising it to unity
// function returns normalisation factor for optional use.
// useful for calculation of normal vector to the plane spaned on (v1, v2)
double normalised_cross_product(double * v1, double * v2, double * result){
  result[0] = v1[1]*v2[2] - v1[2]*v2[1];
  result[1] = v1[2]*v2[0] - v1[0]*v2[2];
  result[2] = v1[0]*v2[1] - v1[1]*v2[0];

  double normalisation = sqrt(result[0]*result[0]
			      +result[1]*result[1]
			      +result[2]*result[2]);

  for(int i=0; i<3; i++)
    result[i]=result[i]/normalisation;

  return normalisation;
}
// scalar product of 3-vectors
double dot_product(double *v1, double *v2){
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}
// scalar product of 3-vectors: components of MC4Vector v1, v2
double dot_product(MC4Vector v1, MC4Vector v2){
  return v1.GetX1()*v2.GetX1()+v1.GetX2()*v2.GetX2()+v1.GetX3()*v2.GetX3();
}

// length of 3-vector v
double magnitude(double *v){
  return sqrt(dot_product(v,v));
}


/** Main function for user provided analysis. Declared in SETUP.C to be  invoked by MC-Tester.
    It assumes the events are `something' -> tau+ tau-, then tau+- -> rho+- nu_tau -> pi+- pi0 nu_tau
    If other events are fed, code will fail (unpredictably).
    NOTE: variable names in the function below can be misleading!!
          (they were copied over from other analysis...) */
int CPtestUserTreeAnalysis(HEPParticle *mother,
			   HEPParticleList *stableDaughters, 
			   int nparams, double *params)
{
   // taking info from event into  local variables:
   // =============================================
   // all information is taken from list of stable dauhters
   // of `something'. We do not have taus, we do not know what is primary mother of each pi0.

    assert(mother!=0);
    assert(stableDaughters!=0);

    HEPParticleListIterator daughters(*stableDaughters);

    // four vectors for the pions.
    // We declare  two copies; for convenience of function types.  
    // Also for  simultaneous use of pion 4-momenta in rho-pair  frame and tau-s rest frames (optional).
    double pi_plus[4]={0};
    double pi_minus[4]={0};
    double pi0_plus[4]={0};
    double pi0_minus[4]={0};

    MC4Vector d_pi0_plus;
    MC4Vector d_pi0_minus;
    MC4Vector d_pi_plus;
    MC4Vector d_pi_minus;

    // variables to store the 4-vector of the rho+ rho- pair center of mass syst.
    double cm_px=0;
    double cm_py=0;
    double cm_pz=0;
    double cm_e=0;

    // variables to store tau+ 4-momentum
    double taup_px=0;
    double taup_py=0;
    double taup_pz=0;
    double taup_e=0;

    // variables to store tau- 4-momentum
    double taum_px=0;
    double taum_py=0;
    double taum_pz=0;
    double taum_e=0;



    // loop over all daughters and sort them by type. Copy into local variables
    // resconstruct tau+- and cm  4-momenta 
    for (HEPParticle *part=daughters.first(); part!=0; part=daughters.next()){

      if(part->GetPDGId()!=16&&part->GetPDGId()!=-16){
	// sum decay product into rho-pair 4-vector (neutrinos are to be excluded)
	cm_px+=part->GetPx();
	cm_py+=part->GetPy();
	cm_pz+=part->GetPz();
	cm_e+=part->GetE();
      }
      MC4Vector p4=part->GetP4();
      // identify decay products ad sum them to the corresponding tau+ tau- momenta  
      switch(part->GetPDGId()){  
         case 211:
	   d_pi_plus.Set(&p4);
	   d_pi_plus.SetM(part->GetM());   
	   taup_px+=part->GetPx();
	   taup_py+=part->GetPy();
	   taup_pz+=part->GetPz();
	   taup_e+=part->GetE();
	   break;
         case  16:
	   taum_px+=part->GetPx();
	   taum_py+=part->GetPy();
	   taum_pz+=part->GetPz();
	   taum_e+=part->GetE();
	   break;
         case -16:
	   taup_px+=part->GetPx();
	   taup_py+=part->GetPy();
	   taup_pz+=part->GetPz();
	   taup_e+=part->GetE();
	   break;
         case -211:
	   d_pi_minus.Set(&p4);
	   d_pi_minus.SetM(part->GetM());
	   taum_px+=part->GetPx();
	   taum_py+=part->GetPy();
	   taum_pz+=part->GetPz();
	   taum_e+=part->GetE();
	   break;
         case 111: //fill the pi0's at this moment we do not know if it belongs to tau+ or tau-
	   if(d_pi0_minus.GetX1()==0&&d_pi0_minus.GetX2()==0&&d_pi0_minus.GetX3()==0){
	     d_pi0_minus.Set(&p4);
	     d_pi0_minus.SetM(part->GetM());
	   }
	   else{
	     d_pi0_plus.Set(&p4);
	     d_pi0_plus.SetM(part->GetM());
	   }	
	   break;
      }
    }

    // Now check which pi0 is associated with
    // which pi+/-. Use the angle to decide.
    double costheta1 = dot_product(d_pi_plus,d_pi0_plus)/(d_pi_plus.Length()*d_pi0_plus.Length());
    double costheta2 = dot_product(d_pi_minus,d_pi0_plus)/(d_pi_minus.Length()*d_pi0_plus.Length());
    if(costheta1<costheta2){ //and if necessary, swap the pi0 vectors
      MC4Vector temp;
      temp.Set(&d_pi0_plus);
      temp.SetM(d_pi0_plus.GetM());
      d_pi0_plus.Set(&d_pi0_minus);
      d_pi0_plus.SetM(d_pi0_minus.GetM());
      d_pi0_minus.Set(&temp);
      d_pi0_minus.SetM(temp.GetM());
    }

    
    taup_px+=d_pi0_plus.GetX1();
    taup_py+=d_pi0_plus.GetX2();
    taup_pz+=d_pi0_plus.GetX3();
    taup_e+=d_pi0_plus.GetX0();        

    taum_px+=d_pi0_minus.GetX1();
    taum_py+=d_pi0_minus.GetX2();
    taum_pz+=d_pi0_minus.GetX3();
    taum_e+=d_pi0_minus.GetX0();        
    
    double taup_mass = sqrt(taup_e*taup_e-taup_px*taup_px-taup_py*taup_py-taup_pz*taup_pz);
    double taum_mass = sqrt(taum_e*taum_e-taum_px*taum_px-taum_py*taum_py-taum_pz*taum_pz);
    //   cout<<"taup_mass ="<<taup_mass<<endl;
    //  cout<<"taum_mass ="<<taum_mass<<endl;

    //Now boost pions into the rho+ rho- center of mass frame.
    double cm_mass = sqrt(cm_e*cm_e-cm_px*cm_px-cm_py*cm_py-cm_pz*cm_pz);

    d_pi0_plus.Boost(cm_px,cm_py,cm_pz,cm_e,cm_mass);
    d_pi0_minus.Boost(cm_px,cm_py,cm_pz,cm_e,cm_mass);
    d_pi_plus.Boost(cm_px,cm_py,cm_pz,cm_e,cm_mass);
    d_pi_minus.Boost(cm_px,cm_py,cm_pz,cm_e,cm_mass);
 
    pi0_plus[0]=d_pi0_plus.GetX1();
    pi0_plus[1]=d_pi0_plus.GetX2();
    pi0_plus[2]=d_pi0_plus.GetX3();
    pi0_plus[3]=d_pi0_plus.GetX0();  

    pi_plus[0]=d_pi_plus.GetX1();
    pi_plus[1]=d_pi_plus.GetX2();
    pi_plus[2]=d_pi_plus.GetX3();
    pi_plus[3]=d_pi_plus.GetX0();

    pi0_minus[0]=d_pi0_minus.GetX1();
    pi0_minus[1]=d_pi0_minus.GetX2();
    pi0_minus[2]=d_pi0_minus.GetX3();
    pi0_minus[3]=d_pi0_minus.GetX0();

    pi_minus[0]=d_pi_minus.GetX1();
    pi_minus[1]=d_pi_minus.GetX2();
    pi_minus[2]=d_pi_minus.GetX3();
    pi_minus[3]=d_pi_minus.GetX0();


    /******* calculate acoplanarity (theta) *****/
    //normal to the plane spanned by pi+ pi0 
    double n_plus[3];
    normalised_cross_product(pi_plus,pi0_plus,n_plus);

    //normal to the plane spanned by pi- pi0
    double n_minus[3];
    normalised_cross_product(pi_minus,pi0_minus,n_minus);
    // for convention like in paper K. Desch et al  hep-ph/0307331 use instead: 
    //    normalised_cross_product(pi0_minus,pi_minus,n_minus);

    // calculate the acoplanarity angle
    double theta = acos(dot_product(n_plus,n_minus));
    
    // extend definition of acoplanarity theta to lie between 0 and 2 PI
    // that is to angle between oriented half-planes
    double pi_minus_n_plus = dot_product(pi_minus,n_plus)/magnitude(pi_minus);    
    if(pi_minus_n_plus>0)
      theta=2*M_PI-theta;


    /*********** calculate C/D reco  (y1y2 in the paper) ***/

    //boost vectors back to the lab frame
    d_pi0_plus.Boost(-cm_px,-cm_py,-cm_pz,cm_e,cm_mass);
    d_pi_plus.Boost(-cm_px,-cm_py,-cm_pz,cm_e,cm_mass);
    d_pi0_minus.Boost(-cm_px,-cm_py,-cm_pz,cm_e,cm_mass);
    d_pi_minus.Boost(-cm_px,-cm_py,-cm_pz,cm_e,cm_mass);

    // For future work
    // ===============
    // construct effective tau 4 vectors from its visible products
    // At present this make no sense for TauSpinner examples,
    // but may be useful for optimization
    // 

    /*
    // taus must be already in  Higgs rest frame
    double e_tau = 125.0/2.0;
    double m_tau = 1.777;
    double p_mag_tau = sqrt(e_tau*e_tau - m_tau*m_tau); 

    MC4Vector p_tau_plus = d_pi_plus + d_pi0_plus;
    MC4Vector p_tau_minus = d_pi_minus + d_pi0_minus;

    double norm_plus = p_mag_tau/p_tau_plus.Length();
    double norm_minus = p_mag_tau/p_tau_minus.Length();

    p_tau_plus.SetX0(e_tau);
    p_tau_plus.SetX1(norm_plus*p_tau_plus.GetX1());
    p_tau_plus.SetX2(norm_plus*p_tau_plus.GetX2());
    p_tau_plus.SetX3(norm_plus*p_tau_plus.GetX3());
    p_tau_plus.SetM(m_tau);

    p_tau_minus.SetX0(e_tau);
    p_tau_minus.SetX1(norm_minus*p_tau_minus.GetX1());
    p_tau_minus.SetX2(norm_minus*p_tau_minus.GetX2());
    p_tau_minus.SetX3(norm_minus*p_tau_minus.GetX3());
    p_tau_minus.SetM(m_tau);
    */

    // Prepare options for calculation of y1, y2 
    // =========================================
    // Boost pions to the mother tau's frames.
    // We do not use first following option as default, because it is
    // difficult/impossible(?) for LHC applications:
    // Second option, where true tau momenta are used, may be useful 
    // for start of future studies.

 
    // Boost to the (reconstructed) tau's frames. It is a first option!
    /*
    d_pi0_plus.BoostP(p_tau_plus);
    d_pi_plus.BoostP(p_tau_plus);
    d_pi0_minus.BoostP(p_tau_minus);
    d_pi_minus.BoostP(p_tau_minus);
    */

    // Boost to the (generated) tau's frames. It is a second option!
    /*
    d_pi0_plus.Boost(taup_px,taup_py,taup_pz,taup_e,taup_mass);
    d_pi0_minus.Boost(taum_px,taum_py,taum_pz,taum_e,taum_mass);
    d_pi_plus.Boost(taup_px,taup_py,taup_pz,taup_e,taup_mass);
    d_pi_minus.Boost(taum_px,taum_py,taum_pz,taum_e,taum_mass);
    */


    // Calculate y1 and y2
    // ===================
    // The y1, y2 variables are the differences of charged/neutral pions energies
    // coming out from rho+- decays.
    // It can be in laboratory, H or corresponding tau rest frame.
    // Depending on the option we gain on sensitivity, but we never loose it completely.
    double y1=(d_pi_plus.GetX0()-d_pi0_plus.GetX0())/(d_pi_plus.GetX0()+d_pi0_plus.GetX0());
    double y2=(d_pi_minus.GetX0()-d_pi0_minus.GetX0())/(d_pi_minus.GetX0()+d_pi0_minus.GetX0());

    // We force using weighted events, 'WT = params[0]' == TauSpinner weight transmitted to user analysis
    //                                                                       like this one
    // ==============================

    double WT = 1.0;
    if( nparams>=1 ) WT = params[0];
    else {
         cout<<"CPtestUserTreeAnalysis: provide weight as parameter by adding"<<endl<<endl
             <<"    Setup::UTA_nparams   = 1;"<<endl
             <<"    Setup::UTA_params[0] = WT;"<<endl<<endl
             <<"to tau-reweight-test.cxx just before call to MC_Analyze"<<endl<<endl;
         exit(-1);
    }

    // Fill histogram
    // ===============
    // Now we have all necessary information, we fill standard histograms
    // of acoplanarity separating events into two categories depending on the 
    // sign of y1*y2. 
    // At this step we have ready information for more sophisticated observables
    // as well. For the start we can use all three variables: y1, y2, theta and weight WT
    // as input eg. for Neural Network discriminator

    char *plotname = new char[30];
    if(y1*y2>0)
      sprintf(plotname,"acoplanarity-angle-C");
    else
      sprintf(plotname,"acoplanarity-angle-D");
    fillUserHisto(plotname,theta,WT,0,2*M_PI);
    delete plotname;

    return 0;
};

