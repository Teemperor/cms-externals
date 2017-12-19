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
      //      printf("user histogram created %s\n", name);
    }
    // fill histogram
    h->Fill(val,weight);

}

// print 3-vector. Note that 'v' is in fact 4-vector with energy= v[3] and momenta= v[0] - v[2] because
// of loading with fortran libraries
void print(double * v){
  cout << "("<<v[0]<<","<<v[1]<<","<<v[2]<<")"<<endl;
}

// calculates vector product of 3-vectors: result = (v1 x v2), normalising it to unity
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

// length of 3-vector v
double magnitude(double *v){
  return sqrt(dot_product(v,v));
}


/** Main function for user provided analysis. Declared in SETUP.C to be  invoked by MC-Tester.
    It assume that the events are H/Z -> tau+ tau-
    and then tau -> pi nu. 
    If other events are fed, code will fail (unpredictably).
**/
int CPbenchPiUserTreeAnalysisA(HEPParticle *mother,
			 HEPParticleList *stableDaughters, 
			 int nparams, double *params)
{

  // taking info from event into  local variables:
  // =============================================
  // all information is taken from list of stable dauhters
  // of `something'. We do not have taus.

    assert(mother!=0);
    assert(stableDaughters!=0);

    HEPParticleListIterator daughters(*stableDaughters);

    // arrays to hold 3-momenta of tau+- and pi+-.
    // also for helpful things
    double tau_plus[3]={0};
    double tau_minus[3]={0};
    double pi_plus[3]={0};
    double pi_minus[3]={0};
    double beam_one[3]={0};
    double plane[3]={0};
    double plane_pi[3]={0};
    //loop over the daughters, boost them into the Higgs (or Z etc.) frame,
    //then fill the 3-momenta arrays, also for tau+-.
    for (HEPParticle *part=daughters.first(); part!=0; part=daughters.next()){
      MC4Vector d4(part->GetE(),part->GetPx(),
		   part->GetPy(),part->GetPz(),part->GetM());
      d4.Boost(mother->GetPx(),mother->GetPy(),mother->GetPz(),
	       mother->GetE(),mother->GetM());

      //tau's 3-momenta are calculated from their daughters
      //(daughters are pi+/- or neutrino)
      if(part->GetPDGId()==211||part->GetPDGId()==-16){
	tau_plus[0]+=d4.GetX1();
	tau_plus[1]+=d4.GetX2();
	tau_plus[2]+=d4.GetX3();
      }
      if(part->GetPDGId()==-211||part->GetPDGId()==16){
	tau_minus[0]+=d4.GetX1();
	tau_minus[1]+=d4.GetX2();
	tau_minus[2]+=d4.GetX3();
      }
      
      //fill pi+ or pi- momenta array
      if(part->GetPDGId()==-211){
	pi_minus[0]=d4.GetX1();
	pi_minus[1]=d4.GetX2();
	pi_minus[2]=d4.GetX3();
      }
      if(part->GetPDGId()==211){
	pi_plus[0]=d4.GetX1();
	pi_plus[1]=d4.GetX2();
	pi_plus[2]=d4.GetX3();
      }
    }
    // construct effective beam direction
      MC4Vector beam(1.0,0.0,0.0,1.0,0.0);
      beam.Boost(mother->GetPx(),mother->GetPy(),mother->GetPz(),
	       mother->GetE(),mother->GetM());

        // take beam space components
	beam_one[0]=beam.GetX1();
	beam_one[1]=beam.GetX2();
	beam_one[2]=beam.GetX3();

	// normal to reaction plane
	plane[0]=beam_one[1]*tau_minus[2]-beam_one[2]*tau_minus[1];
	plane[1]=beam_one[2]*tau_minus[0]-beam_one[0]*tau_minus[2];
	plane[2]=beam_one[0]*tau_minus[1]-beam_one[1]*tau_minus[0];

	// normal to visible products plane
	plane_pi[0]=tau_minus[1]*pi_minus[2]-tau_minus[2]*pi_minus[1];
	plane_pi[1]=tau_minus[2]*pi_minus[0]-tau_minus[0]*pi_minus[2];
	plane_pi[2]=tau_minus[0]*pi_minus[1]-tau_minus[1]*pi_minus[0];

    //calculate the angle between reaction plane and decay plane
    double tr_cos = fabs(dot_product(plane,plane_pi)/(magnitude(plane)*magnitude(plane_pi)));
    //    cout<<"tansverse cos="<<tr_cos<<endl;

    /*** Acollinarity **/
    // calculate the acollinearity angle for the pi+ and pi- (in Higgs rest frame)
    double delta = acos(dot_product(pi_plus,pi_minus)/(magnitude(pi_plus)*magnitude(pi_minus)));

    char *plotname = new char[30];
    double WT = 1.0;
    // We force using weighted events, 'WT = params[0]' == TauSpinner weight transmitted to user analysis
    //                                                                       like this one
    // ==============================
     if( nparams>=1 ) WT = params[0];
     else {
         cout<<"CPtestUserTreeAnalysis: provide weight as parameter by adding"<<endl<<endl
             <<"    Setup::UTA_nparams   = 1;"<<endl
             <<"    Setup::UTA_params[0] = WT;"<<endl<<endl
             <<"to tau-reweight-test.cxx just before call to MC_Analyze"<<endl<<endl;
         exit(-1);
    }
     if(tr_cos>0.5 && mother->GetM()>60.0){    // this cut selects perpendicular configurations otherwise Rxx cancels Ryy and no transverse spin effects are visible.

    // Fill histograms
    // ===============
	sprintf(plotname,"delta");
    fillUserHisto(plotname,delta,WT,0,M_PI);
	sprintf(plotname,"delta2");
    fillUserHisto(plotname,delta,WT,3,M_PI);
     }

    /*** Acoplanarity (theta) **/

    //calculate the angle between  pi+ and tau+
    double projection_plus = dot_product(pi_plus,tau_plus)/(magnitude(tau_plus)*magnitude(tau_plus));
    //calculate the transverse part of pi+ 3-momentum with respect to the tau+ direction
    double pi_plus_pt[3];
    pi_plus_pt[0]  = pi_plus[0]-(projection_plus*tau_plus[0]);
    pi_plus_pt[1]  = pi_plus[1]-(projection_plus*tau_plus[1]);
    pi_plus_pt[2]  = pi_plus[2]-(projection_plus*tau_plus[2]);
    
    //calculate the angle between  pi- and tau-
    double projection_minus = dot_product(pi_minus,tau_minus)/(magnitude(tau_minus)*magnitude(tau_minus));
    //calculate the transverse part of pi- 3-momentum with respect to the tau- direction
    double pi_minus_pt[3];
    pi_minus_pt[0]  = pi_minus[0]-(projection_minus*tau_minus[0]);
    pi_minus_pt[1]  = pi_minus[1]-(projection_minus*tau_minus[1]);
    pi_minus_pt[2]  = pi_minus[2]-(projection_minus*tau_minus[2]);

    //calculate the angle between the pi+ transverse and pi- transverse
    //but this  gives the angle from 0 to pi only.
    double theta = acos(dot_product(pi_plus_pt,pi_minus_pt)/(magnitude(pi_plus_pt)*magnitude(pi_minus_pt)));
    
    //to get the angle between 0 and 2 pi, use the normal to the plane of pi+ pi- transverse momenta
    //if the normal is in the direction of the tau-, set theta = 2pi - theta.
    double n3[3];
    normalised_cross_product(pi_plus_pt,pi_minus_pt,n3);    
    double theta_sign = dot_product(n3,tau_minus)/magnitude(tau_minus);    
    if(theta_sign>0)
      theta=2*M_PI-theta;
 
    // Fill histogram
     // ==============
	sprintf(plotname,"theta");
    fillUserHisto(plotname,theta,WT,0,2*M_PI);
    delete plotname;

    return 0;
};

