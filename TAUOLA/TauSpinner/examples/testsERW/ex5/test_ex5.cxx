// TAUOLA header
#include "Tauola/Tauola.h"

#include "TauSpinner/Particle.h"
#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/tau_reweight_lib.h"
#include "eventReader.h"

// LHAPDF header
#include "LHAPDF/LHAPDF.h"
using namespace std;
using namespace Tauolapp;
using namespace TauSpinner;

//ROOT headers
#include "TH1.h"   
#include "TFile.h"  


void  analMake(SimpleParticle X,  vector<SimpleParticle> tau1_daughters,vector<SimpleParticle> tau2_daughters, TH1D* hist1,  TH1D* hist2, double WT);


int main(int argc, char **argv) {

  char *input_filename = "events.dat";

  if(argc<2)
  {
    cout<<endl<<"Processing all available events from default input file: "<<input_filename<<endl;
    cout<<      "You can change this option by using:   "<<argv[0]<<" <filename> [<events_limit>]"<<endl<<endl;
  }
  else input_filename = argv[1];

  //---------------------------------------------------------------------------
  //- Initialization ----------------------------------------------------------
  //---------------------------------------------------------------------------

  // Limit number of processed events
  int events_limit = 0;

  if(argc>2) events_limit = atoi(argv[2]);

  // Initialize Tauola
  Tauola::initialize();

  string name="cteq6ll.LHpdf";
  LHAPDF::initPDFSetByName(name);

  // Next 3 lines are used to initialize TauSpinner
  //  CMSENE should be adjusted, Ipp    = true should be kept and
  // examples are set to work with unpolarized sammples only, Ipol is not
  // fully checked yet.

  double CMSENE = 13000.; // center of mass system energy used in PDF calculation
  bool   Ipp    = true;   // for pp collisions
  int    Ipol   = 1;      // are input samples polarized?

  // Initialize TauSpinner (flags nonSM and nonSMN are 0 for these tests)
  initialize_spinner(Ipp, Ipol, 0, 0, CMSENE);

  // Initialize transverse spin effects of TauSpinner, these are studied in  CP-tests/H-rho, CP-tests/H-pi
  // setHiggsParametersTR(-1.0, 1.0, 0.0, 0.0); // for scalar H
  // setHiggsParametersTR( 1.0,-1.0, 0.0, 0.0); // for pseudo-scalar H
  double theta=0.0;
  setHiggsParametersTR(-cos(2*theta),cos(2*theta) ,-sin(2*theta),-sin(2*theta)); // for mixed parity case


  // Initialize trasverse spin effects of TauSpinner, these are studied in  CP-tests/Z-pi, CP-tests/Z-rho
  // Multipliers for  components of transverse density matrix of DY
  //                  (Rxx,Ryy,Rxy,Ryx)
  setZgamMultipliersTR(1., 1., 1., 1. );

  // Open I/O files
  HepMC::IO_GenEvent input_file(input_filename,std::ios::in);

  // Initialize MC-Tester, SETUP.C (of executable directory) will configure options of MC-Tester
  // in particular user analysis macro can be  declared. 

  int    events_count = 0;
  double wt_average   = 0.0;

  // initialise histograms 
  TH1D* hist100 = new TH1D("hist100","spint WT weight" ,100,0.0,  10.0);
  TH1D* hist111 = new TH1D("hist111","acollinarity" ,100,0.0,  2*M_PI);
  TH1D* hist112 = new TH1D("hist112","acollinarity" ,100,0.0,  2*M_PI);



  //---------------------------------------------------------------------------
  //- Event loop --------------------------------------------------------------
  //---------------------------------------------------------------------------
  while(true) {
    double WT    = 1.0;
    int    pdgid = 0;

    SimpleParticle         X, tau1, tau2; // SimpleParticle consist of 4 momentum and PDGid.
    vector<SimpleParticle> tau1_daughters, tau2_daughters;

    // Read event from input_file.
    int status = read_HepMC(input_file, X, tau1, tau2, tau1_daughters, tau2_daughters);

    SimpleParticle QQ;
    QQ.setPx( tau1.px()+tau2.px());
    QQ.setPy( tau1.py()+tau2.py());
    QQ.setPz( tau1.pz()+tau2.pz());
    QQ.setE ( tau1.e()+tau2.e());
    QQ.setPdgid(0);

    // Go to next one if there is nothing more to do with this event
    if( status==1 ) continue;

    // Finish if there is nothing more to read from the file
    if( status==0 ) break;

    // pdgid is for our main object X. X can be Higgs Z etc... 
    pdgid = X.pdgid();

    // Declares polarization state of the input sample (sample can have events where spin effects
    // are (are not) taken into account. It can be also included in  part only. 
    //setSpinOfSample(0);

    // TauSpinner calculates spin weight for the event:
    if( abs(X.pdgid())==24 ||  abs(X.pdgid())==37 ) // case of X=W or charged Higgs
    {
      WT = calculateWeightFromParticlesWorHpn(X, tau1, tau2, tau1_daughters); // note that tau2 is tau neutrino
    }
    else if( X.pdgid()==25 || X.pdgid()==36 || X.pdgid()==22 || X.pdgid()==23 ) // case of X=H,A,gamma,Z
    {
      // NOTE: if any component of transverse density matrix for DY has been turned on
      //       using setZgamMultipliersTR, weight WT will contain transverse part
      //       of the spin amplitude (calculated from tables table1-1.txt table2-2.txt
      //       which must be stored in the same directory as the program)
      WT = calculateWeightFromParticlesH(X, tau1, tau2, tau1_daughters,tau2_daughters);
    }
    else
    {
      cout<<"WARNING: Unexpected PDG for tau mother: "<<X.pdgid()<<endl;
    }

    hist100->Fill(WT);

    wt_average += WT;
    events_count++;

    analMake(QQ, tau1_daughters, tau2_daughters, hist111, hist112, WT);

    if(events_count%10000==0) cout<<"EVT: "<<events_count<<endl;
    if(events_limit && events_count>=events_limit) break;
  }

  // write out histograms

  char output_filename[500]="";
  strcat (output_filename,input_filename);
  strcat (output_filename,".root");
  TFile outFile(output_filename, "recreate");

  outFile.cd();
  hist100->Write();
  hist111->Write();
  hist112->Write();
  outFile.Close();

  cout<<endl<<"No of events processed for spin weight: "<<events_count<<endl;
  cout<<      "WT average for these events: "<<wt_average/events_count<<endl;
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

void  analMake(SimpleParticle sp_QQ, vector<SimpleParticle> sp_tau1_daughters, vector<SimpleParticle> sp_tau2_daughters, TH1D* hist111,  TH1D* hist112, double WT){

  // Create Particles from SimpleParticles

  vector<Particle> tau_daughters;

  // Create vector of (all) tau daughters
  for(unsigned int i=0; i<sp_tau1_daughters.size(); i++)
  {
    Particle pp(sp_tau1_daughters[i].px(),
                sp_tau1_daughters[i].py(),
                sp_tau1_daughters[i].pz(),
                sp_tau1_daughters[i].e(),
                sp_tau1_daughters[i].pdgid() );
    //    pp.print();

    tau_daughters.push_back(pp);
  }
  for(unsigned int i=0; i<sp_tau2_daughters.size(); i++)
  {
    Particle pp(sp_tau2_daughters[i].px(),
                sp_tau2_daughters[i].py(),
                sp_tau2_daughters[i].pz(),
                sp_tau2_daughters[i].e(),
                sp_tau2_daughters[i].pdgid() );
    //    pp.print();

    tau_daughters.push_back(pp);
  }


  // taking info from event into  local variables:
  // =============================================
  // all information is taken from list of stable dauhters
  // of `something'. We do not have taus.


    // four vectors for the pions.
    // We declare  two copies; for convenience of function types.  
    // Also for  simultaneous use of pion 4-momenta in rho-pair  frame and tau-s rest frames (optional).
    double pi_plus[4]={0};
    double pi_minus[4]={0};
    double pi0_plus[4]={0};
    double pi0_minus[4]={0};

    double temp_pi0_plus[4]={0};
    double temp_pi0_minus[4]={0};

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


  //loop over the daughters, boost them into the Higgs (or Z etc.) frame,
  //then fill the 3-momenta arrays, also for tau+-.

  for(unsigned int i=0; i<tau_daughters.size();i++){
    if(tau_daughters[i].pdgid()!=-16 && tau_daughters[i].pdgid()!=16){
      // sum decay product into rho-pair 4-vector (neutrinos excluded)
      cm_px+=tau_daughters[i].px();
      cm_py+=tau_daughters[i].py();
      cm_pz+=tau_daughters[i].pz();
      cm_e+=tau_daughters[i].e();
    }

    // identify decay products ad sum them to the corresponding tau+ tau- momenta  
    switch(tau_daughters[i].pdgid()){
    case 211:
	   pi_plus[0]  = tau_daughters[i].px();
	   pi_plus[1] = tau_daughters[i].py();
	   pi_plus[2] = tau_daughters[i].pz();
	   pi_plus[3] = tau_daughters[i].e();
	   taup_px+=tau_daughters[i].px();
	   taup_py+=tau_daughters[i].py();
	   taup_pz+=tau_daughters[i].pz();
	   taup_e+=tau_daughters[i].e();
	   break;
    case 16:
	   taum_px+=tau_daughters[i].px();
	   taum_py+=tau_daughters[i].py();
	   taum_pz+=tau_daughters[i].pz();
	   taum_e+=tau_daughters[i].e();
	   break;
    case -16:
	   taup_px+=tau_daughters[i].px();
	   taup_py+=tau_daughters[i].py();
	   taup_pz+=tau_daughters[i].pz();
	   taup_e+=tau_daughters[i].e();
	   break;
    case -211:
	   pi_minus[0]  = tau_daughters[i].px();
	   pi_minus[1] = tau_daughters[i].py();
	   pi_minus[2] = tau_daughters[i].pz();
	   pi_minus[3] = tau_daughters[i].e();
	   taum_px+=tau_daughters[i].px();
	   taum_py+=tau_daughters[i].py();
	   taum_pz+=tau_daughters[i].pz();
	   taum_e+=tau_daughters[i].e();
	   break;
    case 111: //fill the pi0's at this moment we do not know if it belongs to tau+ or tau-
      if(temp_pi0_minus[0]==0 && temp_pi0_minus[1]==0 &&  temp_pi0_minus[1]==0){
	temp_pi0_minus[0] = tau_daughters[i].px();
	temp_pi0_minus[1] = tau_daughters[i].py();
	temp_pi0_minus[2] = tau_daughters[i].pz();
	temp_pi0_minus[3] = tau_daughters[i].e();
      }
      else{
	temp_pi0_plus[0] = tau_daughters[i].px();
	temp_pi0_plus[1] = tau_daughters[i].py();
	temp_pi0_plus[2] = tau_daughters[i].pz();
	temp_pi0_plus[3] = tau_daughters[i].e();
      }	
      break; 
    }
  }
  // Now check which pi0 is associated with
  // which pi+/-. Use the angle to decide.
  double costheta1 = dot_product(pi_plus,temp_pi0_plus)/(magnitude(pi_plus)*magnitude(temp_pi0_plus));
  double costheta2 = dot_product(pi_minus,temp_pi0_plus)/(magnitude(pi_minus)*magnitude(temp_pi0_plus));
  if(costheta1<costheta2){ //and if necessary, swap the pi0 vectors
    double temp[4]={0};
    temp[0]=temp_pi0_plus[0];
    temp[1]=temp_pi0_plus[1];
    temp[2]=temp_pi0_plus[2];
    temp[3]=temp_pi0_plus[3];
    temp_pi0_plus[0]=temp_pi0_minus[0];
    temp_pi0_plus[1]=temp_pi0_minus[1];
    temp_pi0_plus[2]=temp_pi0_minus[2];
    temp_pi0_plus[3]=temp_pi0_minus[3];
    temp_pi0_minus[0]=temp[0];
    temp_pi0_minus[1]=temp[1];
    temp_pi0_minus[2]=temp[2];
    temp_pi0_minus[3]=temp[3];
  }
  
  double taup_mass = sqrt(taup_e*taup_e-taup_px*taup_px-taup_py*taup_py-taup_pz*taup_pz);
  double taum_mass = sqrt(taum_e*taum_e-taum_px*taum_px-taum_py*taum_py-taum_pz*taum_pz);
  //   cout<<"taup_mass ="<<taup_mass<<endl;
  //  cout<<"taum_mass ="<<taum_mass<<endl;
  
  //Now define rest-frame of rho-rho system and express all as Particle 
  Particle P_rhorho(cm_px,cm_py,cm_pz,cm_e,0);
  
  Particle P_lab_pi0_plus(temp_pi0_plus[0],temp_pi0_plus[1],temp_pi0_plus[2],temp_pi0_plus[3],0);
  Particle P_lab_pi0_minus(temp_pi0_minus[0],temp_pi0_minus[1],temp_pi0_minus[2],temp_pi0_minus[3],0);
  Particle P_lab_pi_plus(pi_plus[0],pi_plus[1],pi_plus[2],pi_plus[3],0);
  Particle P_lab_pi_minus(pi_minus[0],pi_minus[1],pi_minus[2],pi_minus[3],0);
  
  Particle P_pi0_plus(temp_pi0_plus[0],temp_pi0_plus[1],temp_pi0_plus[2],temp_pi0_plus[3],0);
  Particle P_pi0_minus(temp_pi0_minus[0],temp_pi0_minus[1],temp_pi0_minus[2],temp_pi0_minus[3],0);
  Particle P_pi_plus(pi_plus[0],pi_plus[1],pi_plus[2],pi_plus[3],0);
  Particle P_pi_minus(pi_minus[0],pi_minus[1],pi_minus[2],pi_minus[3],0);
  
  //  P_pi_plus.print();
  //  P_pi_minus.print();
  //  P_pi0_plus.print();
  //  P_pi0_minus.print();
  
  //Now boost pions into the rho+ rho- center of mass frame.
  P_pi0_plus.boostToRestFrame(P_rhorho);
  P_pi0_minus.boostToRestFrame(P_rhorho);
  P_pi_plus.boostToRestFrame(P_rhorho);
  P_pi_minus.boostToRestFrame(P_rhorho);
  
  //Now back to table representation
  pi_plus[0]=P_pi_plus.px();
  pi_plus[1]=P_pi_plus.py();
  pi_plus[2]=P_pi_plus.pz();
  pi_plus[3]=P_pi_plus.e();
  
  pi_minus[0]=P_pi_minus.px();
  pi_minus[1]=P_pi_minus.py();
  pi_minus[2]=P_pi_minus.pz();
  pi_minus[3]=P_pi_minus.e();
  
  pi0_plus[0]=P_pi0_plus.px();
  pi0_plus[1]=P_pi0_plus.py();
  pi0_plus[2]=P_pi0_plus.pz();
  pi0_plus[3]=P_pi0_plus.e();
  
  pi0_minus[0]=P_pi0_minus.px();
  pi0_minus[1]=P_pi0_minus.py();
  pi0_minus[2]=P_pi0_minus.pz();
  pi0_minus[3]=P_pi0_minus.e();
  
  
  /******* calculate acoplanarity (theta) *****/
  //normal to the plane spanned by pi+ pi0 
  double n_plus[3];
  normalised_cross_product(pi_plus,pi0_plus,n_plus);
  
  //normal to the plane spanned by pi- pi0
  double n_minus[3];
  normalised_cross_product(pi_minus,pi0_minus,n_minus);

  // calculate the acoplanarity angle
  double acoplanarity = acos(dot_product(n_plus,n_minus));
  
  // extend definition of acoplanarity angle to lie between 0 and 2 PI
  // that is to angle between oriented half-planes
  double pi_minus_n_plus = dot_product(pi_minus,n_plus)/magnitude(pi_minus);    
  if(pi_minus_n_plus>0)
    acoplanarity=2*M_PI-acoplanarity;
  
  // Calculate y1 and y2
  // ===================
  // The y1, y2 variables are the differences of charged/neutral pions energies
  // coming out from rho+- decays.
  // It can be in laboratory, H or corresponding tau rest frame.
  // Depending on the option we gain on sensitivity, but we never loose it completely.
  double y1=(P_lab_pi_plus.e()-P_lab_pi0_plus.e())/(P_lab_pi_plus.e()+P_lab_pi0_plus.e());
  double y2=(P_lab_pi_minus.e()-P_lab_pi0_minus.e())/(P_lab_pi_minus.e()+P_lab_pi0_minus.e());
  
  // Fill histogram
  // ===============
  // Now we have all necessary information, we fill standard histograms
  // of acoplanarity separating events into two categories depending on the 
  // sign of y1*y2. 
  // At this step we have ready information for more sophisticated observables
  // as well. For the start we can use all three variables: y1, y2, theta and weight WT
  // as input eg. for Neural Network discriminator
  
  if(y1*y2>0) {    
    hist111->Fill(acoplanarity,WT);
  } else {
    hist112->Fill(acoplanarity,WT);
  }
  
}

