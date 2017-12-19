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


void  analMake(SimpleParticle X,  vector<SimpleParticle> tau1_daughters,vector<SimpleParticle> tau2_daughters, TH1D* hist1,  TH1D* hist2, TH1D* hist3, double WT);


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
  //  setHiggsParametersTR( 1.0,-1.0, 0.0, 0.0); // for pseudo-scalar H
  double theta=0.2;
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
  TH1D* hist111 = new TH1D("hist111","acollinarity" ,100,0.0,  M_PI);
  TH1D* hist112 = new TH1D("hist112","acollinarity" ,100,3.0,  M_PI);
  TH1D* hist113 = new TH1D("hist113","acomplanarity",100,0.0,2*M_PI);



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

    analMake(QQ, tau1_daughters, tau2_daughters, hist111, hist112, hist113, WT);

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
  hist113->Write();
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

void  analMake(SimpleParticle sp_QQ, vector<SimpleParticle> sp_tau1_daughters, vector<SimpleParticle> sp_tau2_daughters, TH1D* hist111,  TH1D* hist112, TH1D* hist113,double WT){

  // Create Particles from SimpleParticles

  Particle P_QQ     (     sp_QQ.px(),      sp_QQ.py(),      sp_QQ.pz(),      sp_QQ.e(), 0 );

  vector<Particle> tau_daughters;

  // Create vector of (all) tau daughters
  for(unsigned int i=0; i<sp_tau1_daughters.size(); i++)
  {
    Particle pp(sp_tau1_daughters[i].px(),
                sp_tau1_daughters[i].py(),
                sp_tau1_daughters[i].pz(),
                sp_tau1_daughters[i].e(),
                sp_tau1_daughters[i].pdgid() );

    tau_daughters.push_back(pp);
  }
  for(unsigned int i=0; i<sp_tau2_daughters.size(); i++)
  {
    Particle pp(sp_tau2_daughters[i].px(),
                sp_tau2_daughters[i].py(),
                sp_tau2_daughters[i].pz(),
                sp_tau2_daughters[i].e(),
                sp_tau2_daughters[i].pdgid() );

    tau_daughters.push_back(pp);
  }

  // boost to rest frame
  for(unsigned int i=0; i<tau_daughters.size();i++)
    tau_daughters[i].boostToRestFrame(P_QQ);



  // taking info from event into  local variables:
  // =============================================
  // all information is taken from list of stable dauhters
  // of `something'. We do not have taus.

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

  for(unsigned int i=0; i<tau_daughters.size();i++){
      //tau's 3-momenta are calculated from their daughters
      //(daughters are pi+/- or neutrino)
    if(tau_daughters[i].pdgid()==211||tau_daughters[i].pdgid()==-16){
      tau_plus[0]+=tau_daughters[i].px();
      tau_plus[1]+=tau_daughters[i].py();
      tau_plus[2]+=tau_daughters[i].pz();
    }
    if(tau_daughters[i].pdgid()==-211||tau_daughters[i].pdgid()==16){
      tau_minus[0]+=tau_daughters[i].px();
      tau_minus[1]+=tau_daughters[i].py();
      tau_minus[2]+=tau_daughters[i].pz();
    }
    //fill pi+ or pi- momenta array
    if(tau_daughters[i].pdgid()==-211){
      pi_minus[0]=tau_daughters[i].px();
      pi_minus[1]=tau_daughters[i].py();
      pi_minus[2]=tau_daughters[i].pz();
    }
    if(tau_daughters[i].pdgid()==211){
      pi_plus[0]=tau_daughters[i].px();
      pi_plus[1]=tau_daughters[i].py();
      pi_plus[2]=tau_daughters[i].pz();
    }
  }

  // construct effective beam direction  
  Particle P_beam(0, 0, 1, 1, 0);
  // boost to the rest frame of decaying resonance
  P_beam.     boostToRestFrame(P_QQ);

  // take beam space components
  beam_one[0]=P_beam.px();
  beam_one[1]=P_beam.py();
  beam_one[2]=P_beam.pz();
  
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
  //  cout<<"tansverse cos="<<tr_cos<<endl;
  
  /*** Acollinarity **/
  // calculate the acollinearity angle for the pi+ and pi- (in Higgs rest frame)
  double delta = acos(dot_product(pi_plus,pi_minus)/(magnitude(pi_plus)*magnitude(pi_minus)));
  //  std::cout << "delta = " << delta << std::endl;

  if(tr_cos>0.5 && P_QQ.recalculated_mass() > 60.0){    
  // this cut selects perpendicular configurations otherwise Rxx cancels Ryy and no transverse spin effects are visible.

    // Fill histograms
    // ===============
    hist111->Fill(delta,WT);
    hist112->Fill(delta,WT);

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
  
  // Fill histograms
  // ===============
  hist113->Fill(theta,WT);
 
}
