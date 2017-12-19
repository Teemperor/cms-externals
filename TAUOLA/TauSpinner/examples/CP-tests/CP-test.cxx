// MC-TESTER headers
#include "Generate.h"
#include "HepMCEvent.H"
#include "Setup.H"

// TAUOLA header
#include "Tauola/Tauola.h"

#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"

// LHAPDF header
#include "LHAPDF/LHAPDF.h"
using namespace std;
using namespace Tauolapp;
using namespace TauSpinner;

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

  double CMSENE = 14000.; // center of mass system energy used in PDF calculation
  bool   Ipp    = true;   // for pp collisions
  int    Ipol   = 1;      // are input samples polarized?

  // Initialize TauSpinner (flags nonSM and nonSMN are 0 for these tests)
  initialize_spinner(Ipp, Ipol, 0, 0, CMSENE);

  // Initialize transverse spin effects of TauSpinner, these are studied in  CP-tests/H-rho, CP-tests/H-pi
  setHiggsParametersTR(-1.0, 1.0, 0.0, 0.0); // for scalar H
  //setHiggsParametersTR( 1.0,-1.0, 0.0, 0.0); // for pseudo-scalar H
  //double theta=0.2;
  //setHiggsParametersTR(-cos(2*theta),cos(2*theta) ,-sin(2*theta),-sin(2*theta)); // for mixed parity case


  // Initialize trasverse spin effects of TauSpinner, these are studied in  CP-tests/Z-pi, CP-tests/Z-rho
  // Multipliers for  components of transverse density matrix of DY
  //                  (Rxx,Ryy,Rxy,Ryx)
  setZgamMultipliersTR(1., 1., 1., 1. );

  // Open I/O files
  HepMC::IO_GenEvent input_file(input_filename,std::ios::in);

  // Initialize MC-Tester, SETUP.C (of executable directory) will configure options of MC-Tester
  // in particular user analysis macro can be  declared. 
  MC_Initialize();

  int    events_count = 0;
  double wt_average   = 0.0;

  //---------------------------------------------------------------------------
  //- Event loop --------------------------------------------------------------
  //---------------------------------------------------------------------------
  while(true) {
    double WT    = 1.0;
    int    pdgid = 0;

    SimpleParticle         X, tau, tau2; // SimpleParticle consist of 4 momentum and PDGid.
    vector<SimpleParticle> tau_daughters, tau_daughters2;

    // Read event from input_file.
    int status = readParticlesFromTAUOLA_HepMC(input_file, X, tau, tau2, tau_daughters, tau_daughters2);

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
      WT = calculateWeightFromParticlesWorHpn(X, tau, tau2, tau_daughters); // note that tau2 is tau neutrino
    }
    else if( X.pdgid()==25 || X.pdgid()==36 || X.pdgid()==22 || X.pdgid()==23 ) // case of X=H,A,gamma,Z
    {
      // NOTE: if any component of transverse density matrix for DY has been turned on
      //       using setZgamMultipliersTR, weight WT will contain transverse part
      //       of the spin amplitude (calculated from tables table1-1.txt table2-2.txt
      //       which must be stored in the same directory as the program)
      WT = calculateWeightFromParticlesH(X, tau, tau2, tau_daughters,tau_daughters2);
    }
    else
    {
      cout<<"WARNING: Unexpected PDG for tau mother: "<<X.pdgid()<<endl;
    }

    wt_average += WT;
    events_count++;

    // For MC-TESTER analysis we use event defined in TauSpinner/examples/readParticlesFromTAUOLA_HepMC.cxx
    HepMC::GenEvent *mc_tester_event = readParticlesFromTAUOLA_HepMC_getEvent();
    HepMCEvent *temp_event = new HepMCEvent(*mc_tester_event,false);

    // Setup needed for all CP-tests; this is for UserAnalysis macro invoked  by MC-Tester.
    Setup::UTA_nparams   = 1;
    Setup::UTA_params[0] = WT;

    // MC-Tester analysis:  WT is passed for default  MC-Tester distributions.  
    // this weight will be passed to all apearances of the decay branch under
    // consideration of the temp_event. In all cases weight WT will be used.
    // In fact there will be double counting in such a case. 
    MC_Analyze(temp_event,WT);

    // Cleanup
    delete temp_event;

    if(events_count%10000==0) cout<<"EVT: "<<events_count<<endl;
    if(events_limit && events_count>=events_limit) break;
  }

  // Finalize MC-Tester processing
  MC_Finalize();

  cout<<endl<<"No of events processed for spin weight: "<<events_count<<endl;
  cout<<      "WT average for these events: "<<wt_average/events_count<<endl;
}
