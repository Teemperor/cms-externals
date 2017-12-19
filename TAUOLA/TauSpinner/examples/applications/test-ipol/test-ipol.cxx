// HepMC IO_GenEvent header
#include "HepMC/IO_GenEvent.h"

// TAUOLA header
#include "Tauola/Tauola.h"

// LHAPDF header
#include "LHAPDF/LHAPDF.h"

// TauSpinner headers
#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"

// ROOT headers
#include "TH1D.h"
#include "TString.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <cstdio>
using namespace std;
using namespace Tauolapp;
using namespace TauSpinner;

double get_costheta(SimpleParticle &p1, SimpleParticle &p2)
{
// tau+ and tau- in lab frame
  Particle pp1( p1.px(), p1.py(), p1.pz(), p1.e(), p1.pdgid() );
  Particle pp2( p2.px(), p2.py(), p2.pz(), p2.e(), p2.pdgid() );

  double costheta1 =  pp1.pz()/sqrt(pp1.px()*pp1.px()+pp1.py()*pp1.py()+pp1.pz()*pp1.pz());

  double costheta2 = -pp2.pz()/sqrt(pp2.px()*pp2.px()+pp2.py()*pp2.py()+pp2.pz()*pp2.pz());

  double sintheta1 = sqrt(1-costheta1*costheta1);
  double sintheta2 = sqrt(1-costheta2*costheta2);

  return (costheta1*sintheta2 + costheta2*sintheta1) / (sintheta1 + sintheta2);
}

/* Get decay channel
   Returns:
   -1 - unrecognized (unsupported) channel
    0 - no daughters
    1 - electron (tau -> nu_tau nue e [+gamma])
    2 - muon     (tau -> nu_tau nu_mu mu [+gamma])
    3 - pion     (tau -> nu_tau pi)
    4 - rho      (tau -> nu_tau pi0 pi) */
int get_decay_channel(vector<SimpleParticle> &tau_daughters)
{
  int size = tau_daughters.size();

  if(size==0) return 0;

  int id[4] = { 0 };
  if(size>0) id[0] = abs(tau_daughters[0].pdgid());
  if(size>1) id[1] = abs(tau_daughters[1].pdgid());
  if(size>2) id[2] = abs(tau_daughters[2].pdgid());
  if(size>3) id[3] = abs(tau_daughters[3].pdgid());

  // nu_tau pi
  if(size==2 && id[1]==211) return 3;

  // nu_tau rho
  if(size==2 && id[1]==213) return 4;

  if(size==3)
  {
    // nu_tau e nu_e
    if(id[1]==11) return 1;

    // nu_tau mu nu_mu
    if(id[1]==13) return 2;

    // nu_tau pi0 pi (without intermediate rho)
    if( (id[1]==111 && id[2]==211) ||
        (id[1]==211 && id[2]==111)    ) return 4;
  }

  if(size==4 && id[3]==22)
  {
    // nu tau e nu_e gamma
    if(id[1]==11) return 1;

    // nu tau mu nu_mu gamma
    if(id[1]==13) return 2;
  }

  return -1;
}

//-----------------------------------------------------------------------------
//- Main ----------------------------------------------------------------------
//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  if(argc<4)
  {
    cout<<"Usage:   "<<argv[0]<<" <input_file> <output_rootfile> <ipol> [<events_limit>]"<<endl;
    exit(0);
  }

  //---------------------------------------------------------------------------
  //- Initialization ----------------------------------------------------------
  //---------------------------------------------------------------------------

  // Limit number of processed events
  int events_limit = 0;

  if(argc>4) events_limit = atoi(argv[4]);

  // Initialize LHAPDF
  string name="MSTW2008nnlo90cl.LHgrid";
  LHAPDF::initPDFSetByName(name);

  // Initialize Tauola
  Tauola::initialize();

  // Initialize TauSpinner
  double CMSENE = 8000.;      // center of mass system energy (old default was 7000.0).
                              // used in PDF calculation for pp collisions only
  bool Ipp = true;            // 'true' for pp collisions
                              // false otherwise (not prepared yet)
  int Ipol   = atoi(argv[3]); // are input samples polarized? default answer was yes that is  1;

  // NOTE: For this analysis, values below should remain 0.
  //       For details see TauSpinner/examples/tau-reweight-test.cxx
  int nonSM2 = 0; // are we using nonSM calculations?
  int nonSMN = 0; // if we are using nonSM calculations we may want corrections
                  // to shapes only: y/n  (1/0)

  initialize_spinner(Ipp, Ipol, nonSM2, nonSMN,  CMSENE);

  // Open I/O files
  HepMC::IO_GenEvent input_file(argv[1],std::ios::in);

  //---------------------------------------------------------------------------
  //- Plot definitions --------------------------------------------------------
  //---------------------------------------------------------------------------

  TFile out(argv[2],"RECREATE");

  TH1D *pi_e_all    = new TH1D("pi_e_all",   "pi_e all events ",100,0,1);
  TH1D *pi_e_all_WT = new TH1D("pi_e_all_WT","pi_e all events ",100,0,1);
  TH1D *pi_e_lt     = new TH1D("pi_e_lt",    "pi_e costheta'<0",100,0,1);
  TH1D *pi_e_lt_WT  = new TH1D("pi_e_lt_WT", "pi_e costheta'<0",100,0,1);
  TH1D *pi_e_gt     = new TH1D("pi_e_gt",    "pi_e costheta'>0",100,0,1);
  TH1D *pi_e_gt_WT  = new TH1D("pi_e_gt_WT", "pi_e costheta'>0",100,0,1);

  int    events_count = 0;
  double wt_average   = 0.0;

  //---------------------------------------------------------------------------
  //- Event loop --------------------------------------------------------------
  //---------------------------------------------------------------------------
  while(true)
  {
    SimpleParticle X, tau, tau2; // SimpleParticle consist of 4 momentum and PDGid.
    vector<SimpleParticle> tau_daughters, tau_daughters2;

    // Here, we read another event from input_file.

    // NOTE: for W+/- or H+/- tau2 contains neutrino and tau_daughters2 is empty
    // User may introduce his own method of initializing X, tau, tau2, tau_daughters, tau_daughters2.
    // Then it may be convenient to follow format of the
    // 'readParticlesFromTAUOLA_HepMC' example line below
    int status = readParticlesFromTAUOLA_HepMC(input_file, X, tau, tau2, tau_daughters, tau_daughters2);

    // Go to next one if there is nothing more to do with this event
    if( status==1 ) continue;

    // Finish if there is nothing more to read from the file
    if( status==0 ) break;

    // Sets polarization state of the input sample
    //setSpinOfSample(0);

    // BASE_WT can be a weight read from an event record
    // For now it is set to 1.0
    double BASE_WT = 1.0;
    double WT = 1.0;

    // Calculate weight
    if( abs(X.pdgid())==24 ||  abs(X.pdgid())==37 )
    {
      WT = calculateWeightFromParticlesWorHpn(X, tau, tau2, tau_daughters); // note that tau2 is tau neutrino
    }
    else if( X.pdgid()==25 || X.pdgid()==36 || X.pdgid()==22 || X.pdgid()==23 )
    {
      WT = calculateWeightFromParticlesH(X, tau, tau2, tau_daughters,tau_daughters2);
    }
    else
    {
      cout<<"WARNING: Unexpected PDG for tau mother: "<<X.pdgid()<<endl;
    }

    // Returns spin of the taus attributed during reweighting
    //double pol = getTauSpin();

    // At this step we have one weight: WT

    // Options:
    // WT=1.0/WT      // to unspin sample having spin
    // WT=(2-WT)/(WT) // to flip spin correlations between (H+- <-->  W+-)
    // WT=(2-WT)/(WT) // to flip spin correlations between (H <--> Z/gamma)
    //                   (this simple trick is valid with no Z couplings taken)
    //                   for Z calculation of weight for Higgs has to be done
    //                   directly
    // WARNING: these new WT will not be limited from above!

    int channel1 = get_decay_channel(tau_daughters);
    int channel2 = get_decay_channel(tau_daughters2);

    // Analyze only events that belong to supported decay channels
    if(channel1==3 && channel2==3)
    {
      double costheta  = get_costheta(tau_daughters[1],tau_daughters2[1]);
      double costheta2 = (tau.pz() + tau2.pz() < 0.0) ? -costheta : costheta;

      pi_e_all->   Fill(tau_daughters[1].e()/tau.e());
      pi_e_all_WT->Fill(tau_daughters[1].e()/tau.e(),WT);

      if(costheta2<0)
      {
        pi_e_lt->   Fill(tau_daughters[1].e()/tau.e());
        pi_e_lt_WT->Fill(tau_daughters[1].e()/tau.e(),WT);
      }
      else
      {
        pi_e_gt->   Fill(tau_daughters[1].e()/tau.e());
        pi_e_gt_WT->Fill(tau_daughters[1].e()/tau.e(),WT);
      }
    }

    if(events_count%10000==0) cout<<"EVT: "<<events_count<<endl;
    if(events_limit && events_count>=events_limit) break;

    events_count++;
    wt_average += WT;
  }

  out.Write();
  out.Close();

  cout<<endl<<"No of events processed for spin weight: "<<events_count<<endl;
  cout<<      "WT average for these events: "<<wt_average/events_count<<endl;
}
