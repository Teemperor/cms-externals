/**
 * @author T. Przedzinski
 * @date 29 Nov 2013
 *
 * @brief Reweights datafile crom CLEO to RChL currents
 *
 */
// To run with MC-TESTER, uncomment this line or at compilation step add flag:
// -D_USE_MCTESTER_
//#define _USE_MCTESTER_

#ifdef _USE_MCTESTER_

// MC-TESTER headers
#include "Generate.h"
#include "HepMCEvent.H"
#include "Setup.H"

#endif

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

#include <iostream>
#include <fstream>
#include <cstdio>
using namespace std;
using namespace Tauolapp;
using namespace TauSpinner;

//-----------------------------------------------------------------------------
//- Main ----------------------------------------------------------------------
//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  if(argc<2)
  {
    cout<<"Usage:   "<<argv[0]<<" <input_file> [<events_limit>]"<<endl;
    exit(0);
  }

  //---------------------------------------------------------------------------
  //- Initialization ----------------------------------------------------------
  //---------------------------------------------------------------------------

  // Limit number of processed events
  int events_limit = 0;
  
  if(argc>2) events_limit = atoi(argv[2]);

  // Initialize LHAPDF
  string name="MSTW2008nnlo90cl.LHgrid";
  LHAPDF::initPDFSetByName(name);
  
  // Initialize Tauola turning RChL currents on
  Tauola::setNewCurrents(1);
  Tauola::initialize();

  // Initialize TauSpinner
  double CMSENE = 7000; // center of mass system energy (old default was 7000.0).
                        // used in PDF calculation for pp collisions only
  bool Ipp = true;      // 'true' for pp collisions
                        // false otherwise (not prepared yet)
  int Ipol   = 1;       // are input samples polarized? default answer was yes that is  1;

  // NOTE: For this analysis, values below should remain 0.
  //       For details see TauSpinner/examples/tau-reweight-test.cxx
  int nonSM2 = 0; // are we using nonSM calculations?
  int nonSMN = 0; // if we are using nonSM calculations we may want corrections 
                  // to shapes only: y/n  (1/0)

  initialize_spinner(Ipp, Ipol, nonSM2, nonSMN,  CMSENE);

  // Open I/O files
  HepMC::IO_GenEvent input_file(argv[1],std::ios::in);

#ifdef _USE_MCTESTER_
  // Initialize MC-Tester
  MC_Initialize();
#endif

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
    double WTp = 1.0;
    double WTm = 1.0;
    
    // NOTE: here we reweight CLEO currents to RChL using amplitude weight
    //       Methods 'calculateWeightFrom...' return spin weight which
    //       in this case is ignored
    if( abs(X.pdgid())==24 ||  abs(X.pdgid())==37 )
    {
    
      Tauola::setNewCurrents(0);

      calculateWeightFromParticlesWorHpn(X, tau, tau2, tau_daughters); // note that tau2 is tau neutrino 

      WTp = getWtamplitP();
      WTm = getWtamplitM();
      
      Tauola::setNewCurrents(1);
      
      calculateWeightFromParticlesWorHpn(X, tau, tau2, tau_daughters);
      
      WTp = getWtamplitP()/WTp;
      WTm = getWtamplitM()/WTm;
    }
    else if( X.pdgid()==25 || X.pdgid()==36 || X.pdgid()==22 || X.pdgid()==23 )
    {
      Tauola::setNewCurrents(0);

      calculateWeightFromParticlesH(X, tau, tau2, tau_daughters,tau_daughters2);

      WTp = getWtamplitP();
      WTm = getWtamplitM();
      
      Tauola::setNewCurrents(1);
      
      calculateWeightFromParticlesH(X, tau, tau2, tau_daughters,tau_daughters2);
      
      WTp = getWtamplitP()/WTp;
      WTm = getWtamplitM()/WTm;
    }
    else
    {
      cout<<"WARNING: Unexpected PDG for tau mother: "<<X.pdgid()<<endl;
    }
    
    wt_average += WTm;
    events_count++;

#ifdef _USE_MCTESTER_
    
    // We choose the weight depending on whether we look at tau+ or tau-
    // NOTE: this approach won't work if MC-TESTER is configured to
    //       look at e.g. Z or W decays
    if     (Setup::decay_particle== 15) WT = WTm;
    else if(Setup::decay_particle==-15) WT = WTp;
    else WT = 1.0;

    // For MC-TESTER analysis we use event defined in TauSpinner/src/readParticlesFromTAUOLA_HepMC
    HepMC::GenEvent *mc_tester_event = readParticlesFromTAUOLA_HepMC_getEvent();
    HepMCEvent *temp_event = new HepMCEvent(*mc_tester_event,false);
    MC_Analyze(temp_event,WT);

    // Cleanup
    delete temp_event;
#endif

    if(events_count%10000==0) cout<<"EVT: "<<events_count<<endl;
    if(events_limit && events_count>=events_limit) break;
  }

#ifdef _USE_MCTESTER_
  MC_Finalize();
#endif

  cout<<endl<<"No of events processed for reweighting: "<<events_count<<endl;
  cout<<      "WT average for these events: "<<wt_average/events_count<<endl;
}
