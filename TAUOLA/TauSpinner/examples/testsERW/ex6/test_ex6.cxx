/**
 * @author T. Przedzinski
 * @date 7 Jul 2011 22:30 pm
 */

// HepMC IO_GenEvent header
#include "HepMC/IO_GenEvent.h"

// TAUOLA header
#include "Tauola/Tauola.h"
#include "Tauola/TauolaParticlePair.h"

#include "TauSpinner/Particle.h"
#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/tau_reweight_lib.h"
#include "TauSpinner/nonSM.h"

#include "eventReader.h"

// LHAPDF header
#include "LHAPDF/LHAPDF.h"
using namespace std;
using namespace Tauolapp;
using namespace TauSpinner;

//ROOT headers
#include "TH1.h"   
#include "TFile.h"  

double getZPolFixFlavour(SimpleParticle &sp_X, SimpleParticle &sp_B1,  SimpleParticle &sp_B2, SimpleParticle &sp_tau, SimpleParticle &sp_nu_tau);
double getCosThetaStar(SimpleParticle &sp_B1,  SimpleParticle &sp_B2, SimpleParticle &sp_tau, SimpleParticle &sp_nu_tau);
 

// This is an example of the program reweighting unpolarized tau events.
int main(int argc, char **argv) {

  char *input_filename = "events.dat";
  
  if(argc<2)
  {
    cout<<endl<<"Processing all available events from default input file: "<<input_filename<<endl;
    cout<<      "You can change this option by using:   "<<argv[0]<<" <filename> [<events_limit>]"<<endl<<endl;
  }
  else input_filename = argv[1];

  //---------------------------------------------------------------------------
  //- Initialization -----------------d-----------------------------------------
  //---------------------------------------------------------------------------

  // Limit number of processed events
  int events_limit = 0;
  
  if(argc>2) events_limit = atoi(argv[2]);
  
  // Initialize Tauola
  Tauola::initialize();

  string name="MSTW2008nnlo90cl.LHgrid";
  LHAPDF::initPDFSetByName(name);

  double CMSENE = 8000.0; // center of mass system energy.
                           // used in PDF calculation. For pp collisions only
  bool Ipp = true;         // for pp collisions
  int Ipol = 1;            // are input samples polarized?
  int nonSM2 = 0;          // are we using nonSM calculations?
                           // at present we have checked only that for nonSM2 = 0
                           // all works as in the past. nonSM2 = 1 may simply mean 
                           // errors are switched on. More work is needed 
                           // for this option.
  int nonSMN = 0;          // If we are using nonSM calculations we may want corrections 
                           // to shapes only: y/n  (1/0)

  
  // Initialize TauSpinner
  initialize_spinner(Ipp, Ipol, nonSM2, nonSMN,  CMSENE);

  // Open I/O files  (in our example events are taken from "events.dat")
  HepMC::IO_GenEvent input_file(input_filename,std::ios::in);

  int    events_count = 0;
  double wt_average   = 0.0;

  // create histograms
  TH1D *hist100 = new TH1D("hist100"," wtspin weight",100,0.0,10.0);
  TH1D *hist103 = new TH1D("hist103"," tau helicity  ",5,-2.0,2.0);
  TH1D *hist203 = new TH1D("hist203"," tau helicity, (wtspin)  ",5,-2.0,2.0);
  TH1D *hist303 = new TH1D("hist303"," tau helicity, (1/wtspin)  ",5,-2.0,2.0);
  TH1D *hist_wt = new TH1D("hist_wt"," wt A4  ",50,0.9,1.1);
  TH1D *hist103_wt = new TH1D("hist103_wt"," tau helicity  ",5,-2.0,2.0);
  TH1D *hist104    = new TH1D("hist104"," tau fix-flavour helicity  ",5,-2.0,2.0);
  TH1D *hist104_wt = new TH1D("hist104_wt"," tau fix-flavour helicity  ",5,-2.0,2.0);

  TH1D *hist_costhetaCS          = new TH1D("hist_costhetaCS",         " costheta CS          ",50,-1.0,1.0);
  hist_costhetaCS->Sumw2();
  TH1D *hist_costhetaCS_wt       = new TH1D("hist_costhetaCS_wt",      " costheta CS wt       ",50,-1.0,1.0);
  hist_costhetaCS_wt->Sumw2();
  TH1D *hist_costhetaCS_wt_ratio = new TH1D("hist_costhetaCS_wt_ratio"," costheta CS wt ratio ",50,-1.0,1.0);
  hist_costhetaCS_wt_ratio->Sumw2();
  TH1D *hist_costhetaCS_pol       = new TH1D("hist_costhetaCS_pol",      " costheta CS wt       ",50,-1.0,1.0);
  hist_costhetaCS_pol->Sumw2();
  TH1D *hist_costhetaCS_pol_ratio       = new TH1D("hist_costhetaCS_pol_ratio",      " pol ratio      ",50,-1.0,1.0);
  hist_costhetaCS_pol_ratio->Sumw2();
  TH1D *hist_costhetaCS_pol_wt = new TH1D("hist_costhetaCS_pol_wt"," costheta CS pol ratio wt ",50,-1.0,1.0);
  hist_costhetaCS_pol_wt->Sumw2();
  TH1D *hist_costhetaCS_pol_wt_ratio = new TH1D("hist_costhetaCS_pol_wt_ratio"," costheta CS pol ratio wt ",50,-1.0,1.0);
  hist_costhetaCS_pol_wt_ratio->Sumw2();

  TH1D *hist_costhetaStar          = new TH1D("hist_costhetaStar",         " costhetaStar          ",50,-1.0,1.0);
  hist_costhetaStar->Sumw2();
  TH1D *hist_polFixFlavour         = new TH1D("hist_polFixFlavour",        " PolFixFlavour         ",50,-1.0,1.0);
  hist_polFixFlavour->Sumw2();
  TH1D *hist_polFixFlavour_ratio   = new TH1D("hist_polFixFlavour_ratio",  " PolFixFlavour_ratio   ",50,-1.0,1.0);
  hist_polFixFlavour_ratio->Sumw2();
 
  //---------------------------------------------------------------------------
  //- Event loop --------------------------------------------------------------
  //---------------------------------------------------------------------------
  while(true) {
    double WT1    = 1.0, WT2    = 1.0;  // we assume that there may be be at most two bosons decaying into
    int    pdgid1 = 0,   pdgid2 = 0;    // taus (tau neutrinos) and requiring reweight

    SimpleParticle B1, B2, X, tau1, tau2;      // SimpleParticle consist of 4 momentum and PDGid.
    vector<SimpleParticle> tau1_daughters, tau2_daughters;

    std::cout << " new event" << std ::endl;

    // Here, we read another event from input_file.
    
    // In this example, we use event generated by Pythia8 with tau decays from
    // Pythia8 or Tauola++.

    // NOTE: for W+/- or H+/- tau2 contain neutrino and tau_daughters2 is empty
    // User may introduce his own method of initializing X, tau, tau2, tau_daughters, tau_daughters2.
    // Then it may be convenient to follow format of the
    // 'read_HepMC' example line below
    int status = read_HepMCExtended(input_file, B1, B2, X, tau1, tau2, tau1_daughters, tau2_daughters);
    //   int status = read_HepMC(input_file, X, tau1, tau2, tau1_daughters, tau2_daughters2);

    // Go to next one if there is nothing more to do with this event
    if( status==1 ) continue;
    
    // Finish if there is nothing more to read from the file
    if( status==0 ) break;
    
    pdgid1 = X.pdgid();

    // Sets polarization state of the input sample
    //setSpinOfSample(0);
        
    // Calculate weight for first boson
    if( abs(X.pdgid())==24 ||  abs(X.pdgid())==37 )
    {
      WT1 = calculateWeightFromParticlesWorHpn(X, tau1, tau2, tau1_daughters); // note that tau2 is tau neutrino
    }
    else if( X.pdgid()==25 || X.pdgid()==36 || X.pdgid()==22 || X.pdgid()==23 )
    {
      // NOTE: if any component of transverse density matrix for DY has been turned on
      //       using setZgamMultipliersTR, weight WT will contain transverse part
      //       of the spin amplitude (calculated from tables table1-1.txt table2-2.txt
      //       which must be stored in the same directory as the program)
      WT1 = calculateWeightFromParticlesH(X, tau1, tau2, tau1_daughters,tau2_daughters);
    }
    else
    {
      cout<<"WARNING: Unexpected PDG for tau mother: "<<X.pdgid()<<endl;
    }
    
    // Returns spin of the taus attributed during reweighting
    double pol = getTauSpin();
 
    //here select first resonance in the event
    double WT = WT1;
    
    wt_average += WT;
    events_count++;

    hist100->Fill(WT);
    hist103->Fill(pol,1.0);
    hist203->Fill(pol,WT);
    hist303->Fill(pol,1/WT);

    //===================================================================================================    
    SimpleParticle         sp_taup;
    SimpleParticle         sp_taum;
    vector<SimpleParticle> sp_taup_daughters;
    vector<SimpleParticle> sp_taum_daughters;
    
    // First we calculate HH for tau+  
    // We enforce that sp_tau is tau+ so the 'nu_tau' is tau-
    if (tau1.pdgid() == -15 )
      {
	sp_taup           = tau1;
	sp_taum           = tau2;
	sp_taup_daughters = tau1_daughters;
	sp_taum_daughters = tau2_daughters;
      }
    else
      {
	sp_taup           = tau2;
	sp_taum           = tau1;
	sp_taup_daughters = tau2_daughters;
	sp_taum_daughters = tau1_daughters;
      }
    
    double *HHFixp, *HHFixm;

    Particle taup(sp_taup.px(), sp_taup.py(), sp_taup.pz(), sp_taup.e(), sp_taup.pdgid() );
    Particle taum(sp_taum.px(), sp_taum.py(), sp_taum.pz(), sp_taum.e(), sp_taum.pdgid() );


    vector<Particle> taup_daughters;
    // Create list of tau daughters
    for(unsigned int i=0; i<sp_taup_daughters.size(); i++)
    {
      Particle pp(sp_taup_daughters[i].px(),
                  sp_taup_daughters[i].py(),
                  sp_taup_daughters[i].pz(),
                  sp_taup_daughters[i].e(),
                  sp_taup_daughters[i].pdgid() );

      taup_daughters.push_back(pp);
    }

    double phi2 = 0.0, theta2 = 0.0;
    //  To calculate matrix elements TAUOLA need tau decay products in tau rest-frame: we first boost all products
    //  to tau rest frame (tau other tau  of Z/H decay along z-axis, intermediate step for boost is tau-tau pair rest-frame), 
    //  then rotate to have neutrino from tau decay along z axis; 
    //  calculated for that purpose angles phi2, theta2 are stored for rotation back of HHp
    prepareKinematicForHH   (taup, taum, taup_daughters, &phi2, &theta2);

    //  Identify decay channel and then calculate polarimetric vector HH; calculates also WTamplit
    HHFixp = calculateHH(taup.pdgid(), taup_daughters, phi2, theta2);

    Particle taupp(sp_taup.px(), sp_taup.py(), sp_taup.pz(), sp_taup.e(), sp_taup.pdgid() );
    Particle taumm(sp_taum.px(), sp_taum.py(), sp_taum.pz(), sp_taum.e(), sp_taum.pdgid() );


    vector<Particle> taum_daughters;
    // Create list of tau daughters
    for(unsigned int i=0; i<sp_taum_daughters.size(); i++)
    {
      Particle pp(sp_taum_daughters[i].px(),
                  sp_taum_daughters[i].py(),
                  sp_taum_daughters[i].pz(),
                  sp_taum_daughters[i].e(),
                  sp_taum_daughters[i].pdgid() );

      taum_daughters.push_back(pp);
    }
    

    phi2 = 0.0, theta2 = 0.0;
    prepareKinematicForHH   (taumm, taupp, taum_daughters, &phi2, &theta2);

    //  Identify decay channel and then calculate polarimetric vector HH; calculates also WTamplit
    HHFixm = calculateHH(taumm.pdgid(), taum_daughters, phi2, theta2);
 
    double polFixFlavour =  getZPolFixFlavour(X, B1, B2, sp_taup, sp_taum);
    double cosThetaStar  =  getCosThetaStar(B1, B2, sp_taup, sp_taum);

    hist_costhetaStar->Fill(cosThetaStar);
    hist_polFixFlavour->Fill(cosThetaStar,polFixFlavour);
    hist_polFixFlavour_ratio->Fill(cosThetaStar,polFixFlavour);

    //    std::cout << "ERW:cosThetaStar  = " << cosThetaStar << std::endl;
    //    std::cout << "ERW:Hpp, Hmm      = " << HHFixp[2] << "  " <<  HHFixm[2] << std::endl;
    //    std::cout << "ERW:FixFlavourPol = " << polFixFlavour << std::endl;
    // we separate cross section into helicity parts.  From this, we attribute helicity states to taus: ++ and --
    double sign = 1.0; 
    double RRR = Tauola::randomDouble();  
    double Polari=1.0;
    if (RRR<(1.0+polFixFlavour)*(1.0+sign*HHFixp[2]*HHFixm[2]+HHFixp[2]+HHFixm[2])/(2.0+2.0*sign*HHFixp[2]*HHFixm[2]+2.0*polFixFlavour*HHFixp[2]+2.0*polFixFlavour*HHFixm[2])) 
      Polari=-1.0; 
    //===================================================================================================    

    hist104->Fill(Polari,1.0);
 

    // calculate correction weight for EW setting
    //components of the cos(th*) calculation
    float Lplus, Lminus, Pplus, Pminus;
    Lplus  = tau1.e()+tau1.pz();
    Lminus = tau1.e()-tau1.pz();
    Pplus  = tau2.e()+tau2.pz();
    Pminus = tau2.e()-tau2.pz();

    Particle P_tautau( tau1.px()+tau2.px(), tau1.py()+tau2.py(), tau1.pz()+tau2.pz(), tau1.e()+tau2.e(), 0 );
            
    //does the cos(th*) calculation
    float cosThetaCS;
    cosThetaCS  = (Lplus*Pminus - Lminus*Pplus);
    cosThetaCS *= fabs(P_tautau.pz());
    cosThetaCS /= (P_tautau.recalculated_mass()*P_tautau.pz());
    cosThetaCS /= sqrt(  P_tautau.recalculated_mass()*P_tautau.recalculated_mass() 
			 + ( P_tautau.px()*P_tautau.px() + P_tautau.py()*P_tautau.py() ) );

    Particle P_tau1( tau1.px(), tau1.py(), tau1.pz(), tau1.e(), tau1.pdgid() );
    P_tau1.print();
    std::cout << "cosThetaCS = " << cosThetaCS  << endl;

    double wt =  (1 + cosThetaCS*cosThetaCS + 0.08 * cosThetaCS )
                /(1 + cosThetaCS*cosThetaCS + 0.14 * cosThetaCS);

    hist_wt->Fill(wt);

    hist_costhetaCS->Fill(cosThetaCS);
    hist_costhetaCS_wt->Fill(cosThetaCS,wt);
    hist_costhetaCS_wt_ratio->Fill(cosThetaCS,wt);
 
    hist_costhetaCS_pol->Fill(cosThetaCS,pol);
    hist_costhetaCS_pol_ratio->Fill(cosThetaCS,pol);
    hist_costhetaCS_pol_wt->Fill(cosThetaCS,pol*wt);
    hist_costhetaCS_pol_wt_ratio->Fill(cosThetaCS,pol*wt);

    hist103_wt->Fill(pol,wt);
    hist104_wt->Fill(Polari,wt);

    if(events_count%10000==0) cout<<"EVT: "<<events_count<<endl;
    if(events_limit && events_count>=events_limit) break;

  }// end of loop
 
  wt_average = wt_average/events_count;

  cout<<"No of events processed: "<<events_count<<endl;
  cout<<"WT average for these events: "<<wt_average<<endl;

  char output_filename[500]="";
  strcat (output_filename,input_filename);
  strcat (output_filename,".root");
  TFile outFile(output_filename, "recreate");

  outFile.cd();
  hist100->Write();
  hist103->Write();
  hist103_wt->Write();
  hist_wt->Write();
  hist203->Write();
  hist303->Write();
  hist104->Write();
  hist104_wt->Write();

  hist_costhetaCS->Write();
  hist_costhetaCS_wt->Write();
  hist_costhetaCS_wt_ratio->Divide(hist_costhetaCS);
  hist_costhetaCS_wt_ratio->Write();
  hist_costhetaCS_pol->Write();
  hist_costhetaCS_pol_ratio->Divide(hist_costhetaCS);
  hist_costhetaCS_pol_ratio->Write();
  hist_costhetaCS_pol_wt->Write();
  hist_costhetaCS_pol_wt_ratio->Divide(hist_costhetaCS);
  hist_costhetaCS_pol_wt_ratio->Write();

  hist_costhetaStar->Write();
  hist_polFixFlavour->Write();
  hist_polFixFlavour_ratio->Divide(hist_costhetaStar);
  hist_polFixFlavour_ratio->Write();

  outFile.Close();



}

/*******************************************************************************
 Get tau tau polarization in case of Z/gamma*-> tau+ tau- case, fixed incoming 
flavour
 incoming configurations
 ID: flavour of incoming partons
 S: invariant mass^2 of the bozon
 &sp_tau: first tau
 &sp_nu_tau: second tau  (in this case it is misleading name)
 *******************************************************************************/

double getZPolFixFlavour(SimpleParticle &sp_X, SimpleParticle &sp_Q,  SimpleParticle &sp_antiQ, SimpleParticle &sp_taup, SimpleParticle &sp_taum)
{
  // tau+ and tau- in lab frame
  Particle tau_plus ( sp_taup.px(),    sp_taup.py(),    sp_taup.pz(),    sp_taup.e(),    sp_taup.pdgid() );
  Particle tau_minus( sp_taum.px(), sp_taum.py(), sp_taum.pz(), sp_taum.e(), sp_taum.pdgid() );
  //  Particle P_B1(    sp_B1.px(), sp_B1.py(), sp_B1.pz(), sp_B1.e(), sp_B1.pdgid() );
  //  Particle P_B2(    sp_B2.px(), sp_B2.py(), sp_B2.pz(), sp_B2.e(), sp_B2.pdgid() );

  // P_QQ = sum of tau+ and tau- in lab frame
  Particle P_QQ( tau_plus.px()+tau_minus.px(), tau_plus.py()+tau_minus.py(), tau_plus.pz()+tau_minus.pz(), tau_plus.e()+tau_minus.e(), 0 );
     
  Particle P_B1(0, 0,-1, 1, sp_Q.pdgid());
  Particle P_B2(0, 0, 1, 1, sp_antiQ.pdgid());
 

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
  //  std::cout << "ERW:cosThetaStar" << costhe << std::endl;   
  // Invariant mass^2 of, tau+tau- pair plus photons, system!
  Particle P_X(sp_X.px(), sp_X.py(), sp_X.pz(), sp_X.e(), 0);
  double SVAR  = P_X.recalculated_mass()*P_X.recalculated_mass();

  int ID = fabs(P_B1.pdgid());
  if( ID == 3 || ID == 5 ) ID=1;
  if( ID == 4 || ID == 6 ) ID=2;
  
  int tau_pdgid = 15;
  double pol = 0;
  double polp = plzap2(P_B1.pdgid(),tau_pdgid,SVAR,costhe); 
  std::cout << "P_B1.pdgid() = " << P_B1.pdgid() << "  P_B1.pz() = " << P_B1.pz() << "  sp_Q.pz() = " << sp_Q.pz() << std::endl;
  if( sp_Q.pz() < 0 )
    polp = plzap2(P_B1.pdgid(),tau_pdgid,SVAR, -costhe);
  //  std::cout << " ERW: plzap2: " << ID << " " << tau_pdgid << " " << SVAR << " " << costhe << std::endl;
  //  std::cout << " ERW: polp: " << polp <<  std::endl;
  pol = (2*(1-polp)-1);

  //  std::cout << " ERW: pol: " << pol <<   std::endl;



  return pol;
}


/*******************************************************************************
 Get tau tau polarization in case of Z/gamma*-> tau+ tau- case, fixed incoming flavour
 incoming configurations
 ID: flavour of incoming partons
 S: invariant mass^2 of the bozon
 &sp_tau: first tau
 &sp_nu_tau: second tau  (in this case it is misleading name)
 *******************************************************************************/

double getCosThetaStar(SimpleParticle &sp_B1,  SimpleParticle &sp_B2, SimpleParticle &sp_taup, SimpleParticle &sp_taum)
{
  // tau+ and tau- in lab frame
  Particle tau_plus (    sp_taup.px(),    sp_taup.py(),    sp_taup.pz(),    sp_taup.e(),    sp_taup.pdgid() );
  Particle tau_minus( sp_taum.px(), sp_taum.py(), sp_taum.pz(), sp_taum.e(), sp_taum.pdgid() );
  //  Particle P_B1(    sp_B1.px(), sp_B1.py(), sp_B1.pz(), sp_B1.e(), sp_B1.pdgid() );
  //  Particle P_B2(    sp_B2.px(), sp_B2.py(), sp_B2.pz(), sp_B2.e(), sp_B2.pdgid() );

  // P_QQ = sum of tau+ and tau- in lab frame
  Particle P_QQ( tau_plus.px()+tau_minus.px(), tau_plus.py()+tau_minus.py(), tau_plus.pz()+tau_minus.pz(), tau_plus.e()+tau_minus.e(), 0 );
   
  Particle P_B1(0, 0,-1, 1, 0);
  Particle P_B2(0, 0, 1, 1, 0);
 

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

  return costhe;
}

