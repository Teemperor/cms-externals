// HepMC IO_GenEvent header
#include "HepMC/IO_GenEvent.h"

// TAUOLA header
#include "Tauola/Tauola.h"

#include "TauSpinner/Particle.h"
#include "TauSpinner/tau_reweight_lib.h"

#include "read_particles_from_TAUOLA.h"

// LHAPDF header
#include "LHAPDF/LHAPDF.h"

// ROOT header
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
using namespace std;
using namespace Tauolapp;

// We take sigborn and PDF distributions from TauSpinner for this test
namespace TauSpinner{
  extern double sigborn(int ID, double SS, double costhe);
  extern double f(double x, int ID, double SS, double cmsene);
}

const double CMSENE  = 8000.;

/** Calculate cosine of the hard scattering angle */
pair<double,double> get_costheta_and_mass(Particle &q, Particle &anti_q, Particle &taum, Particle &taup)
{
  Particle tau_minus(  taum.px(),  taum.py(),  taum.pz(),  taum.e(),  taum.pdgid());
  Particle tau_plus (  taup.px(),  taup.py(),  taup.pz(),  taup.e(),  taup.pdgid());
  Particle b1       (     q.px(),     q.py(),     q.pz(),     q.e(),     q.pdgid());
  Particle b2       (anti_q.px(),anti_q.py(),anti_q.pz(),anti_q.e(),anti_q.pdgid());
  
  // P_QQ = sum of tau+ and tau- in lab frame
  Particle P_QQ( tau_plus.px()+tau_minus.px(), tau_plus.py()+tau_minus.py(), tau_plus.pz()+tau_minus.pz(), tau_plus.e()+tau_minus.e(), 0 );

  tau_minus.boostToRestFrame(P_QQ);
  tau_plus. boostToRestFrame(P_QQ);
  b1.       boostToRestFrame(P_QQ);
  b2.       boostToRestFrame(P_QQ);
  
  double costheta1 = (tau_minus.px()*b1.px()       +tau_minus.py()*b1.py()       +tau_minus.pz()*b1.pz()       ) /
                 sqrt(tau_minus.px()*tau_minus.px()+tau_minus.py()*tau_minus.py()+tau_minus.pz()*tau_minus.pz()) /
                 sqrt(b1.px()       *b1.px()       +b1.py()       *b1.py()       +b1.pz()       *b1.pz()       );

  double costheta2 = (tau_plus.px()*b2.px()      +tau_plus.py()*b2.py()      +tau_plus.pz()*b2.pz()      ) /
                 sqrt(tau_plus.px()*tau_plus.px()+tau_plus.py()*tau_plus.py()+tau_plus.pz()*tau_plus.pz()) /
                 sqrt(b2.px()      *b2.px()      +b2.py()      *b2.py()      +b2.pz()      *b2.pz()      );
               
  double sintheta1 = sqrt(1-costheta1*costheta1);
  double sintheta2 = sqrt(1-costheta2*costheta2);
  
  //return (costheta1*sintheta2 + costheta2*sintheta1) / (sintheta1 + sintheta2);
  return pair<double,double>( (costheta1 + costheta2) / 2. , P_QQ.recalculated_mass() );
}

/** This is part of the 'getLongitudinalPolarization' from TauSpinner,
    returning its intermediate result - WID and costhe */
pair<double,double*> calculate_costhe_and_WID(double S, Particle &taum, Particle &taup)
{
  Particle tau_minus(  taum.px(),  taum.py(),  taum.pz(),  taum.e(),  taum.pdgid());
  Particle tau_plus (  taup.px(),  taup.py(),  taup.pz(),  taup.e(),  taup.pdgid());

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
  
  double x1x2  = SS/CMSENE/CMSENE;
  double x1Mx2 = P_QQ.pz()/CMSENE*2;
  
  double x1 = (  x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;
  double x2 = ( -x1Mx2 + sqrt(x1Mx2*x1Mx2 + 4*x1x2) )/2;
  
  double *WID = new double[11];
  WID[0] = f(x1, 0,SS,CMSENE)*f(x2, 0,SS,CMSENE) * sigborn(0,SS, costhe);
  WID[1] = f(x1, 1,SS,CMSENE)*f(x2,-1,SS,CMSENE) * sigborn(1,SS, costhe);
  WID[2] = f(x1,-1,SS,CMSENE)*f(x2, 1,SS,CMSENE) * sigborn(1,SS,-costhe);
  WID[3] = f(x1, 2,SS,CMSENE)*f(x2,-2,SS,CMSENE) * sigborn(2,SS, costhe);
  WID[4] = f(x1,-2,SS,CMSENE)*f(x2, 2,SS,CMSENE) * sigborn(2,SS,-costhe);
  WID[5] = f(x1, 3,SS,CMSENE)*f(x2,-3,SS,CMSENE) * sigborn(1,SS, costhe);
  WID[6] = f(x1,-3,SS,CMSENE)*f(x2, 3,SS,CMSENE) * sigborn(1,SS,-costhe);
  WID[7] = f(x1, 4,SS,CMSENE)*f(x2,-4,SS,CMSENE) * sigborn(2,SS, costhe);
  WID[8] = f(x1,-4,SS,CMSENE)*f(x2, 4,SS,CMSENE) * sigborn(2,SS,-costhe);
  WID[9] = f(x1, 5,SS,CMSENE)*f(x2,-5,SS,CMSENE) * sigborn(1,SS, costhe);
  WID[10]= f(x1,-5,SS,CMSENE)*f(x2, 5,SS,CMSENE) * sigborn(1,SS,-costhe);
  
  double sum = 0.0;  // normalize
  for(int i=0;i<=10;i++) sum+=WID[i];

  for(int i=0;i<=10;i++) WID[i]/=sum;
  
  return pair<double,double*>(costhe,WID);
}

int main(int argc, char **argv)
{
  if(argc<3)
  {
    cout<<"Usage: "<<argv[0]<<" <input_file> <output_file> [<events_limit>]"<<endl;
    exit(0);
  }

  int events_limit = 0;
  if(argc>=4) events_limit = atoi(argv[3]);
  
  // Initialize Tauola
  Tauola::initialize();

  // Initialize LHAPDF
  string name="cteq6ll.LHpdf";
  //string name="MRST2004FF4lo.LHgrid";    // Pythia8 default 
  //string name="MSTW2008nnlo90cl.LHgrid"; // TauSpinner examples default
                                           // NOTE: improper choice of PDFs can couse problems
                                           //       when using TauSpinner
  LHAPDF::initPDFSetByName(name);
  
  // Open I/O files
  HepMC::IO_GenEvent input_file(argv[1],std::ios::in);

  TFile out(argv[2],"RECREATE");
  
  // Plot definitions  
  TH1D* hUp[2]   = { 0 };
  TH1D* hDown[2] = { 0 };
  
  hUp[0]   = new TH1D("hUp1",  "up quarks | costheta",100,-1,1);
  hUp[1]   = new TH1D("hUp2",  "up quarks | A",       100,-2,2);
  hDown[0] = new TH1D("hDown1","up quarks | costheta",100,-1,1);
  hDown[1] = new TH1D("hDown2","up quarks | A",       100,-2,2);
  
  TH1D* hCosthetaSpread = new TH1D("hCosthetaSpread","difference in costheta calculation",100,-1,1);
  
  int    print_count  = 10;
  int    events_count = 0;
  int    events_count_up   = 0;
  int    events_count_down = 0;
  double averageAup   = 0.0;
  double averageAdown = 0.0;
  double tauSpinnerUpRatio   = 0.0;
  double tauSpinnerDownRatio = 0.0;
  
  // Begin event loop
  while(true) {

    if(events_count%1000==0) cout<<"EVT: "<<events_count<<endl;
    if(events_limit && events_count>=events_limit) break;
    
    HepMC::GenEvent *HepMCEvt = new HepMC::GenEvent();
    input_file.fill_next_event(HepMCEvt);
    if(input_file.rdstate()) break;
    
    events_count++;
    
    // Get information from HepMC event
    HepMC::GenParticle *hq = NULL, *hanti_q = NULL, *htaum = NULL, *htaup = NULL;
    for(HepMC::GenEvent::particle_const_iterator p=HepMCEvt->particles_begin();p!=HepMCEvt->particles_end();p++)
    {
      if(!(*p)->end_vertex()) continue;
      int pdg  = (*p)->pdg_id();
      int pdg2 = (*(*p)->end_vertex()->particles_out_const_begin())->pdg_id();
 
      if(pdg< 7 && pdg>0 && (*p)->end_vertex()->particles_out_size()==1 && pdg2==23 ) hq      = (*p);
      if(pdg>-7 && pdg<0 && (*p)->end_vertex()->particles_out_size()==1 && pdg2==23 ) hanti_q = (*p);
      if(pdg== 15 && (*p)->end_vertex()->particles_out_size()>1 && abs( (*(*p)->end_vertex()->particles_out_const_begin())->pdg_id() )==16) htaum = (*p);
      if(pdg==-15 && (*p)->end_vertex()->particles_out_size()>1 && abs( (*(*p)->end_vertex()->particles_out_const_begin())->pdg_id() )==16) htaup = (*p);
    }
    
    if(!hq || !hanti_q || !htaum || !htaup) continue;
    
    Particle q     (hq->     momentum().px(),hq->     momentum().py(),hq->     momentum().pz(),hq->     momentum().e(),hq->     pdg_id());
    Particle anti_q(hanti_q->momentum().px(),hanti_q->momentum().py(),hanti_q->momentum().pz(),hanti_q->momentum().e(),hanti_q->pdg_id());
    Particle taum  (htaum->  momentum().px(),htaum->  momentum().py(),htaum->  momentum().pz(),htaum->  momentum().e(),htaum->  pdg_id());
    Particle taup  (htaup->  momentum().px(),htaup->  momentum().py(),htaup->  momentum().pz(),htaup->  momentum().e(),htaup->  pdg_id());

    pair<double,double> ret = get_costheta_and_mass(q,anti_q,taum,taup);
    double costheta = ret.first;
    double SS       = ret.second*ret.second;
    
    // We use sigborn from TauSpinner to calculate A
    // NOTE: incoming quark ID definition is different in sigborn
    int    ID2 = 1;
    if(abs(q.pdgid())==2 || abs(q.pdgid())==4 || abs(q.pdgid())==6 ) ID2 = 2;
    
    double A1  = TauSpinner::sigborn(ID2,SS, 0.5);
    double A2  = TauSpinner::sigborn(ID2,SS,-0.5);
    double A   = 2.5 * (A1 - A2)/(A1 + A2);
    
    if(q.pdgid()==1)
    {
      hDown[0]->Fill(costheta);
      hDown[1]->Fill(A);
      events_count_down++;
      averageAdown += A;
    }
    else if(q.pdgid()==2)
    {
      hUp[0]->Fill(costheta);
      hUp[1]->Fill(A);
      events_count_up++;
      averageAup += A;
    }
    
    pair<double,double*> ret2 = calculate_costhe_and_WID(SS,taum,taup);
    
    double costhe = ret2.first;
    double *WID   = ret2.second;
    
    tauSpinnerDownRatio += WID[1] + WID[2];
    tauSpinnerUpRatio   += WID[3] + WID[4];
    
    delete WID;
    
    // difference between costheta and costhe as calculated by TauSpinner
    // NOTE: abs is needed because the difference in sign might be non-physical
    hCosthetaSpread->Fill( abs(costheta)-abs(costhe) );
    
    
    if(events_count<print_count)
    {
      cout<<"EVENT: "<<events_count<<endl;
      cout<<"In lab frame:"<<endl;
      q.print();
      anti_q.print();
      taum.print();
      taup.print();
      
      // P_QQ = sum of tau+ and tau- in lab frame
      Particle P_QQ( taum.px()+taup.px(), taum.py()+taup.py(), taum.pz()+taup.pz(), taum.e()+taup.e(), 0 );

      q.     boostToRestFrame(P_QQ);
      anti_q.boostToRestFrame(P_QQ);
      taum.  boostToRestFrame(P_QQ);
      taup.  boostToRestFrame(P_QQ);
      
      cout<<"In rest frame of tau pair:"<<endl;
      q.print();
      anti_q.print();
      taum.print();
      taup.print();
      cout<<"cos(theta): "<<costheta<<" SS: "<<SS<<" A: "<<A<<endl<<endl;
    }
    
    delete HepMCEvt;
  }

  out.Write();
  
  // Calcualte A_bar from fits
  TF1 fitUp  ("fitUp",  "1 + x*x + [0]*x",-1,1);
  TF1 fitDown("fitDown","1 + x*x + [0]*x",-1,1);
  fitUp.  SetLineColor(kBlue);
  fitUp.  SetParName  (0,"A");
  fitUp.  SetParameter(0,0.0);
  fitDown.SetLineColor(kBlue);
  fitDown.SetParName  (0,"A");
  fitDown.SetParameter(0,0.0);
  
  // Draw plots
  TCanvas c("c","c",800,600);
  
  hCosthetaSpread->Draw();
  
  c.SaveAs("bornAFB-test-CosthetaSpread.eps");
  
  c.Clear();
  c.SetCanvasSize(1600,1200);
  c.Divide(2,2);
  c.cd(1);
  hDown[0]->SetTitle("quarks down costheta");
  hDown[0]->Draw();
  hDown[0]->Fit("fitDown");
  hDown[0]->Scale(fitDown.Integral(-1,1)/(hDown[0]->Integral()*hDown[0]->GetBinWidth(1)));
  
  c.cd(2);
  hDown[1]->SetTitle("quarks down A");
  hDown[1]->Draw();
  
  c.cd(3);
  hUp[0]->SetTitle("quarks up costheta");
  hUp[0]->Draw();
  hUp[0]->Fit("fitUp");
  hUp[0]->Scale(fitUp.Integral(-1,1)/(hUp[0]->Integral()*hUp[0]->GetBinWidth(1)));

  c.cd(4);
  hUp[1]->SetTitle("quarks up A");
  hUp[1]->Draw();
  
  c.SaveAs("bornAFB-test.eps");
    
  cout<<endl<<"No of events processed:     "<<events_count<<endl;
  cout<<      "Down quarks A_bar from fit: "<<fitDown.GetParameter(0)<<" +/- "<<fitDown.GetParError(0)/sqrt(events_count_down)<<endl;
  cout<<      "Down quarks average A:      "<<averageAdown/events_count_down<<endl;
  cout<<      "Up   quarks A_bar from fit: "<<fitUp.GetParameter(0)<<" +/- "<<fitUp.GetParError(0)/sqrt(events_count_up)<<endl;
  cout<<      "Up   quarks average A:      "<<averageAup/events_count_up<<endl;
  cout<<endl;
  
  cout<<"Rate of down quarks: "<<events_count_down*1.0/events_count<<" from TauSpinner: "<<tauSpinnerDownRatio/events_count<<endl;
  cout<<"Rate of up quarks:   "<<events_count_up  *1.0/events_count<<" from TauSpinner: "<<tauSpinnerUpRatio/events_count<<endl;
}
