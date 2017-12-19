/**
 * @author T. Przedzinski
 * @date 16 Oct 2013
 *
 * @brief Analyses data sample producing plots for the validation of
 *        effects introduced by TauSpinner
 *
 * !!!IMPORTANT!!!
 * See README for a list of assumptions made by this program.
 *
 */
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
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TFile.h"
#include "TPaveStats.h"

#include <iostream>
#include <fstream>
#include <cstdio>
using namespace std;
using namespace Tauolapp;
using namespace TauSpinner;

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

/* Draw original, unweighted plots and ratio plot below them
   Flag isW determines the range for the ratio axis. */
void draw3(TCanvas *c, TH1D *h1, TH1D *h2, TH1D *hdiff, bool isW=0)
{
  c->Clear();
  TPad *c1npad1 = new TPad("c1npad1","",0,0.2,1,1); // create pad for the plots
  TPad *c1npad2 = new TPad("c1napd2","",0,0,1,0.25); // create pad for the ratio

  c->cd();
  c1npad1->Draw(); // draw larger pad in main area
  c1npad1->cd();   // change working area to inside main pad

  double max1=h1->GetMaximum();
  double max2=h2->GetMaximum();
  double max= max1>max2 ? max1:max2;

  h1->GetYaxis()->SetRangeUser(0,max*1.1);

  h1->SetLineColor(2);
  h1->SetStats(0);
  h1->Draw("hist");

  h2->SetLineColor(3);
  h2->SetStats(0);
  h2->Draw("hist same");

  c->cd();        // go back to working area of whole canvas
  c1npad2->Draw(); // draw smaller pad in main canvas
  c1npad2->cd();   // set working area to smaller pad

  max1 = hdiff->GetMinimum();
  max2 = hdiff->GetMaximum();

  hdiff->SetLineColor(kBlack);
  hdiff->SetStats(0); // setup the options for the ratio plot
  hdiff->SetTitle(NULL);
  hdiff->GetYaxis()->SetTitle("Ratio");
  hdiff->GetYaxis()->SetTitleSize(0.12);
  hdiff->GetYaxis()->SetTitleOffset(0.3);
  hdiff->GetYaxis()->SetLabelSize(0.1);
  hdiff->GetXaxis()->SetLabelSize(0.1);
  if(isW) hdiff->GetYaxis()->SetRangeUser(0,2);
  else    hdiff->GetYaxis()->SetRangeUser(0.75,1.25);
  hdiff->GetYaxis()->SetNdivisions(4,kTRUE);
  hdiff->SetTickLength(0.01,"Y");
  hdiff->SetMarkerStyle(7);
  hdiff->DrawCopy("p e"); // draw the ratio plot
  c->cd(); // go back to main canvas
}

/* Struct for storing the configuration options
   parsed from input config file */
struct config
{
  TString input_file;
  TString input_file2;
  TString output_dir;
  long event_limit;
  
  TString LHAPDF_dataset;
  double TauSpinner_CMSENE;
  int    TauSpinner_Ipol;
} conf;

//-----------------------------------------------------------------------------
//- Main ----------------------------------------------------------------------
//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  if(argc<2)
  {
    cout<<"Usage:   "<<argv[0]<<" <config_file_for_2_datafiles>"<<endl;
    exit(0);
  }
  
  //---------------------------------------------------------------------------
  //- Configuration file parsing ----------------------------------------------
  //---------------------------------------------------------------------------
  ifstream in(argv[1]);
  if(!in.is_open())
  {
    cout<<"Config file "<<argv[1]<<" does not exist."<<endl;
    exit(-1);
  }

  int conf_op=0;
  
  cout<<endl;
  cout<<" Tauspinner validation configuration:"<<endl;
  cout<<" ------------------------------------"<<endl;
  while(!in.eof())
  {
    char buf[256];
    in>>buf;
    if(buf[0]=='#')
    {
      in.getline(buf,256);
      continue;
    }
    
    if(conf_op==0) { conf.input_file        = buf;       cout<<" Input file:       "<<conf.input_file       <<endl; }
    if(conf_op==1) { conf.input_file2       = buf;       cout<<" Input file 2:     "<<conf.input_file2      <<endl; }
    if(conf_op==2) { conf.output_dir        = buf;       cout<<" Output directory: "<<conf.output_dir       <<endl; }
    if(conf_op==3) { conf.event_limit       = atoi(buf); cout<<" Events limit:     "<<conf.event_limit      <<endl; }
    if(conf_op==4) { conf.LHAPDF_dataset    = buf;       cout<<" LHAPDF dataset:   "<<conf.LHAPDF_dataset   <<endl; }
    if(conf_op==5) { conf.TauSpinner_CMSENE = atof(buf); cout<<" CMS energy:       "<<conf.TauSpinner_CMSENE<<endl; }
    if(conf_op==6) { conf.TauSpinner_Ipol   = atoi(buf); cout<<" Ipol:             "<<conf.TauSpinner_Ipol  <<endl; }
    conf_op++;
  }
  cout<<" ------------------------------------"<<endl;
  cout<<endl;
  in.close();

  if(conf_op<=7)
  {
    cout<<"Error! Not all options present in config file!"<<endl;
    cout<<"       Are you sure you're using config file for 2-file comparison?"<<endl;
    exit(-1);
  }
  //---------------------------------------------------------------------------
  //- Initialization ----------------------------------------------------------
  //---------------------------------------------------------------------------
  
  // Subdirectory where plots will be stored
  TString plotdir  = conf.output_dir;
  TString prefix   = ""; // e.g. "embedded_"

  // Limit number of processed events
  int events_limit = conf.event_limit;

  // Initialize LHAPDF
  LHAPDF::initPDFSetByName(conf.LHAPDF_dataset.Data());
  
  // Initialize Tauola
  Tauola::initialize();

  // Initialize TauSpinner
  double CMSENE = conf.TauSpinner_CMSENE; // center of mass system energy (old default was 7000.0).
                                          // used in PDF calculation for pp collisions only
  bool Ipp = true;                        // 'true' for pp collisions
                                          // false otherwise (not prepared yet)
  int Ipol   = conf.TauSpinner_Ipol;      // are input samples polarized? default answer was yes that is  1;

  // NOTE: For this analysis, values below should remain 0.
  //       For details see TauSpinner/examples/tau-reweight-test.cxx
  int nonSM2 = 0; // are we using nonSM calculations?
  int nonSMN = 0; // if we are using nonSM calculations we may want corrections 
                  // to shapes only: y/n  (1/0)

  initialize_spinner(Ipp, Ipol, nonSM2, nonSMN,  CMSENE);

  // Open I/O files
  vector<HepMC::IO_GenEvent*> input_files;
  
  HepMC::IO_GenEvent *input_file  = new HepMC::IO_GenEvent(conf.input_file.Data(), std::ios::in);
  HepMC::IO_GenEvent *input_file2 = new HepMC::IO_GenEvent(conf.input_file2.Data(),std::ios::in);
  
  input_files.push_back(input_file);
  input_files.push_back(input_file2);

  TFile rootfile(plotdir + "/" + prefix + "out.root", "RECREATE");
 
  // Copy config file to output directory
  ifstream conf_in (argv[1]);
  ofstream conf_out( (plotdir + "/" + prefix + "input-config-file.txt").Data() );

  conf_out << conf_in.rdbuf();

  conf_out.close();
  conf_in.close();

  //---------------------------------------------------------------------------
  //- Plot definitions --------------------------------------------------------
  //---------------------------------------------------------------------------
  
  // Mass plot definitions
  TH1D *Z_mass = new TH1D("Z","log10(tau tau mass) in Z decays",100,1,5);
  TH1D *W_mass = new TH1D("W","log10(tau tau mass) in W decays",100,1,5);
  
  // Definition of plots for W
  TH1D *h1[3][3][2] = { { { 0 } } };
  for(int i=0;i<3;i++)
  for(int j=0;j<3;j++)
  for(int k=0;k<2;k++)
  {
    TString  title;
    if(k==0) title  = "W-";
    if(k==1) title  = "W+";
    
    if(i==0) title += " | mu, e";
    if(i==1) title += " | pi";
    if(i==2) title += " | rho";
    
    if(j==1) title += " | unweighted";
    if(j==2) title += " | ratio";
    
    TString h_id = "h1_"; h_id +=i; h_id +=j; h_id +=k;

    h1[i][j][k] = new TH1D(h_id,title,10,0,1);
    h1[i][j][k]->SetStats(0);
    h1[i][j][k]->Sumw2();
  }
  
  // Definition of plots for Z
  TH2D *h2[9][3] = { { 0 } };
  TH1D *h2_vis_mass[9][3] = { { 0 } };
  for(int i=0;i<9;i++)
  for(int j=0;j<3;j++)
  {
    TString  title;
    if(i==0) title  = "mu-, e- vs mu+, e+";
    if(i==1) title  = "mu-, e- vs pi+";
    if(i==2) title  = "pi- vs mu+, e+";
    if(i==3) title  = "pi- vs pi+";
    if(i==4) title  = "mu-, e- vs rho+";
    if(i==5) title  = "rho- vs mu+, e+";
    if(i==6) title  = "pi- vs rho+";
    if(i==7) title  = "rho- vs pi+";
    if(i==8) title  = "rho- vs rho+";
    
    if(j==1) title += " | unweighted";
    if(j==2) title += " | ratio";
    
    TString h_id  = "h2_"; h_id +=i; h_id +="0"; h_id +=j;
   
    h2[i][j] = new TH2D(h_id, title,10,0,1,10,0,1);
    h2[i][j]->SetStats(0);
    h2[i][j]->Sumw2();
        
    h_id += "_vis_mass";
    
    h2_vis_mass[i][j] = new TH1D(h_id, title,24,0,120);
    h2_vis_mass[i][j]->SetStats(0);
    h2_vis_mass[i][j]->Sumw2();
  }

  double wt_average   = 0.0;
  int    events_count = 0;

  //---------------------------------------------------------------------------
  //- Event loops -------------------------------------------------------------
  //---------------------------------------------------------------------------
  for(int input_file_no=0;input_file_no<2;input_file_no++)
  {
    if(input_file_no==0) cout<<"Processing first datafile..."<<endl;
    if(input_file_no==1) cout<<"Processing second datafile..."<<endl;
    
    wt_average   = 0.0;
    events_count = 0;
    
    while(true)
    {
      SimpleParticle X, tau, tau2; // SimpleParticle consist of 4 momentum and PDGid.
      vector<SimpleParticle> tau_daughters, tau_daughters2;

      // Here, we read another event from input_file.

      // NOTE: for W+/- or H+/- tau2 contains neutrino and tau_daughters2 is empty
      // User may introduce his own method of initializing X, tau, tau2, tau_daughters, tau_daughters2.
      // Then it may be convenient to follow format of the
      // 'readParticlesFromTAUOLA_HepMC' example line below
      int status = readParticlesFromTAUOLA_HepMC(*input_files[input_file_no], X, tau, tau2, tau_daughters, tau_daughters2);

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
      if(channel1>0)
      {
        // Tau pair mass
        Particle tau_pair(tau.px()+tau2.px(),tau.py()+tau2.py(),tau.pz()+tau2.pz(),tau.e()+tau2.e(),100);
        double tau_pair_mass = log10(tau_pair.recalculated_mass());
        
        // Visible tau pair mass
        double v_px=0.0,v_py=0.0,v_pz=0.0,v_e=0.0;
        
        for(unsigned int j=0;j<tau_daughters.size();j++)
        {
          int pdgid = abs(tau_daughters[j].pdgid());

          // CONTROVERSIAL LINE (NAMING INCONSISTENCY):
          if(pdgid==22 ) continue; // In principle photons are visible
                                   // and should be added to visible tau mass

          if(pdgid==12 || pdgid==14 || pdgid==16) continue;
          
          v_px += tau_daughters[j].px();
          v_py += tau_daughters[j].py();
          v_pz += tau_daughters[j].pz();
          v_e  += tau_daughters[j].e();
        }
        
        for(unsigned int j=0;j<tau_daughters2.size();j++)
        {
          int pdgid = abs(tau_daughters2[j].pdgid());

          // CONTROVERSIAL LINE (NAMING INCONSISTENCY):
          if(pdgid==22 ) continue; // In principle photons are visible
                                   // and should be added to visible tau mass

          if(pdgid==12 || pdgid==14 || pdgid==16) continue;
          
          v_px += tau_daughters2[j].px();
          v_py += tau_daughters2[j].py();
          v_pz += tau_daughters2[j].pz();
          v_e  += tau_daughters2[j].e();
        }
        
        double vis_tau_pair_mass = sqrt(v_e*v_e - v_px*v_px - v_py*v_py - v_pz*v_pz);
        
        // Energy spectra of charged pion/lepton
        double e1         = tau_daughters[1].e()/tau.e();
        
        // For rho - if no rho is explicitly written in the event record,
        // sum of the energy of pi-/+ and pi0 is used
        if(channel1==4 && tau_daughters.size()==3) e1 = (tau_daughters[1].e() + tau_daughters[2].e())/tau.e();

        // Z/H - 2D scatter plots
        if(X.pdgid()==23 && channel2>0)
        {
          // Fill tau pair mass plot
          if(input_file_no==0) Z_mass->Fill(tau_pair_mass);

          // Energy spectra of charged pion/lepton
          double e2         = tau_daughters2[1].e()/tau2.e();
          
          // For rho - if no rho is explicitly written in the event record,
          // sum of the energy of pi-/+ and pi0 is used
          if(channel2==4 && tau_daughters2.size()==3) e2 = (tau_daughters2[1].e() + tau_daughters2[2].e())/tau2.e();
          
          // Select plot id based on channels combination
          int plot_id=-1;

          if(channel1<3)
          {
            if(channel2<3)  plot_id = 0; // mu,e | mu,e
            if(channel2==3) plot_id = 1; // mu,e | pi
            if(channel2==4) plot_id = 4; // mu,e | rho
          }
          else if(channel1==3)
          {
            if(channel2<3)  plot_id = 2; // pi | mu,e
            if(channel2==3) plot_id = 3; // pi | pi
            if(channel2==4) plot_id = 6; // pi | rho
          }
          else if(channel1==4)
          {
            if(channel2<3)  plot_id = 5; // rho | mu,e
            if(channel2==3) plot_id = 7; // rho | pi
            if(channel2==4) plot_id = 8; // rho | rho
          }
          
          // Plot
          if(input_file_no==0)
          {
            h2[plot_id][0]->Fill(e1,e2,BASE_WT);
            // 3rd plot is ratio plot; will by divided by h2[plot_id][1] later
            h2[plot_id][2]->Fill(e1,e2,BASE_WT);
            
            h2_vis_mass[plot_id][0]->Fill(vis_tau_pair_mass,BASE_WT);
            // 3rd plot is ratio plot; will by divided by h2_vis_mass[plot_id][1] later
            h2_vis_mass[plot_id][2]->Fill(vis_tau_pair_mass,BASE_WT);
          }
          else if(input_file_no==1)
          {
            if(conf.TauSpinner_Ipol==1) h2[plot_id][1]->Fill(e1,e2,BASE_WT/WT);
            else                        h2[plot_id][1]->Fill(e1,e2,BASE_WT*WT);
            
            if(conf.TauSpinner_Ipol==1) h2_vis_mass[plot_id][1]->Fill(vis_tau_pair_mass,BASE_WT/WT);
            else                        h2_vis_mass[plot_id][1]->Fill(vis_tau_pair_mass,BASE_WT*WT);
          }
        }
        // W - 1D plots
        else if(abs(X.pdgid())==24)
        {
          // Fill tau pair mass plot
          if(input_file_no==0) W_mass->Fill(tau_pair_mass);
          
          // Select plot id
          int plot_id=-1;
          if     (channel1<3 ) plot_id = 0; // mu,e
          else if(channel1==3) plot_id = 1; // pi
          else if(channel1==4) plot_id = 2; // rho

          int k=0;
          if(X.pdgid()==24) k=1;
          
          if(input_file_no==0)
          {
            h1[plot_id][0][k]->Fill(e1,BASE_WT);
            
            // 3rd plot is ratio plot; will by divided by  h1[plot_id][1][k] later
            h1[plot_id][2][k]->Fill(e1,BASE_WT);
          }
          else if(input_file_no==1)
          {
            if(conf.TauSpinner_Ipol==1) h1[plot_id][1][k]->Fill(e1,BASE_WT/WT);
            else                        h1[plot_id][1][k]->Fill(e1,BASE_WT*WT);
          }
        }
      }

      if(events_count%10000==0) cout<<"EVT: "<<events_count<<endl;
      if(events_limit && events_count>=events_limit) break;
      
      events_count++;
      wt_average += WT;
    }
  }
  
  // Save plots
  TCanvas c("c","c",800,600);
  c.SetFillColor(kWhite);
  //gPad->SetPhi(210); // Rotate scatter plot
  //gPad->Update();

  // W
  if(W_mass->GetEntries()>0)
  {
    W_mass->Draw();
    TPaveStats *st = (TPaveStats*)W_mass->FindObject("stats");
    if(st) st->SetFillColor(kWhite);
    c.SaveAs(plotdir + "/" + prefix + "W_mass.eps");
    
    for(int i=0;i<3;i++)
    {
      for(int k=0;k<2;k++)
      {
        h1[i][0][k]->Scale(1./h1[i][0][k]->Integral());
        h1[i][1][k]->Scale(1./h1[i][1][k]->Integral());
        h1[i][2][k]->Scale(1./h1[i][2][k]->Integral());
        
        h1[i][2][k]->Divide(h1[i][1][k]);
        
        draw3(&c,h1[i][0][k],h1[i][1][k],h1[i][2][k],1);
        
        TString plotname = h1[i][2][k]->GetName();
        c.SaveAs(plotdir + "/" + prefix + plotname + ".eps");
      }
    }
    
    // Print summary and save it to output file
    ofstream out( (plotdir + "/" + prefix + "input-W-event-count.txt").Data() );
    
    char buf[128];
    cout<<"            W Events:"<<endl;
    out <<"            W Events:"<<endl;
    
    sprintf(buf,"%10i  Total\n",(int)W_mass->GetEntries());
    cout<<buf;
    out <<buf;
    
    for(int k=0;k<2;k++)
    for(int i=0;i<3;i++)
    {
      sprintf(buf,"%10i  %s\n",(int)h1[i][0][k]->GetEntries(),h1[i][0][k]->GetTitle());
      cout<<buf;
      out <<buf;
    }
    
    out.close();
  }

  if(Z_mass->GetEntries()>0)
  {
    Z_mass->Draw();
    TPaveStats *st = (TPaveStats*)Z_mass->FindObject("stats");
    if(st) st->SetFillColor(kWhite);
    c.SaveAs(plotdir + "/" + prefix + "Z_mass.eps");
    
    // Z plots
    for(int i=0;i<9;i++)
    {
      for(int j=0;j<3;j++)
      {
        // Add histogram to MC-TESTER booklet
        //Setup::user_histograms->Add(h2[i][j]);
        
        if(j==2)
        {
          h2[i][0]->Scale(1./h2[i][0]->Integral());
          h2[i][1]->Scale(1./h2[i][1]->Integral());
          h2[i][2]->Scale(1./h2[i][2]->Integral());
        
          h2[i][2]->Divide(h2[i][1]);
          
          h2_vis_mass[i][0]->Scale(1./h2_vis_mass[i][0]->Integral());
          h2_vis_mass[i][1]->Scale(1./h2_vis_mass[i][1]->Integral());
          h2_vis_mass[i][2]->Scale(1./h2_vis_mass[i][2]->Integral());
          
          h2_vis_mass[i][2]->Divide(h2_vis_mass[i][1]);
        }
        
        h2[i][j]->GetXaxis()->SetTitle("Ch-");
        h2[i][j]->GetXaxis()->SetTitleSize(0.06);
        h2[i][j]->GetXaxis()->SetTitleOffset(1.1);
        h2[i][j]->GetYaxis()->SetTitle("Ch+");
        h2[i][j]->GetYaxis()->SetTitleSize(0.06);
        h2[i][j]->GetYaxis()->SetTitleOffset(1.25);
        
        h2[i][j]->Draw("LEGO");
        
        TString plotname = h2[i][j]->GetName();
        c.SaveAs(plotdir + "/" + prefix + plotname + ".eps");
      }
       
      TString plotname = h2[i][0]->GetName();
      TString title    = h2[i][0]->GetTitle();
      
      TH1D *px_0 = h2[i][0]->ProjectionX();
      TH1D *px_1 = h2[i][1]->ProjectionX();
      TH1D *px_2 = (TH1D*)h2[i][0]->ProjectionX()->Clone(plotname + "_px_ratio");
      
      px_2->Divide(px_1);
      
      TH1D *py_0 = h2[i][0]->ProjectionY();
      TH1D *py_1 = h2[i][1]->ProjectionY();
      TH1D *py_2 = (TH1D*)h2[i][0]->ProjectionY()->Clone(plotname + "_py_ratio");
      
      py_2->Divide(py_1);

      
      px_2->GetXaxis()->SetTitle("");
      px_2->SetTitle(title + " | Product charge -");
      
      py_2->GetXaxis()->SetTitle("");
      py_2->SetTitle(title + " | Product charge +");
      
      draw3(&c,px_0,px_1,px_2);
      c.SaveAs(plotdir + "/" + prefix + plotname + "_projectionX.eps");

      draw3(&c,py_0,py_1,py_2);
      c.SaveAs(plotdir + "/" + prefix + plotname + "_projectionY.eps");
      
      draw3(&c,h2_vis_mass[i][0],h2_vis_mass[i][1],h2_vis_mass[i][2]);
      c.SaveAs(plotdir + "/" + plotname + "_vis_mass.eps");
    }
    
    // Print summary and save it to output file
    ofstream out( (plotdir + "/" + prefix + "input-Z-event-count.txt").Data() );
    
    char buf[128];
    cout<<"            Z Events:"<<endl;
    out <<"            Z Events:"<<endl;
    
    sprintf(buf,"%10i  Total\n",(int)Z_mass->GetEntries());
    cout<<buf;
    out <<buf;
    
    for(int i=0;i<9;i++)
    {
      sprintf(buf,"%10i  %s\n",(int)h2[i][0]->GetEntries(),h2[i][0]->GetTitle());
      cout<<buf;
      out <<buf;
    }
    
    out.close();
  }

  rootfile.Write();
  rootfile.Close();
  cout<<endl<<"No of events processed for spin weight: "<<events_count<<endl;
  cout<<      "WT average for these events: "<<wt_average/events_count<<endl;
}
