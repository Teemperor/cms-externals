{

  gROOT->LoadMacro("PlotStyle.C");
  gROOT->LoadMacro("PlotUtils.C");
  SetAtlasStyle();

  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();


  TFile myfile_Z("out.Z.root");
  TFile myfile_H("out.H.root");

  TH2F *myh2_Z = (TH2F*)myfile_Z.Get("h2_300");
  TH2F *myh2_H = (TH2F*)myfile_H.Get("h3_300");

  TH1F *myh1_mvisZ1 = (TH1F*)myfile_Z.Get("h2_300_vis_mass");
  TH1F *myh1_mvisZ2 = (TH1F*)myfile_Z.Get("h2_301_vis_mass");
  TH1F *myh1_mvisZr = (TH1F*)myfile_Z.Get("h2_302_vis_mass");

  TH1F *myh1_mvisH1 = (TH1F*)myfile_H.Get("h3_300_vis_mass");
  TH1F *myh1_mvisH2 = (TH1F*)myfile_H.Get("h3_301_vis_mass");
  TH1F *myh1_mvisHr = (TH1F*)myfile_H.Get("h3_302_vis_mass");



  //        
          
  TCanvas *c1 = new TCanvas("c1","Stack",700,600);
  c1->cd();

  myh2_Z->GetXaxis()->SetTitle("x_{2}=E_{#pi^{-}}/E_{#tau^{-}}");
  myh2_Z->GetXaxis()->SetTitleSize(0.06);
  myh2_Z->GetXaxis()->SetTitleOffset(1.1);
  myh2_Z->GetXaxis()->SetLabelSize(0.03);
  myh2_Z->GetYaxis()->SetTitle("x_{1}=E_{#pi^{+}}/E_{#tau^{+}}");
  myh2_Z->GetYaxis()->SetTitleSize(0.06);
  myh2_Z->GetYaxis()->SetTitleOffset(1.25);
  myh2_Z->GetYaxis()->SetLabelSize(0.03);
  myh2_Z->GetZaxis()->SetLabelSize(0.03);
  myh2_Z->Draw("LEGO");
    
  c1.SaveAs("h2_300_note.eps");

  //

  TCanvas *c2 = new TCanvas("c2","Stack",700,600);
  c2->cd();
  myh2_H->GetXaxis()->SetTitle("x_{2}=E_{#pi^{-}}/E_{#tau^{-}}");
  myh2_H->GetXaxis()->SetTitleSize(0.06);
  myh2_H->GetXaxis()->SetTitleOffset(1.1);
  myh2_H->GetXaxis()->SetLabelSize(0.03);
  myh2_H->GetYaxis()->SetTitle("x_{1}=E_{#pi^{+}}/E_{#tau^{+}}");
  myh2_H->GetYaxis()->SetTitleSize(0.06);
  myh2_H->GetYaxis()->SetTitleOffset(1.25);
  myh2_H->GetYaxis()->SetLabelSize(0.03);
  myh2_H->GetZaxis()->SetLabelSize(0.03);
  myh2_H->Draw("LEGO");

  c2.SaveAs("h3_300_note.eps");

  //

  TCanvas *c3 = new TCanvas("c3","Stack",700,600);
  c3->Clear();
  TPad *c1npad1 = new TPad("c1npad1","",0,0.2,1,1); // create pad for the plots
  TPad *c1npad2 = new TPad("c1napd2","",0,0,1,0.25); // create pad for the ratio

  c3->cd();
  c1npad1->Draw(); // draw larger pad in main area
  c1npad1->cd();   // change working area to inside main pad

  myh1_mvisZ1->SetLineColor(2);
  myh1_mvisZ1->SetStats(0);
  myh1_mvisZ1->GetXaxis()->SetTitle("M_{vis} [GeV]");
  myh1_mvisZ1->GetXaxis()->SetTitleOffset(1.2);
  myh1_mvisZ1->GetXaxis()->SetTitleSize(0.04);
  myh1_mvisZ1->SetMaximum(0.021);
  myh1_mvisZ1->Draw("hist");

  myh1_mvisZ2->SetLineColor(3);
  myh1_mvisZ2->SetLineStyle(2);
  myh1_mvisZ2->SetMarkerStyle(23);
  myh1_mvisZ2->SetMarkerColor(3);
  //myh1_mvisH2->SetStats(0);
  myh1_mvisZ2->Draw("p e same");
  myh1_mvisZ1->Draw("hist same");

  leg = new TLegend(0.63,0.9,0.95,0.8);
  leg->AddEntry(myh1_mvisZ1,"polarized","l");
  leg->AddEntry(myh1_mvisZ2,"unpolarized","lp");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  c3->cd();        // go back to working area of whole canvas
  c1npad2->Draw(); // draw smaller pad in main canvas
  c1npad2->cd();   // set working area to smaller pad

  myh1_mvisZr->SetLineColor(kBlack);
  myh1_mvisZr->SetStats(0); // setup the options for the ratio plot
  myh1_mvisZr->SetTitle(NULL);
  myh1_mvisZr->GetYaxis()->SetTitle("Ratio");
  myh1_mvisZr->GetYaxis()->SetTitleSize(0.14);
  myh1_mvisZr->GetYaxis()->SetTitleOffset(0.3);
  myh1_mvisZr->GetYaxis()->SetLabelSize(0.1);
  myh1_mvisZr->GetXaxis()->SetLabelSize(0.1);
  myh1_mvisZr->GetYaxis()->SetRangeUser(0.75,1.25);
  myh1_mvisZr->GetYaxis()->SetNdivisions(4,kTRUE);
  myh1_mvisZr->SetMaximum(3.1);
  myh1_mvisZr->SetMinimum(-0.4);
  myh1_mvisZr->SetTickLength(0.01,"Y");
  myh1_mvisZr->SetMarkerStyle(20);
  myh1_mvisZr->SetMarkerSize(0.1);
  myh1_mvisZr->DrawCopy("p e"); // draw the ratio plot

  float xmin = myh1_mvisZr->GetBinLowEdge(1); 
  float xmax = xmin + myh1_mvisZr->GetBinWidth(1)* myh1_mvisZr->GetNbinsX();
  TLine *line = new TLine(  xmin,1,xmax,1);
  line->Draw("p e same");

  c3->cd(); // go back to main canvas

  c3.SaveAs("h2_300_vis_mass_note.eps");

  //

  TCanvas *c4 = new TCanvas("c4","Stack",700,600);
  c4->Clear();
  TPad *c4npad1 = new TPad("c4npad1","",0,0.2,1,1); // create pad for the plots
  TPad *c4npad2 = new TPad("c4napd2","",0,0,1,0.25); // create pad for the ratio

  c4->cd();
  c4npad1->Draw(); // draw larger pad in main area
  c4npad1->cd();   // change working area to inside main pad

  myh1_mvisH1->SetLineColor(2);
  myh1_mvisH1->SetStats(0);
  //myh1_mvisH1->SetMaximum(0.11);
  myh1_mvisH1->GetXaxis()->SetTitle("M_{vis} [GeV]");
  myh1_mvisH1->GetXaxis()->SetTitleOffset(1.2);
  myh1_mvisH1->GetXaxis()->SetTitleSize(0.04);
  myh1_mvisH1->SetMaximum(0.021);
  myh1_mvisH1->Draw("hist");

  myh1_mvisH2->SetLineColor(3);
  myh1_mvisH2->SetLineStyle(2);
  //myh1_mvisH2->SetStats(0);
  myh1_mvisH2->SetMarkerStyle(23);
  myh1_mvisH2->SetMarkerColor(3);
  myh1_mvisH2->Draw("p e same");
  myh1_mvisH1->Draw("hist same");

  leg1 = new TLegend(0.63,0.9,0.95,0.8);
  leg1->AddEntry(myh1_mvisH1,"polarized","l");
  leg1->AddEntry(myh1_mvisH2,"unpolarized"," pl");
  leg1->SetTextSize(0.05);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->Draw();


  c4->cd();        // go back to working area of whole canvas
  c4npad2->Draw(); // draw smaller pad in main canvas
  c4npad2->cd();   // set working area to smaller pad

  myh1_mvisHr->SetLineColor(kBlack);
  myh1_mvisHr->SetStats(0); // setup the options for the ratio plot
  myh1_mvisHr->SetTitle(NULL);
  myh1_mvisHr->GetYaxis()->SetTitle("Ratio");
  myh1_mvisHr->GetYaxis()->SetTitleSize(0.14);
  myh1_mvisHr->GetYaxis()->SetTitleOffset(0.3);
  myh1_mvisHr->GetYaxis()->SetLabelSize(0.1);
  myh1_mvisHr->GetXaxis()->SetLabelSize(0.1);
  myh1_mvisHr->GetYaxis()->SetRangeUser(0.75,1.25);
  myh1_mvisHr->GetYaxis()->SetNdivisions(4,kTRUE);
  myh1_mvisHr->SetMaximum(3.1);
  myh1_mvisHr->SetMinimum(-0.4);
  myh1_mvisHr->SetTickLength(0.01,"Y");
  myh1_mvisHr->SetMarkerStyle(20);
  myh1_mvisHr->SetMarkerSize(0.1);
  myh1_mvisHr->DrawCopy("p e"); // draw the ratio plot

  float xmin = myh1_mvisHr->GetBinLowEdge(1); 
  float xmax = xmin + myh1_mvisHr->GetBinWidth(1)* myh1_mvisHr->GetNbinsX();
  TLine *line = new TLine(  xmin,1,xmax,1);
  line->Draw("p e same");

  c4->cd(); // go back to main canvas
      
  c4.SaveAs("h3_300_vis_mass_note.eps");    

}
