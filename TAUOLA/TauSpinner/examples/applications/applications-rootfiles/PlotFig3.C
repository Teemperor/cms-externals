{
  gROOT->LoadMacro("PlotStyle.C");
  gROOT->LoadMacro("PlotUtils.C");
  SetAtlasStyle();

  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();


  TFile f("out.lepton.root");

  TH1D *h[4] = { 0 };
  h[0] = (TH1D*) f.Get("electron0");
  h[1] = (TH1D*) f.Get("electron1");
  h[2] = (TH1D*) f.Get("electron2");
  h[3] = (TH1D*) f.Get("electron3");

  h[0]->Sumw2();
  h[1]->Sumw2();
  h[2]->Sumw2();


  h[0]->Scale(1./h[0]->Integral());
  h[1]->Scale(1./h[1]->Integral());
  h[2]->Scale(1./h[2]->Integral());

  h[0]->Rebin(10);
  h[1]->Rebin(10);
  h[2]->Rebin(10);

  double max0 = h[0]->GetMaximum();
  double max1 = h[1]->GetMaximum();
  double max2 = h[2]->GetMaximum();

  double max = max0;
  if(max1>max) max = max1;
  if(max2>max) max = max2;

  TCanvas *c1 = new TCanvas("c1","Stack",800,600);
  c1->cd();

  h[0]->SetLineColor(kRed);
  h[0]->SetLineStyle(1);
  h[0]->SetMarkerStyle(23);
  h[0]->SetMarkerColor(kRed);
  h[1]->SetLineColor(3);
  h[1]->SetLineStyle(2);
  h[1]->SetMarkerStyle(22);
  h[1]->SetMarkerColor(3);
  h[2]->SetLineColor(kBlue);
  h[2]->SetLineStyle(5);  

  h[2]->SetTitle("");
  h[2]->SetStats(NULL);
  h[2]->GetXaxis()->SetTitle("x_{1,2}=E(e)/E(#tau)");
  h[2]->SetMaximum(0.18);
  h[2]->Draw(" hist");
  h[1]->Draw("hist SAME");
  h[0]->Draw("hist SAME");

  leg = new TLegend(0.20,0.20,0.60,0.40);
  leg->AddEntry(h[1],"unpolarised, no QED brem","l");
  leg->AddEntry(h[0],"polarised, no QED brem","l");
  leg->AddEntry(h[2],"polarised, with QED brem","l"); 
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();


  //h[0]->GetYaxis()->SetRangeUser(0,max*1.1);

  c1->SaveAs("Z0electronQED_note.eps");
  //c.SaveAs("Z0electronQED.png");


  TH1D *h2[4] = { 0 };
  h2[0] = (TH1D*) f.Get("muon0");
  h2[1] = (TH1D*) f.Get("muon1");
  h2[2] = (TH1D*) f.Get("muon2");
  h2[3] = (TH1D*) f.Get("muon3");

  h2[0]->Scale(1./h2[0]->Integral());
  h2[1]->Scale(1./h2[1]->Integral());
  h2[2]->Scale(1./h2[2]->Integral());

  h2[0]->Rebin(10);
  h2[1]->Rebin(10);
  h2[2]->Rebin(10);

  max0 = h2[0]->GetMaximum();
  max1 = h2[1]->GetMaximum();
  max2 = h2[2]->GetMaximum();

  max = max0;
  if(max1>max) max = max1;
  if(max2>max) max = max2;

  TCanvas *c2 = new TCanvas("c2","Stack",800,600);
  c2->cd();

  h2[0]->SetLineColor(kRed);
  h2[1]->SetLineColor(3);
  h2[1]->SetLineStyle(2);
  h2[2]->SetLineColor(kBlue);
  h2[2]->SetLineStyle(5);  

  h2[0]->SetTitle("");
  h2[0]->SetStats(NULL);
  h2[0]->GetXaxis()->SetTitle("x_{1,2}=E(#mu)/E(#tau)");
  h2[0]->SetMaximum(0.18);
  h2[0]->Draw();
  h2[1]->Draw("SAME");
  h2[2]->Draw("SAME");

  leg2 = new TLegend(0.20,0.20,0.60,0.40);
  leg2->AddEntry(h[1],"unpolarised, no QED brem","l");
  leg2->AddEntry(h[0],"polarised, no QED brem","l");
  leg2->AddEntry(h[2],"polarised, with QED brem","l"); 
  leg2->SetTextSize(0.04);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->Draw();


  // h2[0]->GetYaxis()->SetRangeUser(0,max*1.1);

  c2->SaveAs("Z0muonQED_note.eps");
  //c.SaveAs("Z0muonQED.png");
}
