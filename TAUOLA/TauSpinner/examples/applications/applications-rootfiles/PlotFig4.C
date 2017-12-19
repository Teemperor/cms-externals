{
  gROOT->LoadMacro("PlotStyle.C");
  gROOT->LoadMacro("PlotUtils.C");
  SetAtlasStyle();

  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();

  TFile f("out.W.a1.root");

  TH1D *h[4] = { 0 };
  h[0] = (TH1D*) f.Get("a1_0");
  h[1] = (TH1D*) f.Get("a1_1");

  //TCanvas c("c","c",800,600);
  h[0]->Sumw2();
  h[1]->Sumw2();
  h[0]->Scale(1./h[0]->Integral());
  h[1]->Scale(1./h[1]->Integral());

  //h[0]->Rebin(5);
  //h[1]->Rebin(5);

  double max0 = h[0]->GetMaximum();
  double max1 = h[1]->GetMaximum();

  double max = max0;
  if(max1>max) max = max1;

  TCanvas *c1 = new TCanvas("c1","Stack",800,600);
  c1->cd();

  h[1]->SetLineColor(kGreen);
  h[1]->SetMarkerColor(kGreen);
  h[1]->SetMarkerStyle(22);
  h[1]->SetLineStyle(2);
  h[0]->SetLineColor(kRed);
  h[0]->SetLineStyle(1);

  h[0]->SetTitle("");
  h[0]->SetStats(NULL);
  h[0]->GetXaxis()->SetTitle("E_{#pi-} / E_{vis}");
  h[0]->SetMaximum(0.02);
  h[0]->Draw("hist");
  h[1]->Draw("p e SAME");
  h[0]->Draw("hist same");

  leg = new TLegend(0.63,0.9,0.95,0.8);
  leg->AddEntry(h[0],"polarized","l");
  leg->AddEntry(h[1],"unpolarized","lp");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  // h[0]->GetYaxis()->SetRangeUser(0,max*1.1);

  c1->SaveAs("Wa1pion_note.eps");
  //c.SaveAs("Wa1pion.png");

  TFile f2("out.Z.a1.root");

  TH1D *h2[0] = (TH1D*) f2.Get("a1_0");
  TH1D *h2[1] = (TH1D*) f2.Get("a1_1");

  h2[0]->Sumw2();
  h2[0]->Scale(1./h2[0]->Integral());
  h2[1]->Sumw2();
  h2[1]->Scale(1./h2[1]->Integral());

  //h2[0]->Rebin(5);
  //h2[1]->Rebin(5);

  double max0 = h2[0]->GetMaximum();
  double max1 = h2[1]->GetMaximum();

  double max = max0;
  if(max1>max) max = max1;

  TCanvas *c2 = new TCanvas("c2","Stack",800,600);
  c2->cd();

  h2[1]->SetLineColor(kGreen);
  h2[1]->SetMarkerColor(kGreen);
  h2[1]->SetMarkerStyle(22);
  h2[1]->SetLineStyle(2);
  h2[0]->SetLineColor(kRed);
  h2[0]->SetLineStyle(1);

  h2[0]->SetTitle("");
  h2[0]->SetStats(NULL);
  h2[0]->GetXaxis()->SetTitle("E_{#pi-} / E_{vis}");
  h2[0]->SetMaximum(0.02);
  h2[0]->Draw("hist");
  h2[1]->Draw("p e SAME");
  h2[0]->Draw("hist same");

  //h2[0]->GetYaxis()->SetRangeUser(0,max*1.1);

  leg = new TLegend(0.63,0.9,0.95,0.8);
  leg->AddEntry(h2[0],"polarized","l");
  leg->AddEntry(h2[1],"unpolarized","lp");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  c2->SaveAs("Z0a1pion_note.eps");
  //c.SaveAs("Z0a1pion.png");
}
