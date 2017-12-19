{

  gROOT->LoadMacro("PlotStyle.C");
  gROOT->LoadMacro("PlotUtils.C");
  SetAtlasStyle();

  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();


  TFile myfile_Z("out.Z.root");


  TH1F *myh1_mvisZ1 = (TH1F*)myfile_Z.Get("h2_000_px");
  TH1F *myh1_mvisZ2 = (TH1F*)myfile_Z.Get("h2_001_px");
  TH1F *myh1_mvisZr = (TH1F*)myfile_Z.Get("h2_000_px_ratio");

  TH1F *myh1_mvisZ1a = (TH1F*)myfile_Z.Get("h2_200_px");
  TH1F *myh1_mvisZ2a = (TH1F*)myfile_Z.Get("h2_201_px");
  TH1F *myh1_mvisZra = (TH1F*)myfile_Z.Get("h2_200_px_ratio");

  TH1F *myh1_mvisZ1b = (TH1F*)myfile_Z.Get("h2_800_px");
  TH1F *myh1_mvisZ2b = (TH1F*)myfile_Z.Get("h2_801_px");
  TH1F *myh1_mvisZrb = (TH1F*)myfile_Z.Get("h2_800_px_ratio");

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
  myh1_mvisZ1->GetXaxis()->SetTitle("x_{1,2}=E(l)/E(#tau)");
  myh1_mvisZ1->GetXaxis()->SetTitleOffset(1.2);
  myh1_mvisZ1->GetXaxis()->SetTitleSize(0.04);
  myh1_mvisZ1->Draw("hist");

  myh1_mvisZ2->SetLineColor(3);
  myh1_mvisZ2->SetMarkerColor(3);
  myh1_mvisZ2->SetMarkerStyle(22);
  myh1_mvisZ2->SetLineStyle(2);
  //myh1_mvisZ2->SetStats(0);
  myh1_mvisZ2->Draw("p e same");
  myh1_mvisZ1->Draw("hist same");

  leg = new TLegend(0.68,0.9,1.,0.8);
  leg->AddEntry(myh1_mvisZ1,"polarized","l");
  leg->AddEntry(myh1_mvisZ2,"unpolarized","lp");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  TLatex *t = new TLatex();
  t->SetNDC();
  t->DrawLatex(0.22,0.35,"P_{pol}^{fit} = -0.142 #pm 0.003");
  t->DrawLatex(0.22,0.27,"P_{unpol}^{fit} = 0.0006 #pm 0.003");
  c1npad1->Update();

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
  myh1_mvisZr->SetTickLength(0.01,"Y");
  myh1_mvisZr->SetMarkerStyle(20);
  myh1_mvisZr->SetMarkerSize(0.1);
  myh1_mvisZr->DrawCopy("p e"); // draw the ratio plot

  float xmin = myh1_mvisZr->GetBinLowEdge(1); 
  float xmax = xmin + myh1_mvisZr->GetBinWidth(1)* myh1_mvisZr->GetNbinsX();
  TLine *line = new TLine(  xmin,1,xmax,1);
  line->Draw("p e same");

  c3->cd(); // go back to main canvas

  c3.SaveAs("h2_000_projectionX_note.eps");


  //

  TCanvas *c4 = new TCanvas("c4","Stack",700,600);
  c4->Clear();
  TPad *c4npad1 = new TPad("c4npad1","",0,0.2,1,1); // create pad for the plots
  TPad *c4npad2 = new TPad("c4napd2","",0,0,1,0.25); // create pad for the ratio

  c4->cd();
  c4npad1->Draw(); // draw larger pad in main area
  c4npad1->cd();   // change working area to inside main pad

  myh1_mvisZ1a->SetLineColor(2);
  myh1_mvisZ1a->SetStats(0);
  myh1_mvisZ1a->GetXaxis()->SetTitle("x_{1,2}=E(#pi^{#pm})/E(#tau)");
  myh1_mvisZ1a->GetXaxis()->SetTitleOffset(1.2);
  myh1_mvisZ1a->GetXaxis()->SetTitleSize(0.04); 
  myh1_mvisZ1a->Draw("hist");

  myh1_mvisZ2a->SetLineColor(3);
  myh1_mvisZ2a->SetMarkerColor(3);
  myh1_mvisZ2a->SetMarkerStyle(22);
  myh1_mvisZ2a->SetLineStyle(2);
  //myh1_mvisZ2a->SetStats(0);
  myh1_mvisZ2a->Draw("p e same");

  leg1 = new TLegend(0.68,0.5,1.,0.6);
  leg1->AddEntry(myh1_mvisZ1a,"polarised","l");
  leg1->AddEntry(myh1_mvisZ2a,"unpolarised","lp");
  leg1->SetTextSize(0.05);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->Draw();

  TLatex *t4 = new TLatex();
  t4->SetNDC();
  t4->DrawLatex(0.56,0.35,"P_{pol}^{fit} = -0.145 #pm 0.002");
  t4->DrawLatex(0.56,0.27,"P_{unpol}^{fit} = 0.0003 #pm 0.002");
  c4npad1->Update();

  c4->cd();        // go back to working area of whole canvas
  c4npad2->Draw(); // draw smaller pad in main canvas
  c4npad2->cd();   // set working area to smaller pad

  myh1_mvisZra->SetLineColor(kBlack);
  myh1_mvisZra->SetStats(0); // setup the options for the ratio plot
  myh1_mvisZra->SetTitle(NULL);
  myh1_mvisZra->GetYaxis()->SetTitle("Ratio");
  myh1_mvisZra->GetYaxis()->SetTitleSize(0.14);
  myh1_mvisZra->GetYaxis()->SetTitleOffset(0.3);
  myh1_mvisZra->GetYaxis()->SetLabelSize(0.1);
  myh1_mvisZra->GetXaxis()->SetLabelSize(0.1);
  myh1_mvisZra->GetYaxis()->SetRangeUser(0.75,1.25);
  myh1_mvisZra->GetYaxis()->SetNdivisions(4,kTRUE);
  myh1_mvisZra->SetTickLength(0.01,"Y");
  myh1_mvisZra->SetMarkerStyle(20);
  myh1_mvisZra->SetMarkerSize(0.1);
  myh1_mvisZra->DrawCopy("p e"); // draw the ratio plot

  float xmin = myh1_mvisZra->GetBinLowEdge(1); 
  float xmax = xmin + myh1_mvisZra->GetBinWidth(1)* myh1_mvisZra->GetNbinsX();
  TLine *line = new TLine(  xmin,1,xmax,1);
  line->Draw("p e same");

  c4->cd(); // go back to main canvas

  c4.SaveAs("h2_300_projectionX_note.eps");

  //

  TCanvas *c5 = new TCanvas("c5","Stack",700,600);
  c5->Clear();
  TPad *c5npad1 = new TPad("c5npad1","",0,0.2,1,1); // create pad for the plots
  TPad *c5npad2 = new TPad("c5napd2","",0,0,1,0.25); // create pad for the ratio

  c5->cd();
  c5npad1->Draw(); // draw larger pad in main area
  c5npad1->cd();   // change working area to inside main pad

  myh1_mvisZ1b->SetLineColor(2);
  myh1_mvisZ1b->SetStats(0);
  myh1_mvisZ1b->GetXaxis()->SetTitle("x_{1,2}=E(#rho)/E(#tau)");
  myh1_mvisZ1b->GetXaxis()->SetTitleOffset(1.2);
  myh1_mvisZ1b->GetXaxis()->SetTitleSize(0.04);
  myh1_mvisZ1b->SetMaximum(0.016);
  myh1_mvisZ1b->Draw("hist");

  myh1_mvisZ2b->SetLineColor(3);
  myh1_mvisZ2b->SetLineStyle(2);
  myh1_mvisZ2b->SetMarkerColor(3);
  myh1_mvisZ2b->SetMarkerStyle(22);
  myh1_mvisZ2b->SetStats(0);
  myh1_mvisZ2b->Draw("p e same");
  myh1_mvisZ1b->Draw("hist same");

  leg2 = new TLegend(0.68,0.5,1.,0.6);
  leg2->AddEntry(myh1_mvisZ1b,"polarized","l");
  leg2->AddEntry(myh1_mvisZ2b,"unpolarized","lp");
  leg2->SetTextSize(0.05);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->Draw();

  c5->cd();        // go back to working area of whole canvas
  c5npad2->Draw(); // draw smaller pad in main canvas
  c5npad2->cd();   // set working area to smaller pad

  myh1_mvisZrb->SetLineColor(kBlack);
  myh1_mvisZrb->SetStats(0); // setup the options for the ratio plot
  myh1_mvisZrb->SetTitle(NULL);
  myh1_mvisZrb->GetYaxis()->SetTitle("Ratio");
  myh1_mvisZrb->GetYaxis()->SetTitleSize(0.14);
  myh1_mvisZrb->GetYaxis()->SetTitleOffset(0.3);
  myh1_mvisZrb->GetYaxis()->SetLabelSize(0.1);
  myh1_mvisZrb->GetXaxis()->SetLabelSize(0.1);
  myh1_mvisZrb->GetYaxis()->SetRangeUser(0.75,1.25);
  myh1_mvisZrb->GetYaxis()->SetNdivisions(4,kTRUE);
  myh1_mvisZrb->SetTickLength(0.01,"Y");
  myh1_mvisZrb->SetMarkerStyle(20);
  myh1_mvisZrb->SetMarkerSize(0.1);
  myh1_mvisZrb->SetMaximum(1.25);
  myh1_mvisZrb->SetMinimum(0.75);
  myh1_mvisZrb->DrawCopy("p e"); // draw the ratio plot

  float xmin = myh1_mvisZr->GetBinLowEdge(1); 
  float xmax = xmin + myh1_mvisZr->GetBinWidth(1)* myh1_mvisZr->GetNbinsX();
  TLine *line = new TLine(  xmin,1,xmax,1);
  line->Draw("p e same");

  c5->cd(); // go back to main canvas

  c5.SaveAs("h2_800_projectionX_note.eps");

}
