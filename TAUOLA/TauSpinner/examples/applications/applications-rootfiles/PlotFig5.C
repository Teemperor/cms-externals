{
 gROOT->LoadMacro("PlotStyle.C");
 gROOT->LoadMacro("PlotUtils.C");
  SetAtlasStyle();

 gROOT->SetStyle("ATLAS");
 gROOT->ForceStyle();


 TFile myfile("out.test.bornAFB.root");


 TH1F *my_hDown1 = (TH1F*)myfile.Get("hDown1");
 TH1F *my_hUp1   = (TH1F*)myfile.Get("hUp1");

 TCanvas *c1 = new TCanvas("c1","Stack",800,600);
 c1->cd();

 TF1 *f1 = new TF1("myfunc","1+x*x+1.689*x",-1,1);

 my_hDown1->SetTitle("");
 my_hDown1->SetStats(NULL);
 my_hDown1->GetXaxis()->SetTitle("cos#theta");
 my_hDown1->Sumw2();
 my_hDown1->Scale(f1->Integral(-1,1)/(my_hDown1->Integral()*my_hDown1->GetBinWidth(1)));
 my_hDown1->SetMarkerStyle(20);
 my_hDown1->SetMarkerSize(0.8);
 my_hDown1->Draw("p e");

 TF1 *f1 = new TF1("myfunc","1+x*x+1.692*x",-1,1);
 f1->SetLineColor(2);
 f1->Draw("SAME");
 my_hDown1->Draw("p e same");

 TLatex *t = new TLatex();
 t->SetNDC();
 t->DrawLatex(0.22,0.65,"A^{fit}= 1.692#pm0.003");

 c1->Update();
 c1->SaveAs("quarkdown.eps");

 TCanvas *c2 = new TCanvas("c2","Stack",800,600);
 c2->cd();

 TF1 *f2 = new TF1("myfunc","1+x*x+1.617*x",-1,1);

 my_hUp1->SetTitle("");
 my_hUp1->SetStats(NULL);
 my_hUp1->GetXaxis()->SetTitle("cos#theta");
 my_hUp1->Sumw2();
 my_hUp1->Scale(f2->Integral(-1,1)/(my_hUp1->Integral()*my_hUp1->GetBinWidth(1)));;
 my_hUp1->SetMarkerStyle(20);
 my_hUp1->SetMarkerSize(0.8);
 my_hUp1->Draw("p e");


 f2->SetLineColor(2);
 f2->Draw("SAME");
 my_hUp1->Draw("p e same");

 TLatex *t1 = new TLatex();
 t1->SetNDC();
 t1->DrawLatex(0.22,0.65,"A^{fit}= 1.617#pm0.002");

 c2->Update();
 c2->SaveAs("quarkup.eps");
}
