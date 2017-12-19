/* Draw original, unweighted plots and ratio plot below them */
void draw2(TCanvas *c, TH1D *h1, TH1D *h2)
{
  h1->Scale(h2->Integral()/h1->Integral());
  TH1D *hdiff = (TH1D*)h1->Clone("_d");
  hdiff->Divide(h2);

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
  hdiff->GetYaxis()->SetRangeUser(0.75,1.25);
  hdiff->GetYaxis()->SetNdivisions(3,kTRUE);
  hdiff->SetTickLength(0.01,"Y");
  hdiff->SetMarkerStyle(7);
  hdiff->DrawCopy("p e"); // draw the ratio plot
  c->cd(); // go back to main canvas
}

test_ipol()
{
  TString path = "test-ipol-plots";
  TFile* f[4] = { 0 };

  f[0] = new TFile("ipol.polarized.root ");       // Ipol = 1 (benchmark. Do not use WT plot from this file)
  f[1] = new TFile("ipol.unpolarized.root");      // Ipol = 0
  f[2] = new TFile("ipol.modified.plzap0.root");  // Ipol = 2
  f[3] = new TFile("ipol.modified.costhe0.root"); // Ipol = 3

  TH1D *all[4][2] = { { 0 } };
  TH1D *lt [4][2] = { { 0 } };
  TH1D *gt [4][2] = { { 0 } };

  // Get all available histograms
  for(int i=0;i<4;i++)
  for(int j=0;j<2;j++)
  {
    TString plotname1="pi_e_all";
    TString plotname2="pi_e_lt";
    TString plotname3="pi_e_gt";

    if(j==1)
    {
      plotname1 += "_WT";
      plotname2 += "_WT";
      plotname3 += "_WT";
    }

    all[i][j] = (TH1D*)f[i]->Get(plotname1);
    lt[i][j]  = (TH1D*)f[i]->Get(plotname2);
    gt[i][j]  = (TH1D*)f[i]->Get(plotname3);

    all[i][j]->Sumw2();
    lt[i][j] ->Sumw2();
    gt[i][j] ->Sumw2();

    all[i][j]->SetStats(0);
    lt[i][j] ->SetStats(0);
    gt[i][j] ->SetStats(0);
  }

  TCanvas c("c","c",800,600);

  for(int j=1;j<4;j++)
  {
    TString test  = "ipol 1 test | ";
    if(j==2) test = "ipol 2 test | ";
    if(j==3) test = "ipol 3 test | ";

    TString title  = "samples without weights | All events";
    TString filename = "nowt_all_"; filename += j;
    all[0][0]->SetTitle(test + title);
    draw2(&c,all[0][0],all[j][0]);
    c.SaveAs(path + "/" + filename + ".eps");

    TString title  = "samples without weights | costheta' < 0.0";
    TString filename = "nowt_lt_"; filename += j;
    lt[0][0]->SetTitle(test + title);
    draw2(&c,lt[0][0],lt[j][0]);
    c.SaveAs(path + "/" + filename + ".eps");

    TString title  = "samples without weights | costheta' > 0.0";
    TString filename = "nowt_gt_"; filename += j;
    gt[0][0]->SetTitle(test + title);
    draw2(&c,gt[0][0],gt[j][0]);
    c.SaveAs(path + "/" + filename + ".eps");

    TString title  = "green line is weighted | All events";
    TString filename = "wt_all_"; filename += j;
    all[0][0]->SetTitle(test + title);
    draw2(&c,all[0][0],all[j][1]);
    c.SaveAs(path + "/" + filename + ".eps");

    TString title  = "green line is weighted | costheta' < 0.0";
    TString filename = "wt_lt_"; filename += j;
    lt[0][0]->SetTitle(test + title);
    draw2(&c,lt[0][0],lt[j][1]);
    c.SaveAs(path + "/" + filename + ".eps");

    TString title  = "green line is weighted | costheta' > 0.0";
    TString filename = "wt_gt_"; filename += j;
    gt[0][0]->SetTitle(test + title);
    draw2(&c,gt[0][0],gt[j][1]);
    c.SaveAs(path + "/" + filename + ".eps");
  }
}
