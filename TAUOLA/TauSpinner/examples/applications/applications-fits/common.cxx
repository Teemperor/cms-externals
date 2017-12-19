#include "common.h"

// ROOT headers
#include "TF1.h"
#include "TPad.h"
#include "TPaveText.h"

//RooFit headers
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"

/* Draw original, unweighted plots and ratio plot below them
   Flag isW determines the range for the ratio axis. */
void draw3(TCanvas *c, TH1D *h1, TH1D *h2, TH1D *hdiff, bool isW, TString fit, Double_t x0, Double_t x1) {
    c->Clear();

    TPad *c1npad1, *c1npad2, *c1npad3 = NULL;
    if (fit == "") {
        c1npad1 = new TPad("c1npad1", "", 0, 0.2, 1, 1, kWhite); // create pad for the plots
        c1npad2 = new TPad("c1napd2", "", 0, 0, 1, 0.25, kWhite); // create pad for the ratio
    } else {
        c1npad1 = new TPad("c1npad1", "", 0, 0.4, 1, 1, kWhite); // create pad for the plots
        c1npad2 = new TPad("c1napd2", "", 0, 0.2, 1, 0.4, kWhite); // create pad for the ratio
        c1npad3 = new TPad("c1napd3", "", 0, 0, 1, 0.2, kWhite);
    }

    c->cd();
    c1npad1->SetFillColor(kWhite);
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

    RooRealVar y1("y1", "y1", x0, x1);
    RooRealVar P1("P1", "P1", -2, 2);
    TString fit1 = fit;
    fit1 = fit1.ReplaceAll("x", "y1");
    fit1 = fit1.ReplaceAll("P", "P1");
    RooDataHist data1("data1", "data1", y1, h1);

    RooRealVar y2("y2", "y2", x0, x1);
    RooRealVar P2("P2", "P2", -2, 2);
    TString fit2 = fit;
    fit2 = fit2.ReplaceAll("x", "y2");
    fit2 = fit2.ReplaceAll("P", "P2");
    RooDataHist data2("data2", "data2", y2, h2);

    if (fit != "") {
        RooGenericPdf fit_pdf1("fit_pdf1", fit1, RooArgSet(y1, P1));
        fit_pdf1.chi2FitTo(data1, RooFit::SumW2Error(kTRUE), RooFit::Range(x0, x1));
        P1.Print();

        RooGenericPdf fit_pdf2("fit_pdf2", fit2, RooArgSet(y2, P2));
        fit_pdf2.chi2FitTo(data2, RooFit::SumW2Error(kTRUE), RooFit::Range(x0, x1));
        P2.Print();
    }

    c->cd();        // go back to working area of whole canvas
    c1npad2->SetFillColor(kWhite);
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

    TF1 *line_1 = new TF1("line_1", "1", 0, 1);
    line_1->SetLineColor(kBlack);
    line_1->Draw("SAME");

    if (fit != "") {
        TString fit1 = fit;
        fit1 = fit1.ReplaceAll("P", TString("(") + TString::Format("%f", P1.getVal()) + ")");

        TString fit2 = fit;
        fit2 = fit2.ReplaceAll("P", TString("(") + TString::Format("%f", P2.getVal()) + ")");

        TString fitratio = TString("(") + fit1 + ")/(" + fit2 + ")";

        TF1 *ffitratio = new TF1("fitration", fitratio, x0, x1);
        ffitratio->SetLineColor(kBlue);
        ffitratio->SetLineWidth(1);
        ffitratio->Draw("SAME");
        
        /*
        // Show fits on histograms, not just on ratio:

        c1npad1->cd();
        TString normf1 = fit1 + "*"; normf1 += h1->GetBinWidth(1);
        TString normf2 = fit2 + "*"; normf2 += h2->GetBinWidth(1);
        TF1 *ff1 = new TF1("ff1",normf1,x0,x1);
        ff1->SetLineColor(kBlue);
        TF1 *ff2 = new TF1("ff2",normf2,x0,x1);
        ff2->SetLineColor(kMagenta);
        ff1->Draw("SAME");
        ff2->Draw("SAME");
        */
    }

    c->cd(); // go back to main canvas

    if (fit != "") {
        c1npad3->SetFillColor(kWhite);
        c1npad3->Draw(); // draw smaller pad in main canvas
        c1npad3->cd();   // set working area to smaller pad

        TString p_orig = "P_{orig} = " + TString::Format("%.2e", P1.getVal()) + " #pm " + TString::Format("%.2e", P1.getError());
        TString p_unweighted = "P_{unweighted} = " + TString::Format("%.2e", P2.getVal()) + " #pm " + TString::Format("%.2e", P2.getError());

        TPaveText *pt = new TPaveText(.05,.1,.95,.8);
        pt->SetTextSize(0.18);
        pt->SetFillColor(kWhite);
        pt->AddText(p_orig);
        pt->AddText(p_unweighted);
        
        pt->SetBorderSize(0);
        pt->Draw();

        c->cd();
    }
}
