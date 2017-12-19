#include "common.h"

// ROOT headers
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

#include <fstream>
#include <iostream>
#include <vector>

#include <cstdlib>

using namespace std;

struct FitRange {
    Double_t x0;
    Double_t x1;
};

TString get_next_line(ifstream &in) {
    char buf[256];
    in.getline(buf, 256);
    while (((buf[0] == '#') || (TString(buf) == ""))&& !(in.eof())) {
        in.getline(buf, 256);
    }
    TString result(buf);
    return result;
}

void process_rootfile(TString rootfile_name, TString plotdir, const vector<FitRange> &fit_ranges) {
    TFile rootfile(rootfile_name);

    const int lepton = 0;
    const int pion = 1;
    const int rho = 2;

    TString fits[] = {"(2 - 6 * x**2 + 4 * x**3 + 1/3 * (-1 + 9 * x**2 - 8 * x**3)"
            " - P * ( - 2/3 + 4 * x - 6 * x**2 + 8/3 * x**3 + 1/3 * (1 - 12 * x + 27 * x**2 - 16 * x**3)))",
            "(1 + P * (2 * x - 1))",
            ""};

    int Z_m[] = {lepton, lepton,
                      pion, pion, 
                      lepton, rho,
                      pion, rho, rho};
    
    int Z_p[] = {lepton, pion,
                      lepton, pion,
                      rho, lepton,
                      rho, pion, rho};

    int H_m[] = {lepton, lepton,
                      pion, pion,
                      lepton, rho,
                      pion, rho, rho};
    
    int H_p[] = {lepton, pion,
                      lepton, pion,
                      rho, lepton,
                      rho, pion, rho};

    int W[] = {lepton, pion, rho};

    TCanvas c("c","c",800,750);
    c.SetFillColor(kWhite);

    //W plots
    TH1D *W_mass = (TH1D *)rootfile.Get("W");
    if(W_mass && W_mass->GetEntries()>0) {
        c.SetLogy(1);
        W_mass->SetStats(kFALSE);
        W_mass->Draw();
        c.SaveAs(plotdir + "/W_mass.eps");
        c.SetLogy(0);

        TH1D *h1[3][3][2];
        for(int i=0;i<3;i++) {
            TString prefix = "h1_" + TString::Format("%d", i);
            for(int k=0;k<2;k++) {
                for (int j=0;j<3;j++) {
                    h1[i][j][k] = (TH1D *)rootfile.Get(prefix +
                                                 TString::Format("%d", j) + TString::Format("%d", k) + ";1");
                }

                draw3(&c, h1[i][0][k], h1[i][1][k], h1[i][2][k], 1, fits[W[i]], fit_ranges[W[i]].x0, fit_ranges[W[i]].x1);
                        
                TString plotname = h1[i][2][k]->GetName();
                c.SaveAs(plotdir + "/" + plotname + ".eps");
            }
        }

        // Print summary and save it to output file
        ofstream out( (plotdir + "/" + "input-W-event-count.txt").Data() );
    
        char buf[128];
        cout<<"            W Events:"<<endl;
        out <<"            W Events:"<<endl;
    
        sprintf(buf,"%10i  Total\n",(int)W_mass->GetEntries());
        cout<<buf;
        out <<buf;
    
        for(int k=0;k<2;k++)
            for(int i=0;i<3;i++) {
                sprintf(buf,"%10i  %s\n",(int)h1[i][0][k]->GetEntries(),h1[i][0][k]->GetTitle());
                cout<<buf;
                out <<buf;
            }
    
        out.close();
    }

    //Z plots
    TH1D *Z_mass = (TH1D *)rootfile.Get("Z");
    if(Z_mass && Z_mass->GetEntries()>0) {
        c.SetLogy(1);
        Z_mass->SetStats(kFALSE);
        Z_mass->Draw();
        c.SaveAs(plotdir + "/Z_mass.eps");
        c.SetLogy(0);


        TH2D *h2[9];
        TH1D *px_0[9], *px_1[9], *px_2[9];
        TH1D *py_0[9], *py_1[9], *py_2[9];
        for(int i = 0; i < 9; i++) {
            TString prefix = "h2_" + TString::Format("%d", i) + "0";
            h2[i] = (TH2D *)rootfile.Get(prefix + "0;1");
            px_0[i] = (TH1D *)rootfile.Get(prefix + "0_px;1");
            px_1[i] = (TH1D *)rootfile.Get(prefix + "1_px;1");
            px_2[i] = (TH1D *)rootfile.Get(prefix + "0_px_ratio;1");
            py_0[i] = (TH1D *)rootfile.Get(prefix + "0_py;1");
            py_1[i] = (TH1D *)rootfile.Get(prefix + "1_py;1");
            py_2[i] = (TH1D *)rootfile.Get(prefix + "0_py_ratio;1");

            TString plotname_prefix = "h2_" + TString::Format("%d", i) + "00_projection";
            draw3(&c, px_0[i], px_1[i], px_2[i], 0, fits[Z_m[i]], fit_ranges[Z_m[i]].x0, fit_ranges[Z_m[i]].x1);
            c.SaveAs(plotdir + "/" + plotname_prefix + "X.eps");

            draw3(&c, py_0[i], py_1[i], py_2[i], 0, fits[Z_p[i]], fit_ranges[Z_p[i]].x0, fit_ranges[Z_p[i]].x1);
            c.SaveAs(plotdir + "/" + plotname_prefix + "Y.eps");
        }

        // Print summary and save it to output file
        ofstream out((plotdir + "/" + "input-Z-event-count.txt").Data());
    
        char buf[128];
        cout<<"            Z Events:"<<endl;
        out <<"            Z Events:"<<endl;
    
        sprintf(buf,"%10i  Total\n",(int)Z_mass->GetEntries());
        cout<<buf;
        out <<buf;
    
        for(int i=0; i<9; i++) {
            sprintf(buf,"%10i  %s\n",(int)h2[i]->GetEntries(),h2[i]->GetTitle());
            cout<<buf;
            out <<buf;
        }
        out.close();
    }

    //H plots
    TH1D *H_mass = (TH1D *)rootfile.Get("H");
    if(H_mass && H_mass->GetEntries()>0) {
        c.SetLogy(1);
        H_mass->SetStats(kFALSE);
        H_mass->SetTitle("log10(tau tau mass) in Phi (Z) decays");
        H_mass->Draw();
        c.SaveAs(plotdir + "/Phi_mass.eps");
        c.SetLogy(0);

        TH1D *px_0[9], *px_1[9], *px_2[9];
        TH1D *py_0[9], *py_1[9], *py_2[9];
        for(int i = 0; i < 9; i++) {
            TString prefix = "h3_" + TString::Format("%d", i) + "0";
            px_0[i] = (TH1D *)rootfile.Get(prefix + "0_px;1");
            px_1[i] = (TH1D *)rootfile.Get(prefix + "1_px;1");
            px_2[i] = (TH1D *)rootfile.Get(prefix + "0_px_ratio;1");
            py_0[i] = (TH1D *)rootfile.Get(prefix + "0_py;1");
            py_1[i] = (TH1D *)rootfile.Get(prefix + "1_py;1");
            py_2[i] = (TH1D *)rootfile.Get(prefix + "0_py_ratio;1");

            TString plotname_prefix = "h3_" + TString::Format("%d", i) + "00_projection";
            draw3(&c, px_0[i], px_1[i], px_2[i], 0, fits[H_m[i]], fit_ranges[H_m[i]].x0, fit_ranges[H_m[i]].x1);
            c.SaveAs(plotdir + "/" + plotname_prefix + "X.eps");

            draw3(&c, py_0[i], py_1[i], py_2[i], 0, fits[H_p[i]], fit_ranges[H_p[i]].x0, fit_ranges[H_p[i]].x1);
            c.SaveAs(plotdir + "/" + plotname_prefix + "Y.eps");
        }
    }

}


int main(int argc, char **argv) {
    if(argc<2) {
        cout << "Usage:   " << argv[0] << " <config_file>" << endl;
        exit(0);
    }
  
    //---------------------------------------------------------------------------
    //- Configuration file parsing ----------------------------------------------
    //---------------------------------------------------------------------------
    ifstream in(argv[1]);
    if(!in.is_open()) {
        cout << "Config file " << argv[1] << " does not exist." << endl;
        exit(-1);
    }

    int conf_op=0;

    cout << endl;
    cout << " Tauspinner fit configuration:" << endl;
    cout << " ------------------------------------" << endl;

    TString plotdir;
    vector<TString> rootfile_names;
    vector<FitRange> fit_ranges;


    while (true) {
        TString buf = get_next_line(in);
        if (in.eof()) {
            break;
        }

        FitRange fit_range;

        switch (conf_op) {
        case 0:
            plotdir = buf;
            cout << " Output directory: " << plotdir << endl;
            break;
        case 1:
            fit_range.x0 = buf.Atof();
            buf = get_next_line(in);
            fit_range.x1 = buf.Atof();
            cout << " Next fitting range limits: " << TString::Format("%.2e", fit_range.x0) << " "
                 << TString::Format("%.2e", fit_range.x1) << endl;
            fit_ranges.push_back(fit_range);
            for (int i = 0; i < 2; i++) {
                buf = get_next_line(in);
                fit_range.x0 = buf.Atof();
                buf = get_next_line(in);
                fit_range.x1 = buf.Atof();
                cout << " Next fitting range limits: " << TString::Format("%.2e", fit_range.x0) << " "
                     << TString::Format("%.2e", fit_range.x1) << endl;
                fit_ranges.push_back(fit_range);
            }
            break;
        default:
            rootfile_names.push_back(buf);
            cout << " Input file #" << conf_op - 2 << ": " << buf << endl;
        }

        conf_op++;
    }
    cout << " ------------------------------------" << endl;
    cout << endl;
    in.close();

    // Copy config file to output directory
    ifstream conf_in (argv[1]);
    ofstream conf_out((plotdir + "/" + "fit-config-file.txt").Data());
    conf_out << conf_in.rdbuf();
    conf_out.close();
    conf_in.close();

    for (vector<TString>::iterator i = rootfile_names.begin(); i != rootfile_names.end(); ++i) {
        cout << endl;
        cout << "===========================================" << endl;
        cout << "Processing " << *i << endl;
        cout << "===========================================" << endl;
        cout << endl;
        process_rootfile(*i, plotdir, fit_ranges);
    }
}
