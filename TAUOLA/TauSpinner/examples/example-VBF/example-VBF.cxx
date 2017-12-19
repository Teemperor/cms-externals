#include <iostream>
#include <vector>

// HepMC IO_GenEvent header
#include "HepMC/IO_GenEvent.h"

// TAUOLA header
#include "Tauola/Tauola.h"

// TauSpinner headers
#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/tau_reweight_lib.h"
#include "TauSpinner/vbfdistr.h"
// For SPIN2 code placed in directory example-VBF:
// #include "spin2distr.h"  // will work once `export SPIN2LIB=...' 

// LHAPDF header
#include "LHAPDF/LHAPDF.h"

#include "read_particles_for_VBF.h"
using namespace std;
using namespace Tauolapp;
using namespace TauSpinner;

/** Example function that can be used to modify/replace alpha_s calculation of vbfdist 
    WARNING: 4-momenta must be in GeV */
void  alphasModif(double Q2,int scalePDFOpt, int KEY){
  //any textbook LL calculation
    const double PI=3.14159265358979324;
  // number of flavors
    const double Nf = 5;
  // alphas_s at mZ
    const double AlphasMZ = 0.118;
    const double MZ =91.1876;
  // alphas_s at scale Q2 (GeV^2)
    double alfas=AlphasMZ / ( 1 + AlphasMZ/(4*PI) * (11 - 2./3 * Nf) * log(Q2/(MZ*MZ)));
    if(scalePDFOpt==0) alfas = 0.118;
  // test Alphas ( Q2 = 1000^2 ) = 0.0877445

  //reasonable choice can be:
  //  alfas = LHAPDF::alphasPDF(sqrt(Q2));

  //---------------------------------------
    if(params_r_.as != alfas){ // we pass alpha_s to calculation of amplitudes
       params_r_.as = alfas;
	//  reinitialize constants couplings for amplitude calculations
	vbf_reinit_(&KEY);
      }
  }




/** Example function that can be used to modify/replace matrix element calculation of vbfdist */
double vbfdistrModif(int I1, int I2, int I3, int I4, int H1, int H2, double P[6][4], int KEY, double vbfdistr_result)
{
  //       vbfdistr_result is the result of library calculation. It use initialization set with 
  //       vbfinit_(&ref,&variant);
  //       choice  of &variant
    return vbfdistr_result;
}

//-----------------------------------------------------------------
//replacement of default (not best quality) random number generator
// #include <TRandom.h>
// TRandom gen;
// double randomik(){
// return gen.Rndm();
// }
//-----------------------------------------------------------------
int main(int argc, char **argv) {

    if(argc<2) {
        cout<<"Usage:    "<<argv[0]<<" <input_file> [<events_limit>]" << endl;
        cout<<"Consider: "<<argv[0]<<" events-VBF.dat 10" << endl;
        exit(-1);
    }

    char *input_filename = argv[1];
    int   events_limit   = 0;
    if(argc>2) events_limit = atoi(argv[2]);

    //---------------------------------------------------------------------------
    //- Initialization ----------------------------------------------------------
    //---------------------------------------------------------------------------

    // Initialize Tauola
    Tauola::initialize();
    Tauola::spin_correlation.setAll(false);

    // Initialize random numbers:
    // ##1##
    // Important when you re-decay taus: set seed fortauola-fortran random number generator RANMAR
    // int ijklin=..., int ntotin=..., int ntot2n=...; /
    // Tauola::setSeed(ijklin,ntotin,ntot2n);
    // Tauola::setSeed(time(NULL), 0, 0);

    // ##2##
    // Important when you use attributed by TauSpinner  helicities
    // Replace C++ Tauola Random generator with your own (take care of seeds). Prepared method: 
    // gen.SetSeed(time(NULL));
    // Tauola::setRandomGenerator( randomik );

    // Initialize LHAPDF
    // string name="MSTW2008nnlo90cl.LHgrid";
    string name="cteq6ll.LHpdf";
    // choice used for events-VBF.lhe which is tiny, thus it is not 
    // string name="MSTW2008nlo68cl.LHgrid"; //      statistically important
    LHAPDF::initPDFSetByName(name);

    double CMSENE = 13000.0;  // 14000.0;
    bool   Ipp    = true;
    int    Ipol   = 1;
    int    nonSM2 = 0;
    int    nonSMN = 0;

    // Initialize TauSpinner
    initialize_spinner(Ipp, Ipol, nonSM2, nonSMN,  CMSENE);

    int ref=4;      // EW scheme to be used for default vbf calculation.
                    // WARNING: electroweak scheme adopted to typical choice
                    // made in madgraph 2016/2017 is 1. 
    int variant =4; // EW scheme to be used in optional matrix element reweighting (nonSM2=1). Then
                    //   for vbf calculation, declared above prototype method vbfdistrModif (or user function)
                    //   will be used. At its disposal result of calculation with variant of  EW scheme will be available.
    vbfinit_(&ref,&variant);

    int QCDdefault=1; // QCD scheme to be used for default vbf calculation.
    int QCDvariant=1; // QCD scheme to be used in optional matrix element reweighting (nonSM2=1).
    setPDFOpt(QCDdefault,QCDvariant);

    // Set function that modifies/replaces Matrix Element calculation of vbfdistr
    // TauSpinner::set_vbfdistrModif(vbfdistrModif);

    // Set function that modifies/replaces Matrix Element calculation of vbfdistr with code of SPIN2
    // spin2init_(&ref,&variant);
    // TauSpinner::set_vbfdistrModif(SPIN2::spin2distr);
 
 // Set function that modifies/replaces alpha_s calculation of vbfdistr
    // TauSpinner::set_alphasModif(alphasModif);

    // Open I/O files  (in our example events are taken from "events.dat")
    HepMC::IO_GenEvent input_file(input_filename,std::ios::in);

    if(input_file.rdstate()) {
        cout<<endl<<"ERROR: file "<<input_filename<<" not found."<<endl<<endl;
        exit(-1);
    }
    int events_read =0;

    int    events_count = 0;
    double wt_sum       = 0.0;

    //---------------------------------------------------------------------------
    //- Event loop --------------------------------------------------------------
    //---------------------------------------------------------------------------
    while( !input_file.rdstate() ) {
        double    WT      = 1.0;
        double    W[2][2] = { { 0.0 } };

        SimpleParticle p1, p2, X, p3, p4, tau1, tau2;
        vector<SimpleParticle> tau1_daughters, tau2_daughters;

        int status = read_particles_for_VBF(input_file,p1,p2,X,p3,p4,tau1,tau2,tau1_daughters,tau2_daughters);
	++events_read;

// WARNING: meaning of status may depend on the variant of  read_particles_for_VBF().
        if( status == 1 ) break;                      // variant A
	//        if     ( status == 0 ) break;       // variant B
	//        else if( status == 1 ) { continue;} // variant B

	// At present we do not need X-boson  4-momentum at all but it may be useful for care of QED FSR
	// to construct bare taus for production matrix element calculation.
	//  read_particles_for_VBF() may constructs X but it should not be taken as granted.
	// X Pdgid defines if Higgs production matrix element has to be used.
        if(X.pdgid()==NULL) {
	
          X.setE(tau1.e()+tau2.e());
          X.setPx(tau1.px()+tau2.px());
	  X.setPy(tau1.py()+tau2.py());
	  X.setPz(tau1.pz()+tau2.pz());
          X.setPdgid(23);
	  
	
	};


        // Fix mass of jet partons (must be 0)
        // We do it naively by changing energy of the quarks
        double p3_pp = p3.px()*p3.px() + p3.py()*p3.py() + p3.pz()*p3.pz();
        double p3_m2 = p3.e()*p3.e() - p3_pp;
        double p4_pp = p4.px()*p4.px() + p4.py()*p4.py() + p4.pz()*p4.pz();
        double p4_m2 = p4.e()*p4.e() - p4_pp;

        if(fabs(p3_m2) > 2.0e-5) p3.setE( sqrt(p3_pp) );
        if(fabs(p4_m2) > 2.0e-5) p4.setE( sqrt(p4_pp) );

	// setNonSMkey(1);  // e.g. for use of SPIN2 ME
        WT = calculateWeightFromParticlesVBF(p3, p4, X, tau1, tau2, tau1_daughters, tau2_daughters);

	// Advanced use of matrix element calculation.  
        //     KEY=1 for Higgs hypothesis, KEY=3 Higgs hypothesis with user ME, 
        //     KEY=0 for DY-like process , KEY=2 DY-process with  user ME
	// int KEY=0;
        // getME2VBF(p3, p4, X, tau1, tau2, W, KEY);


        wt_sum += WT;
        ++events_count;
        if( events_count%100 == 0 ) cout << "EVT: " << events_count << endl;
        if( events_limit && events_count >= events_limit ) break;
    }
    
    cout<<endl<<"No of events read from the file: "<<events_read<<endl;
    cout<<endl<<"No of events processed for spin weight: "<<events_count<<endl;
    cout<<      "WT average for these processed events: "<<wt_sum/events_count<<endl;
}
