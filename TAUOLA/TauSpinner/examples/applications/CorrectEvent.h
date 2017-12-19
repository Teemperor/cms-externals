/**
 * @author T. Przedzinski
 * @date 19 Nov 2013
 *
 * @brief Corrects energy of the particles based on 
 *        generated_mass() information from HepMC::GenEvent
 *
 * Introduces class FourVector for manipulating momentum
 * of the particles stored in HepMC file.
 *
 * Included by hepmc-tauola-redecay.cxx but the function call
 * is commented out by default
 * This algorithm may create severe problems if kinematical impurities of
 * the events are not the consequence of rounding errors only.
 *
 * Algorithm may be adaptable to use in `de-brem' cases. User should take care 
 * that hardest (of largest p_T with respect to lepton) photon is left for
 * use in TauSpinner for leptonic channels only, 
 * its energy must be above Tauola++ threshold.
 * 
 * Also final state pi- pi0 gamma nu should not be de-bremmed as it 
 * encapsulates omega -> pi0 gamma decay.
 *
 * Algorithm is not balanced for bremsstrahlung. It corrects kinematic 
 * re-scaling momentum of tau-neutrino. Even though it is consistent with 
 * LL of QED it is far form best choice. For non-collinear photons
 * algorithm will have a tendency to shift charged momentum to smaller
 * and tau neutrino to higher than optimal values. For collinear
 * emissions it should perform better.
 *
 * There is plenty of room for improvements, including rather straightforward
 * and well understood steps.
 */
#include "HepMC/GenParticle.h"
#include "HepMC/GenEvent.h"
#include <vector>
#include <iostream>
#include "Tauola/f_Variables.h"
using namespace std;

/** Threshold above which the correction is considered too large */
const double CORRECTION_THRESHOLD = 0.0001;

/** Class for storing information on particle momentum, including mass */
class FourVector
{
public:
  FourVector();
  FourVector(double _x, double _y, double _z, double _e, double _m);
  FourVector(const HepMC::FourVector& x);
  FourVector(const FourVector& x);
  FourVector operator+(const FourVector& y);
  FourVector operator-(const FourVector& y);
  FourVector operator+=(const FourVector& x);
  FourVector operator-=(const FourVector& x);
  double recalculated_mass();
  void print();
public:
  double px,py,pz,e,m;
};

FourVector::FourVector():px(0.),py(0.),pz(0.),e(0.),m(0.)                                                         {}
FourVector::FourVector(double _x, double _y, double _z, double _e, double _m=0.):px(_x),py(_y),pz(_z),e(_e),m(_m) {}
FourVector::FourVector(const HepMC::FourVector& x):px(x.px()),py(x.py()),pz(x.pz()),e(x.e()),m(0.)                {}
FourVector::FourVector(const FourVector& x):px(x.px),py(x.py),pz(x.pz),e(x.e),m(x.m)                              {}

FourVector FourVector::operator+(const FourVector& y)
{
  return FourVector(px+y.px,py+y.py,pz+y.pz,e+y.e,0.);
}

FourVector FourVector::operator-(const FourVector& y)
{
  return FourVector(px-y.px,py-y.py,pz-y.pz,e-y.e,0.);
}

FourVector FourVector::operator+=(const FourVector& x)
{
  return (*this)=(*this)+x;
}

FourVector FourVector::operator-=(const FourVector& x)
{
  return (*this)=(*this)-x;
}

double FourVector::recalculated_mass()
{
  double mm = e*e-px*px-py*py-pz*pz;
  if(mm>0.0)      return  sqrt(mm);
  else if(mm<0.0) return -sqrt(-mm);
  else return 0.0;
}

void FourVector::print()
{
  printf("%+13.6e %+13.6e %+13.6e %+13.6e m: %+13.6e recalculated_m: %+13.6e\n",
           px,    py,     pz,     e,         m,                      recalculated_mass() );
}

/*******************************************************************************
  Get daughters of HepMC::GenParticle

  Recursively searches for final-state daughters of 'x'
*******************************************************************************/
inline vector<HepMC::GenParticle*> getDaughters(HepMC::GenParticle *x)
{
  vector<HepMC::GenParticle*> daughters;
  if(!x->end_vertex()) return daughters;

  // Check decay products of 'x'
  for(HepMC::GenVertex::particles_out_const_iterator p = x->end_vertex()->particles_out_const_begin(); p!=x->end_vertex()->particles_out_const_end(); ++p)
  {
    HepMC::GenParticle *pp = *p;
    HepMC::FourVector   mm = pp->momentum();

    // If the daughter of 'x' has its end vertex - recursively read
    // all of its daughters. Ignore decays of pi0 if present
    if( pp->end_vertex() && pp->pdg_id()!=111)
    {
      vector<HepMC::GenParticle*> sub_daughters = getDaughters(pp);
      daughters.insert(daughters.end(),sub_daughters.begin(),sub_daughters.end());
    }
    // Otherwise - add this particle to the list of daughters.
    else daughters.push_back(pp);
  }

  return daughters;
}

/*****************************************************************************
  Correct event
  Performs following operation:
  - for each stable tau decay product its energy is recalculated from the mass
    and momentum
  - if method = 1 (default), the four-momentum of tau is corrected using
    tau mass. The difference is compensated by modifying tau neutrino
    four-momentum
    if method != 1, the four-momentum of tau is recalculated from the sum of
    four-momenta of its decay products
  - checks if resulting operation did not have introduced sizeable
    modifications. If it does, warning message is printed. This may indicate
    problems in the production file other than rounding error
*****************************************************************************/
void CorrectEvent(HepMC::GenEvent *evt, int method = 1)
{
  HepMC::GenParticle *pX    = NULL;
  HepMC::GenParticle *pTaup = NULL;
  HepMC::GenParticle *pTaum = NULL;
  HepMC::GenParticle *pTaup_nu = NULL;
  HepMC::GenParticle *pTaum_nu = NULL;
  
  // Get tau pair and intermediate boson
	for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();
	      p != evt->particles_end(); ++p )
	{
    if( (*p)->pdg_id() == 15 && (*p)->end_vertex() && abs((*(*p)->end_vertex()->particles_out_const_begin())->pdg_id())==16) pTaup = *p;
    if( (*p)->pdg_id() ==-15 && (*p)->end_vertex() && abs((*(*p)->end_vertex()->particles_out_const_begin())->pdg_id())==16) pTaum = *p;
	}

  if(!pTaup || !pTaum) return;
  
  // Get tau daughters
  vector<HepMC::GenParticle*> taumDaughters = getDaughters(pTaum);
  vector<HepMC::GenParticle*> taupDaughters = getDaughters(pTaup);

  FourVector taum_sum,taup_sum;
  FourVector taum_nu,taup_nu;

  // Correct tau+ daughters
  for(unsigned int i=0;i<taumDaughters.size();i++)
  {
    HepMC::FourVector m = taumDaughters[i]->momentum();
    FourVector buf = m;
    buf.m = taumDaughters[i]->generated_mass();
    
    // Correct energy
    buf.e = sqrt(buf.m*buf.m + buf.px*buf.px + buf.py*buf.py + buf.pz*buf.pz);
    
    m.setE(buf.e);
    
    taumDaughters[i]->set_momentum(m);
    
    taum_sum += buf;
    
    // Store tau neutrino 4-vector
    if( abs(taumDaughters[i]->pdg_id())==16)
    {
      pTaum_nu = taumDaughters[i];
      taum_nu  = buf;
    }
  }
  
  // Correct tau- daughters
  for(unsigned int i=0;i<taupDaughters.size();i++)
  {
    HepMC::FourVector m = taupDaughters[i]->momentum();
    FourVector buf = m;
    buf.m = taupDaughters[i]->generated_mass();

    // Correct energy
    buf.e = sqrt(buf.m*buf.m + buf.px*buf.px + buf.py*buf.py + buf.pz*buf.pz);
    
    m.setE(buf.e);
    
    taupDaughters[i]->set_momentum(m);
    
    taup_sum += buf;
    
    // Store tau neutrino 4-vector
    if( abs(taupDaughters[i]->pdg_id())==16)
    {
      pTaup_nu = taupDaughters[i];
      taup_nu  = buf;
    }
  }
  
  // Check the size of the corrections
  FourVector taum_diff = pTaum->momentum();  
  FourVector taup_diff = pTaup->momentum();
  
  taum_diff -= taum_sum;
  taup_diff -= taup_sum;
  
  double diff1 = sqrt(taum_diff.px*taum_diff.px + taum_diff.py*taum_diff.py + taum_diff.pz*taum_diff.pz + taum_diff.e*taum_diff.e);
  double diff2 = sqrt(taup_diff.px*taup_diff.px + taup_diff.py*taup_diff.py + taup_diff.pz*taup_diff.pz + taup_diff.e*taup_diff.e);
  
  //cout<<"DIFF: "<<diff1<<" "<<diff2<<endl;
  
  if( diff1>CORRECTION_THRESHOLD*taum_sum.e || diff2>CORRECTION_THRESHOLD*taup_sum.e )
  {
    cout<<"CorrectEvent: WARNING! Large correction: "<<diff1/taum_sum.e<<" "<<diff2/taup_sum.e<<endl;
  }
  
  if(method==1)
  {
    //cout.precision(12); 
    //cout<<"BEFORE TAU 4-VECTOR CORRECTION:"<<endl;
    //cout<<"NU TAU+:        "<<taum_nu.px<<" "<<taum_nu.py<<" "<<taum_nu.pz<<" "<<taum_nu.e<<endl;
    //cout<<"NU TAU-:        "<<taup_nu.px<<" "<<taup_nu.py<<" "<<taup_nu.pz<<" "<<taup_nu.e<<endl;
    //cout<<"TAU+:        "<<taum_sum.px<<" "<<taum_sum.py<<" "<<taum_sum.pz<<" "<<taum_sum.e<<endl;
    //cout<<"TAU-:        "<<taup_sum.px<<" "<<taup_sum.py<<" "<<taup_sum.pz<<" "<<taup_sum.e<<endl;

    // Correct tau momentum
    bool furtherp=false;
    bool furtherm=false;
    double tau_mass    = Tauolapp::parmas_.amtau;
    double taum_mass   = taum_sum.recalculated_mass();
    double taum_buf    = taum_sum.e*taum_nu.e - taum_sum.px*taum_nu.px - taum_sum.py*taum_nu.py - taum_sum.pz*taum_nu.pz;
    double taum_lambda = -(taum_mass*taum_mass - tau_mass*tau_mass)/2/taum_buf;
    if ( taum_lambda < -1.0) {
      cout<<"CorrectEvent.h failed, lambda- =" <<taum_lambda<<endl;
      // typically that means, the neutrino momentum was too small
      // to correct energy-momentum non-conservation. It may happen
      // neutrino momentum can be arbitrary small.
      // algorithm acting on other decay products, massive and decay 
      // channel dependent may be nedded and installed in this place.
      taum_lambda=-0.95;
      furtherm=true;  // alternatively corrections can be installed later.
    }
    taum_sum.px = taum_sum.px + taum_lambda*taum_nu.px;
    taum_sum.py = taum_sum.py + taum_lambda*taum_nu.py;
    taum_sum.pz = taum_sum.pz + taum_lambda*taum_nu.pz;
    taum_sum.e  = taum_sum.e  + taum_lambda*taum_nu.e;

    taum_nu.px  = taum_nu.px * (1+taum_lambda);
    taum_nu.py  = taum_nu.py * (1+taum_lambda);
    taum_nu.pz  = taum_nu.pz * (1+taum_lambda);
    taum_nu.e   = taum_nu.e  * (1+taum_lambda);

    double taup_mass   = taup_sum.recalculated_mass();
    double taup_buf    = taup_sum.e*taup_nu.e - taup_sum.px*taup_nu.px - taup_sum.py*taup_nu.py - taup_sum.pz*taup_nu.pz;
    double taup_lambda = -(taup_mass*taup_mass - tau_mass*tau_mass)/2/taup_buf;
    if ( taup_lambda < -1.0) {
      cout<<"CorrectEvent.h failed, lambda+ =" <<taup_lambda<<endl;
      // typically that means, the neutrino momentum was too small
      // to correct energy-momentum non-conservation. It may happen
      // neutrino momentum can be arbitrary small.
      // algorithm acting on other decay products, massive and decay 
      // channel dependent may be nedded and installed in this place.
      taup_lambda=-0.95;
      furtherp=true;  // alternatively corrections can be installed later.
    }


    taup_sum.px = taup_sum.px + taup_lambda*taup_nu.px;
    taup_sum.py = taup_sum.py + taup_lambda*taup_nu.py;
    taup_sum.pz = taup_sum.pz + taup_lambda*taup_nu.pz;
    taup_sum.e  = taup_sum.e  + taup_lambda*taup_nu.e;

    taup_nu.px  = taup_nu.px * (1+taup_lambda);
    taup_nu.py  = taup_nu.py * (1+taup_lambda);
    taup_nu.pz  = taup_nu.pz * (1+taup_lambda);
    taup_nu.e   = taup_nu.e  * (1+taup_lambda);


    //cout<<"AFTER TAU 4-VECTOR CORRECTION:"<<endl;
    //cout<<"LAMBDA TAU+: "<<taum_lambda<<endl;
    //cout<<"LAMBDA TAU-: "<<taup_lambda<<endl;
    //cout<<"NU TAU+:        "<<taum_nu.px<<" "<<taum_nu.py<<" "<<taum_nu.pz<<" "<<taum_nu.e<<endl;
    //cout<<"NU TAU-:        "<<taup_nu.px<<" "<<taup_nu.py<<" "<<taup_nu.pz<<" "<<taup_nu.e<<endl;
    //cout<<"TAU+:        "<<taum_sum.px<<" "<<taum_sum.py<<" "<<taum_sum.pz<<" "<<taum_sum.e<<endl;
    //cout<<"TAU-:        "<<taup_sum.px<<" "<<taup_sum.py<<" "<<taup_sum.pz<<" "<<taup_sum.e<<endl;
    
    // Save modified nu_tau 4-vectors
    HepMC::FourVector fv_taum_nu(taum_nu.px,taum_nu.py,taum_nu.pz,taum_nu.e);
    HepMC::FourVector fv_taup_nu(taup_nu.px,taup_nu.py,taup_nu.pz,taup_nu.e);
    
    pTaum_nu->set_momentum(fv_taum_nu);
    pTaup_nu->set_momentum(fv_taup_nu);
    if(furtherm){
      // correcting with tau-neutrino is insufficient for Taum,
      // first we implement whatever is available.
      // HepMC::FourVector fv_taum(taum_sum.px,taum_sum.py,taum_sum.pz,taum_sum.e);
      // pTaum->set_momentum(fv_taum);
    }
    if(furtherp){
      // correcting with tau-neutrino is insufficient for Taup,
      // first we implement whatever is available.
      // HepMC::FourVector fv_taup(taup_sum.px,taup_sum.py,taup_sum.pz,taup_sum.e);
      // pTaup->set_momentum(fv_taup);
    }
  }
  
  // Save modified tau 4-vectors
  HepMC::FourVector fv_taum(taum_sum.px,taum_sum.py,taum_sum.pz,taum_sum.e);
  HepMC::FourVector fv_taup(taup_sum.px,taup_sum.py,taup_sum.pz,taup_sum.e);
  
  pTaum->set_momentum(fv_taum);
  pTaup->set_momentum(fv_taup);
}
