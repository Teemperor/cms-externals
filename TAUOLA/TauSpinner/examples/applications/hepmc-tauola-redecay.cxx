/**
 * @author T. Przedzinski
 * @date 19 Nov 2013
 *
 * @brief Analyses data sample replacing tau decays with new ones
 *        generated using Tauola++. Output is saved to new file
 *
 * !!!IMPORTANT!!!
 * See README for a list of assumptions made by this program.
 *
 */
// HepMC IO_GenEvent header
#include "HepMC/IO_GenEvent.h"

// TAUOLA header
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"

#include <iostream>
#include <fstream>
#include <cstdio>

#include "CorrectEvent.h"
using namespace std;
using namespace Tauolapp;

/* This routine is used to introdcce experimental cuts.
   It is used in main() to generate new decays if the previously generated one
   does not pass these cuts */
int passed_cuts(HepMC::GenEvent *event)
{
  // Find tau pair
  HepMC::GenParticle *t1=NULL,*t2=NULL;
  for(HepMC::GenEvent::particle_iterator i = event->particles_begin();i!=event->particles_end();i++)
  {
    if(abs((*i)->pdg_id())==15 && (*i)->end_vertex() && abs((*(*i)->end_vertex()->particles_out_const_begin())->pdg_id())==16)
    {     
      if(!t1)  t1 = *i;
      else if(!t2) t2 = *i;
      else cout<<"Err! p1 & p2 already filled!"<<endl;
    }
  }
  if(!t1 || !t2) return -1;

  // Apply cuts on leptons
  // Calculate sum of pion momentum
  double px=0,py=0,pz=0,e=0;
  for(HepMC::GenVertex::particles_out_const_iterator i=t1->end_vertex()->particles_out_const_begin();i!=t1->end_vertex()->particles_out_const_end();i++)
  {
    HepMC::FourVector m1 = (*i)->momentum();
    double pt1 = sqrt(m1.px()*m1.px() + m1.py()*m1.py());

    if(abs((*i)->pdg_id())==11 && pt1<12.0) return 0;
    else if(abs((*i)->pdg_id())==13 && pt1<8.0 ) return 0;
    else if(abs((*i)->pdg_id())==211 || abs((*i)->pdg_id())==111)
    {
      px += m1.px();
      py += m1.py();
      pz += m1.pz();
      e  += m1.e ();
    }
  }

  double e_l = sqrt(px*px+py*py);

  // Apply cuts on leptons
  // Calculate sum of pion momentum
  double px2=0,py2=0,pz2=0,e2=0;
  for(HepMC::GenVertex::particles_out_const_iterator i=t2->end_vertex()->particles_out_const_begin();i!=t2->end_vertex()->particles_out_const_end();i++)
  {
    HepMC::FourVector m1 = (*i)->momentum();
    double pt1 = sqrt(m1.px()*m1.px() + m1.py()*m1.py());

    if(abs((*i)->pdg_id())==11 && pt1<12.0) return 0;
    else if(abs((*i)->pdg_id())==13 && pt1<8.0 ) return 0;
    else if(abs((*i)->pdg_id())==211 || abs((*i)->pdg_id())==111)
    {
      px2 += m1.px();
      py2 += m1.py();
      pz2 += m1.pz();
      e2  += m1.e ();
    }
  }

  double e_l2 = sqrt(px2*px2+py2*py2);

  // Apply p_t cuts on pi- pi0 pair: first tau p_t > 20, second tau p_t > 15
  if(e_l >e_l2 && (e_l <20.0 || e_l2<15.0) ) return 0;
  if(e_l2>e_l  && (e_l2<20.0 || e_l <15.0) ) return 0;

  return 1;
}

//----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  if(argc<3)
  {
    cout<<"Usage: "<<argv[0]<<" <input_file> <output_file> [<max_events>]"<<endl;
    exit(0);
  }
  
  int events_limit = 0;

  if(argc>3) events_limit = atoi(argv[3]);

  // Initialize Tauola
  Tauola::initialize();
  
  // Uncomment this line if producing sample
  // without polarization or spin correlations
  //Tauola::spin_correlation.setAll(false);

  // Modify branching ratios so that only channels 1-4 are generated 
  for(int i=5;i<30;i++) Tauola::setTauBr(i,0.0);

  HepMC::IO_GenEvent  input (argv[1],std::ios::in);
  HepMC::IO_GenEvent  output(argv[2],std::ios::out);
  
  int events_count     = 0;
  int events_generated = 0;
  int events_invalid   = 0;

  // Begin event loop
  while(true)
  {
    if(events_limit && events_count>=events_limit) break;

    // Read event from datafile
    HepMC::GenEvent *event = new HepMC::GenEvent();
    input.fill_next_event(event);
    if(input.rdstate()) break;

    // Convert units to GEV for Tauola++
    event->use_units(HepMC::Units::GEV,HepMC::Units::CM);
    
    // Uncomment this line to correct tau momentum and energy of tau decay products 
    // from effects of rounding errors. This function require checks if rounding errors
    // are indeed the problem and not some other effects. Responsability is to the user. 
    //CorrectEvent(event);

    // Generate events as many times as needed until one passess cuts
    while(true)
    {
      if(events_generated%1000==0) cout<<"EVT generated: "<<events_generated<<" accepted: "<<events_count<<" limit: "<<events_limit<<endl;

      TauolaHepMCEvent *tevent = new TauolaHepMCEvent(event);

      tevent->undecayTaus();
      tevent->decayTaus();
      
      delete tevent;

      events_generated++;

      // Apply cuts
      // Insert 'break;' before this line to ignore cuts
      int res = passed_cuts(event);

      // Faulty event - exit loop
      if(res==-1)
      {
        events_invalid++;
        break;
      }

      // Passed cuts - exit loop
      if(res>0) break;
    }
   
    output.write_event(event);
    events_count++;

    delete event;
  }
  
  Tauola::summary();
  cout<<"Invalid events in input file:           "<<events_invalid<<endl;
  cout<<"Events generated:                       "<<events_generated<<endl;
  cout<<"Passed cuts (and saved to output file): "<<events_count<<endl;
}

