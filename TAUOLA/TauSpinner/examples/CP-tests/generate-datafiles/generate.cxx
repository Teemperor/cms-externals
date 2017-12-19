#include "Tauola/Log.h"
#include "Tauola/Tauola.h"
#include "Tauola/TauolaHepMCEvent.h"

#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"

#include "HepMC/IO_GenEvent.h"
using namespace std;
using namespace Pythia8; 
using namespace Tauolapp;

void recursive_copy_daughters(HepMC::GenEvent *evt, HepMC::GenVertex *v)
{
  if(!v) return;
  
  evt->add_vertex(v);
  
  for(HepMC::GenVertex::particles_out_const_iterator i = v->particles_out_const_begin(); i!=v->particles_out_const_end(); i++)
  {
    recursive_copy_daughters(evt,(*i)->end_vertex());
  }
}

int main(int argc,char **argv){

  if(argc<5)
  {
    cout<<"Usage: "<<argv[0]<<" <output_file> <pythia_conf> <tau_decay_mode> <events>"<<endl;
    cout<<"       (note that this program is prepared for generation of H or Z/gamma decays only)"<<endl;
    exit(-1);
  }
  
  char *filename  = argv[1];
  int   events    = atoi(argv[4]);
  int   decayMode = atoi(argv[3]);

  Log::SummaryAtExit();

  // Initialization of pythia
  HepMC::Pythia8ToHepMC ToHepMC;
  Pythia pythia;
  Event& event = pythia.event;

  pythia.particleData.readString("15:mayDecay = off");
  pythia.readFile(argv[2]);

  pythia.init( 2212, 2212, 14000.0); // proton proton collisions

  Tauola::setSameParticleDecayMode    (decayMode);
  Tauola::setOppositeParticleDecayMode(decayMode);
  Tauola::setRadiation(false);
  Tauola::initialize();
  Tauola::spin_correlation.setAll(false);

  HepMC::IO_GenEvent file(filename,ios::out);

  // Begin event loop
  for (int iEvent = 0; iEvent < events; ++iEvent) {

    if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<endl;
    if(!pythia.next()) continue;

    // Convert event record to HepMC
    HepMC::GenEvent * HepMCEvt      = new HepMC::GenEvent();
    HepMC::GenEvent * HepMCEvt_copy = new HepMC::GenEvent();

    // Conversion needed if HepMC use different momentum units
    // than Pythia. However, requires HepMC 2.04 or higher.
    HepMCEvt->use_units(HepMC::Units::GEV,HepMC::Units::MM);

    ToHepMC.fill_next_event(event, HepMCEvt);

    // Run TAUOLA on the event
    TauolaHepMCEvent * t_event = new TauolaHepMCEvent(HepMCEvt);
    t_event->decayTaus();
    delete t_event; 

    // Cut smaller sub-event
    HepMC::GenParticle *beam1 = NULL, *beam2 = NULL;
    for(HepMC::GenEvent::particle_const_iterator p=HepMCEvt->particles_begin();p!=HepMCEvt->particles_end();p++)
    {
      if(!(*p)->end_vertex()) continue;
      int abs_pdg = abs( (*p)->pdg_id() );
      int pdg2 = (*(*p)->end_vertex()->particles_out_const_begin())->pdg_id();
 
      if( ((abs_pdg>0 && abs_pdg<7) || abs_pdg==21) && (pdg2==22 || pdg2==23 || pdg2==25 || pdg2==36) ) {
        if(!beam1) beam1 = (*p);
        else if(!beam2) beam2 = (*p);
        else {
          cout<<"Error: found more than two mothers of H/Z"<<endl;
          exit(-1);
        }
      }
    }
    
    if(!beam1 || !beam2) exit(-1);
    
    HepMC::GenVertex *v1 = new HepMC::GenVertex();
    v1->add_particle_out(beam1);

    HepMC::GenVertex *v2 = new HepMC::GenVertex();
    v2->add_particle_out(beam2);

    HepMCEvt_copy->add_vertex(v1);
    HepMCEvt_copy->add_vertex(v2);

    HepMCEvt_copy->set_beam_particles(beam1,beam2);
    recursive_copy_daughters(HepMCEvt_copy,beam1->end_vertex());
    recursive_copy_daughters(HepMCEvt_copy,beam2->end_vertex());
    
    if(iEvent==0) HepMCEvt_copy->print();

    file.write_event(HepMCEvt_copy);
    
    // Clean up
    delete HepMCEvt;
    delete HepMCEvt_copy;
  }

  pythia.statistics();
  Tauola::summary();
}

