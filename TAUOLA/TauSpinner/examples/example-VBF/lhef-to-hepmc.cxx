#include <iostream>
#include <vector>

// HepMC IO_GenEvent header
#include "HepMC/IO_GenEvent.h"

// LHE reader (Copyright (C) 2009-2013 Leif Lonnblad)
#include "LHEF.h"

using namespace std;

/** @brief Convert an event in LHEF format to HepMC */
HepMC::GenEvent* lhef2hepmc(LHEF::Reader &reader, int event_no)
{
    HepMC::GenEvent *evt = new HepMC::GenEvent();

    evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

    const double weight = reader.hepeup.weight();

    evt->weights().push_back(weight);
    evt->set_event_number(event_no);
    evt->set_event_scale (reader.hepeup.SCALUP);
    evt->set_alphaQCD    (reader.hepeup.AQCDUP);
    evt->set_alphaQED    (reader.hepeup.AQEDUP);

    // Create particles and vertices
    vector<HepMC::GenParticle*> particles;

    for (int i = 0; i < reader.hepeup.NUP; ++i)
    {
        HepMC::FourVector   momentum(reader.hepeup.PUP[i][0], reader.hepeup.PUP[i][1],reader.hepeup.PUP[i][2], reader.hepeup.PUP[i][3]);
        HepMC::GenParticle *p = new HepMC::GenParticle(momentum,  reader.hepeup.IDUP[i], reader.hepeup.ISTUP[i]);
        p->set_generated_mass(reader.hepeup.PUP[i][4]);

        particles.push_back(p);

        int mom1 = reader.hepeup.MOTHUP[i].first;
        int mom2 = reader.hepeup.MOTHUP[i].second;

        if( mom1 > 0 || mom2 > 0 ) {
            HepMC::GenVertex *vertex = particles[mom1-1]->end_vertex();

            // Create new vertex if not already present
            if( !vertex ) {
                vertex = new HepMC::GenVertex();
                evt->add_vertex(vertex);

                // Add mothers to new vertex
                if( mom1 > 0 ) vertex->add_particle_in(particles[mom1-1]);
                if( mom2 > 0 ) vertex->add_particle_in(particles[mom2-1]);
            }

            // Add daughter to new or previously existing vertex
            vertex->add_particle_out(p);
        }
    }

    if( particles.size()>2 ) evt->set_beam_particles(particles[0],particles[1]);

    // Convert to MEV
    evt->use_units(HepMC::Units::MEV, HepMC::Units::MM);

    //evt->print();

    return evt;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main(int argc, char **argv) {

    if(argc<3) {
        cout<<"Usage:    "<<argv[0]<<" <input_lhe_file> <output_hepmc_file>" << endl;
        cout<<"Consider: "<<argv[0]<<" events-VBF.lhe converted.hepmc" << endl;
        exit(-1);
    }

    char *input_filename  = argv[1];
    char *output_filename = argv[2];

    //-------------------------------------------------------------------------
    //- Initialization --------------------------------------------------------
    //-------------------------------------------------------------------------

    int event_no = 0;

    // Open I/O files
    LHEF::Reader       reader(input_filename);
    HepMC::IO_GenEvent writer(output_filename,ios::out);

    //-------------------------------------------------------------------------
    //- Event loop ------------------------------------------------------------
    //-------------------------------------------------------------------------
    while( reader.readEvent() ) {

        HepMC::GenEvent *evt = lhef2hepmc(reader,event_no);

        if( !evt ) break;

        writer.write_event(evt);

        ++event_no;

        if( event_no%1000 == 0 ) cout << "EVT: " << event_no << endl;
    }
  
    cout << "Done." << endl;
}
