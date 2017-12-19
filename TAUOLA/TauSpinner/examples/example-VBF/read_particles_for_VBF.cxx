#include <iostream>
#include <vector>
#include <list>

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"

#include "TauSpinner/SimpleParticle.h"

#include "read_particles_for_VBF.h"
using namespace std;
using namespace TauSpinner;
const int pdgid_higgs=25;

SimpleParticle SimpleParticleFromHepMC(HepMC::GenParticle *p) {
    return SimpleParticle( p->momentum().px(), p->momentum().py(), p->momentum().pz(), p->momentum().e(), p->pdg_id() );
}

HepMC::GenParticle* findLastSelf(HepMC::GenParticle *part) {
    if( !part ) return NULL;

    HepMC::GenParticle *ret = part;
    while( ret->end_vertex() ) {
        HepMC::GenParticle *p = *ret->end_vertex()->particles_out_const_begin();
        if( ret->pdg_id() == p->pdg_id() ) ret = p;
        else break;
    }

    return ret;
}

bool findBeams(HepMC::GenEvent *evt, HepMC::GenParticle *&beam1, HepMC::GenParticle *&beam2) {
    if( !evt ) return false;
    for(HepMC::GenEvent::particle_const_iterator it = evt->particles_begin();
                                                 it!= evt->particles_end();
                                               ++it) {
        if(!(*it)->end_vertex()) continue;

        int pdg  = abs( (*it)->pdg_id() );
        int pdg2 = (*(*it)->end_vertex()->particles_out_const_begin())->pdg_id();

        if( ((pdg< 7 && pdg>0) || pdg==21) && (pdg2==21||pdg2==22||pdg2==23||abs(pdg2)<15|| abs(pdg2)<7||pdg2==pdgid_higgs ) ) {
            if(!beam1) beam1 = (*it);
            else {
                beam2 = (*it);
		break;
            }
        }
    }
    if( !beam1 || !beam2 ) {
        cout << "ERROR: Beams not found" << endl;
        return false;
    }
    return true;
}

/******************************************************* 
  findBosonAndJets: search for the Higgs boson and jets 
********************************************************/
bool findBosonAndJets(HepMC::GenParticle *beam1, HepMC::GenParticle *&boson, HepMC::GenParticle *&jet1, HepMC::GenParticle *&jet2) {
    if( !beam1 ) return false;

    //we search for jets originating from beams directly
    for(HepMC::GenVertex::particles_out_const_iterator it  = beam1->end_vertex()->particles_out_const_begin();
                                                       it != beam1->end_vertex()->particles_out_const_end();
                                                     ++it) {
        if( (*it)->pdg_id() == pdgid_higgs ) boson = *it;
		else if(abs((*it)->pdg_id()) >21) {}
		else if(abs((*it)->pdg_id()) ==15) {}
        else {
	  if     (!jet1&& abs((*it)->pdg_id()) !=15 && abs((*it)->pdg_id()) <22 ) jet1 = (*it);
	  else if(!jet2&& abs((*it)->pdg_id()) !=15 && abs((*it)->pdg_id()) <22 ) jet2 = (*it);
            else {
                cout << " findBosonAndJets: ERROR section A: more than 2 jets" << endl;
                return false;
            }
        }
	// taus mean that there is end of work
	//	if( abs((*it)->pdg_id()) ==15) break; 
	// if any of beam daughters is intermediate we search its daughters for jets as well 
	// skip self decays
	if(*it ) {
	    HepMC::GenParticle *Interm=*it;

	    	    
	    while(true) {
	      if(Interm->end_vertex()){ 	    
		HepMC::GenParticle *p = *Interm->end_vertex()->particles_out_const_begin();
		// second alternative is prepared for Higgstrahlung or  triple boson coupling.
		if( p->pdg_id() == Interm->pdg_id() || abs(p->pdg_id())> 20 ) {
		  // But triple boson coupling==triple boson vertex is not present in our events 
		  Interm = p;
		}
		else break;
	      }
	      else break;
	    }
	    if (Interm->end_vertex()&& abs(Interm->pdg_id()) !=15 ){	    	  
	for(HepMC::GenVertex::particles_out_const_iterator itt  = Interm->end_vertex()->particles_out_const_begin();
                                                	   itt != Interm->end_vertex()->particles_out_const_end();
	                                                 ++itt) {
	  //  cout << "Internal loop:  it->" << abs((*it)->pdg_id())<<endl;
	    if(    (*itt)->pdg_id() == pdgid_higgs) ;
		   else if( abs((*itt)->pdg_id()) ==15 ){break;} 
                  else if(!jet1&& abs((*itt)->pdg_id()) !=15  ) jet1 = (*itt);
	  	  else if(!jet2&& abs((*itt)->pdg_id()) !=15  ) jet2 = (*itt);
	  	  else {
	                cout << "findBosonAndJets: ERROR section B: more than 2 jets?" << endl;
                        exit(0);
	              return false;
		  }

	}
	    	  
	}
	}
    }

    if(  !jet1 || !jet2 ) {
        cout << "ERROR:  jets not found" << endl;
        return false;
    }
    return true;
}

bool findTaus(HepMC::GenParticle *boson, HepMC::GenParticle *&tau1, HepMC::GenParticle *&tau2) {
    if( !boson ) return false;
    while(true) {
      if (boson->end_vertex()){
        HepMC::GenParticle *p = *boson->end_vertex()->particles_out_const_begin();
	// second alternative is prepared for Higgstrahlung or  triple boson coupling.
	if( (p->pdg_id() == boson->pdg_id()&&(p->end_vertex()))|| abs(p->pdg_id())> 21){
	  //  if( (p->pdg_id() == boson->pdg_id()&& abs(p->pdg_id())>21)||(abs(boson->pdg_id())>21 && abs(p->pdg_id())> 21)){
	  // But triple boson coupling==triple boson vertex is not present in our events.
	  boson = p;
	}

        else{	 break;}
      }
     }
    for(HepMC::GenVertex::particles_out_const_iterator it  = boson->end_vertex()->particles_out_const_begin();
                                                       it != boson->end_vertex()->particles_out_const_end();
                                                     ++it) {
        if( abs((*it)->pdg_id()) == 15 ) {
            if(!tau1) tau1 = *it;
            else {
                tau2 = *it;
                break;
            }
        }
    }

    if( !tau1 || !tau2 ) {
        cout << "ERROR: tau pair not found" << endl;
        return false;
    }

    tau1 = findLastSelf(tau1);
    tau2 = findLastSelf(tau2);

    return true;
}

/*******************************************************************************
 Get daughters of HepMC::GenParticle
 
 Recursively searches for final-state daughters of 'x'
 *******************************************************************************/
inline vector<SimpleParticle> *getDaughters(HepMC::GenParticle *x)
{
	vector<SimpleParticle> *daughters = new vector<SimpleParticle>();
	if(!x->end_vertex()) return daughters;
	
	// Check decay products of 'x'
	for(HepMC::GenVertex::particles_out_const_iterator p = x->end_vertex()->particles_out_const_begin(); p!=x->end_vertex()->particles_out_const_end(); ++p)
	{
		HepMC::GenParticle *pp = *p;
		HepMC::FourVector   mm = pp->momentum();
		
		// If the daughter of 'x' has its end vertex - recursively read
		// all of its daughters.
		if( pp->end_vertex() && pp->pdg_id()!=111)
		{
			vector<SimpleParticle> *sub_daughters = getDaughters(pp);
			daughters->insert(daughters->end(),sub_daughters->begin(),sub_daughters->end());
			
			delete sub_daughters;
		}
		// Otherwise - add this particle to the list of daughters.
		else
		{
			SimpleParticle tp( mm.px(), mm.py(), mm.pz(), mm.e(), pp->pdg_id() );
			daughters->push_back(tp);
		}
	}
	
	return daughters;
}


/**********************************************************************************
 getDaughters_nonrecursive is a non recursive versione of getDaughters.
 It does not explore the complete tree in multi decays.
 **********************************************************************************/
 void getDaughters_nonrecursive(HepMC::GenParticle *p, vector<SimpleParticle> &daughters) {
	if(!p || !p->end_vertex()) return;
	
	if (!p) {
		cout << "in read_particles_for_VBF.cxx: getDaughters finds !p == TRUE" << endl;
	}
	if (!p->end_vertex()) {
		cout << "in read_particles_for_VBF.cxx: getDaughters finds !p->end_vertex() == TRUE" << endl;
	}

	while(true) {
		HepMC::GenParticle *pp = *p->end_vertex()->particles_out_const_begin();
		if( pp->pdg_id() == p->pdg_id() ) {
			p = pp;
			
			cout << "pp->pdg_id() == p->pdg_id()" << endl;
			
		}
		else break;
	}

	for(HepMC::GenVertex::particles_out_const_iterator it  = p->end_vertex()->particles_out_const_begin();
													   it != p->end_vertex()->particles_out_const_end();
													 ++it) {
		daughters.push_back( SimpleParticleFromHepMC(*it) );
	}
		
}

/**************************************************************************************************
This routine read from the data file necessary information for Tauspinner to work in mode with 
two reconstructed jets. This information is then passed to the user main program.
Routine  is not expected to be universal, but a proptotype. It assumes that hard event is written as:
   parton parton --> tau tau jet jet, tau --> decay products,
there may be intermediate states (bosons) predecessing tau-pair and or jet-jet pair too.
Incoming partons and intermediate bosons will not be used, exept, may be pdg_id of Higgs.
**************************************************************************************************/ 
int read_particles_for_VBF(HepMC::IO_GenEvent     &file,
                           SimpleParticle         &p1,
                           SimpleParticle         &p2,
                           SimpleParticle         &X,
                           SimpleParticle         &p3,
                           SimpleParticle         &p4,
                           SimpleParticle         &tau1,
                           SimpleParticle         &tau2,
                           vector<SimpleParticle> &tau1_daughters,
                           vector<SimpleParticle> &tau2_daughters) {

    HepMC::GenParticle *boson  = NULL;
    HepMC::GenParticle *beam1  = NULL, *beam2  = NULL;
    HepMC::GenParticle *jet1   = NULL, *jet2   = NULL;
    HepMC::GenParticle *tau1_h = NULL, *tau2_h = NULL;

    // Get next event from file
    HepMC::GenEvent *evt = new HepMC:: GenEvent();

    file.fill_next_event(evt);
    // if necessary convert units to GeV 
    if(evt->momentum_unit() !=HepMC::Units::GEV){
      evt->use_units(HepMC::Units::GEV,HepMC::Units::MM);
    }

    // Run TAUOLA on the event to decay TAUS because this example is prepared
    // for madgraph files without decays
    Tauolapp::TauolaHepMCEvent * t_event =  new Tauolapp::TauolaHepMCEvent(evt);

    //  we may want to undecay them first, but as default we do not do that.
    // t_event->undecayTaus();
    t_event->decayTaus();
     delete t_event; 

    if( file.rdstate() ) {
        delete evt;
        return 1;
    }



    // Find all particles needed to compute the weight
    bool ok = true;

    ok &= findBeams(evt, beam1, beam2);
    ok &= findBosonAndJets(beam1, boson, jet1, jet2);
      if(boson) ok &= findTaus(boson, tau1_h, tau2_h);  
      else      ok &= findTaus(beam1, tau1_h, tau2_h);
    //    ok &= findTaus(boson, tau1_h, tau2_h);

    if( !ok ) {
        evt->print();
        delete evt;
        return 1;
    }

    // Now fill simple particles
    p1   = SimpleParticleFromHepMC(beam1);
    p2   = SimpleParticleFromHepMC(beam2);
    if(boson) X    = SimpleParticleFromHepMC(boson);
    // X    = SimpleParticleFromHepMC(boson);
    p3   = SimpleParticleFromHepMC(jet1);
    p4   = SimpleParticleFromHepMC(jet2);
    tau1 = SimpleParticleFromHepMC(tau1_h);
    tau2 = SimpleParticleFromHepMC(tau2_h);
	
	/* chose beetween getDaughters and getDaughters_nonrecursive */
	
//    getDaughters_nonrecursive(tau1_h,tau1_daughters);
//    getDaughters_nonrecursive(tau2_h,tau2_daughters);
	tau1_daughters = *getDaughters(tau1_h);
	tau2_daughters = *getDaughters(tau2_h);	

    delete evt;
    return 0;
}
