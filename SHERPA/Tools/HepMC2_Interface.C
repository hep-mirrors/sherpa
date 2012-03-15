#include "SHERPA/Tools/HepMC2_Interface.H"

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/PdfInfo.h"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__HEPMC2__UNITS
#include "HepMC/Units.h"
#endif

using namespace SHERPA;
using namespace ATOOLS;

HepMC2_Interface::HepMC2_Interface():
  p_event(new HepMC::GenEvent())
{
}

HepMC2_Interface::~HepMC2_Interface()
{
  delete p_event;
  DeleteGenSubEventList();
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs,
                                    HepMC::GenEvent& event, double weight)
{
#ifdef USING__HEPMC2__UNITS
  event.use_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
#endif
  event.set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
  m_blob2genvertex.clear();
  m_particle2genparticle.clear();
  HepMC::GenVertex * vertex;
  std::vector<HepMC::GenParticle*> beamparticles;
  for (ATOOLS::Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    if (Sherpa2HepMC(*(blit),vertex)) {
      event.add_vertex(vertex);
      if ((*blit)->Type()==ATOOLS::btp::Signal_Process) {
        if ((**blit)["NLO_subeventlist"]) {
          THROW(fatal_error,"Events containing correlated subtraction events"
                +std::string(" cannot be translated into the full HepMC event")
                +std::string(" format.\n")
                +std::string("   Try 'EVENT_MODE=HepMC_Short' instead."));
        }
        event.set_signal_process_vertex(vertex);
        if((*blit)->NInP()==2) {
          kf_code fl1=(*blit)->InParticle(0)->Flav().HepEvt();
          kf_code fl2=(*blit)->InParticle(1)->Flav().HepEvt();
          double x1=(*blit)->InParticle(0)->Momentum()[0]/rpa->gen.PBeam(0)[0];
          double x2=(*blit)->InParticle(1)->Momentum()[0]/rpa->gen.PBeam(1)[0];
          double q(0.0), p1(0.0), p2(0.0);
          Blob_Data_Base *facscale((**blit)["Factorisation_Scale"]);
          if (facscale) q=sqrt(facscale->Get<double>());
          Blob_Data_Base *xf1((**blit)["XF1"]);
          Blob_Data_Base *xf2((**blit)["XF2"]);
          if (xf1) p1=xf1->Get<double>();
          if (xf1) p2=xf2->Get<double>();
          HepMC::PdfInfo pdfinfo(fl1, fl2, x1, x2, q, p1, p2);
          event.set_pdf_info(pdfinfo);
        }
      }
      else if ((*blit)->Type()==ATOOLS::btp::Beam) {
        if (vertex->particles_in_size()==1) {
          beamparticles.push_back(*vertex->particles_in_const_begin());
        }
      }
    }
  }
  if (beamparticles.size()==2) {
    event.set_beam_particles(beamparticles[0],beamparticles[1]);
  }
  std::vector<double> weights; weights.push_back(weight);
  event.weights()=weights;
  return true;
}

bool HepMC2_Interface::Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs,
                                         HepMC::GenEvent& event, double weight)
{
#ifdef USING__HEPMC2__UNITS
  event.use_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
#endif
  event.set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
  HepMC::GenVertex * vertex=new HepMC::GenVertex();
  std::vector<HepMC::GenParticle*> beamparticles;
  std::vector<std::pair<HepMC::FourVector,int> > beamparts,
                                                 remnantparts1, remnantparts2;
  Blob *sp(blobs->FindFirst(btp::Signal_Process));
  NLO_subevtlist* subevtlist(NULL);
  if (sp) {
    Blob_Data_Base * bdb((*sp)["NLO_subeventlist"]);
    if (bdb) subevtlist=bdb->Get<NLO_subevtlist*>();
  }
  for (ATOOLS::Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    Blob* blob=*blit;
    for (int i=0;i<blob->NInP();i++) {
      if (blob->InParticle(i)->ProductionBlob()==NULL) {
        Particle* parton=blob->InParticle(i);
        ATOOLS::Vec4D mom  = parton->Momentum();
        HepMC::FourVector momentum(mom[1],mom[2],mom[3],mom[0]);
        HepMC::GenParticle* inpart = new HepMC::GenParticle(momentum,parton->Flav().HepEvt(),2);
        vertex->add_particle_in(inpart);
        // distinct because SHRIMPS has no bunches for some reason
        if (blob->Type()==btp::Beam || blob->Type()==btp::Bunch) {
          beamparticles.push_back(inpart);
          beamparts.push_back(std::make_pair(momentum,parton->Flav().HepEvt()));
        }
      }
    }
    for (int i=0;i<blob->NOutP();i++) {
      if (blob->OutParticle(i)->DecayBlob()==NULL) {
        Particle* parton=blob->OutParticle(i);
        ATOOLS::Vec4D mom  = parton->Momentum();
        HepMC::FourVector momentum(mom[1],mom[2],mom[3],mom[0]);
        HepMC::GenParticle* outpart = new HepMC::GenParticle(momentum,parton->Flav().HepEvt(),1);
        vertex->add_particle_out(outpart);
        if (blob->Type()==btp::Beam) {
          if      (mom[3]>0) remnantparts1.push_back(std::make_pair(momentum,parton->Flav().HepEvt()));
          else if (mom[3]<0) remnantparts2.push_back(std::make_pair(momentum,parton->Flav().HepEvt()));
          else THROW(fatal_error,"Ill defined beam remnants.");
        }
      }
    }

    if ((*blit)->Type()==ATOOLS::btp::Signal_Process) {
      if((*blit)->NInP()==2) {
        kf_code fl1=(*blit)->InParticle(0)->Flav().HepEvt();
        kf_code fl2=(*blit)->InParticle(1)->Flav().HepEvt();
        double x1=(*blit)->InParticle(0)->Momentum()[0]/rpa->gen.PBeam(0)[0];
        double x2=(*blit)->InParticle(1)->Momentum()[0]/rpa->gen.PBeam(1)[0];
        double q(0.0), p1(0.0), p2(0.0);
        Blob_Data_Base *facscale((**blit)["Factorisation_Scale"]);
        if (facscale) q=sqrt(facscale->Get<double>());
        Blob_Data_Base *xf1((**blit)["XF1"]);
        Blob_Data_Base *xf2((**blit)["XF2"]);
        if (xf1) p1=xf1->Get<double>();
        if (xf1) p2=xf2->Get<double>();
        HepMC::PdfInfo pdfinfo(fl1, fl2, x1, x2, q, p1, p2);
        event.set_pdf_info(pdfinfo);
      }
    }
  }
  event.add_vertex(vertex);
  if (beamparticles.size()==2) {
    event.set_beam_particles(beamparticles[0],beamparticles[1]);
  }
  std::vector<double> weights; weights.push_back(weight);
  event.weights()=weights;

  if (subevtlist) {
    // build GenEvent for all subevts (where only the signal is available)
    // use stored beam & remnant particles from above
    // sub->m_flip==1 -> momenta need to be inverted for analyses
    for (size_t i(0);i<subevtlist->size();++i) {
      NLO_subevt * sub((*subevtlist)[i]);
      if (sub->m_result==0.) continue;
      HepMC::GenVertex * subvertex(new HepMC::GenVertex());
      HepMC::GenEvent * subevent(new HepMC::GenEvent());
      // assume that only 2->(n-2) processes
      HepMC::GenParticle *beam[2];
      if (beamparts.size()==2) {
        for (size_t j(0);j<beamparts.size();++j) {
          beam[j] = new HepMC::GenParticle(beamparts[j].first,
                                           beamparts[j].second,2);
          subvertex->add_particle_in(beam[j]);
        }
      }
      else {
        const ATOOLS::Vec4D *mom(sub->p_mom);
        const ATOOLS::Flavour *fl(sub->p_fl);
        for (size_t j(0);j<2;++j) {
          HepMC::FourVector momentum;
          if (sub->m_flip) momentum.set(-mom[j][1],-mom[j][2],-mom[j][3],mom[j][0]);
          else             momentum.set( mom[j][1], mom[j][2], mom[j][3],mom[j][0]);
          ATOOLS::Flavour flc(fl[j]);
          HepMC::GenParticle* inpart
              = new HepMC::GenParticle(momentum,flc.HepEvt(),1);
          subvertex->add_particle_in(inpart);
        }
      }
      const ATOOLS::Vec4D *mom(sub->p_mom);
      const ATOOLS::Flavour *fl(sub->p_fl);
      for (size_t j(2);j<(*subevtlist)[i]->m_n;++j) {
        HepMC::FourVector momentum;
        if (sub->m_flip) momentum.set(-mom[j][1],-mom[j][2],-mom[j][3],mom[j][0]);
        else             momentum.set( mom[j][1], mom[j][2], mom[j][3],mom[j][0]);
        ATOOLS::Flavour flc(fl[j]);
        HepMC::GenParticle* outpart
            = new HepMC::GenParticle(momentum,flc.HepEvt(),1);
        subvertex->add_particle_out(outpart);
      }
      // if beamremnants are present:
      // scale beam remnants of real event for energy momentum conservation :
      //   flavours might not add up properly for sub events,
      //   but who cares. they go down the beam pipe.
      // also assume they are all massless :
      //   this will give momentum conservation violations
      //   which are collider dependent only
      //
      //   for real event (as constructed in BRH) :
      //
      //     P = p + \sum r_i  with  P^2 = m^2  and  r_i^2 = p^2 = 0
      //
      //     further  P  = ( P^0 , 0 , 0 , P^z )
      //              p  = (  p  , 0 , 0 ,  p  )
      //             r_i = ( r_i , 0 , 0 , r_i )
      //
      //     => P^0 = p + \sum r_i
      //        P^z = \sqrt{(P^0)^2 - m^2} <  p + \sum r_i
      //
      // in a mass-symmetric collider, the excess is the same for both beams
      // ensuring momentum conservation
      //
      //   for sub event (constructed here):
      //
      //     P = p~ + \sum r_i~ = p~ + \sum u*r_i
      //
      //     where u = ( P^0 - p^0~ ) / \sum r_i^0
      //
      //     again, the r_i~ = u*r_i are constructed such that
      //
      //     => P^0 = p~ + \sum u*r_i
      //        P^z = \sqrt{(P^0)^2 - m^2} <  p~ + \sum u*r_i =  p + \sum r_i
      //
      //     leading to the same momentum conservation violations per beam
      //
      if (remnantparts1.size()!=0 && remnantparts2.size()!=0) {
        double res1(0.),res2(0.);
        for (size_t j(0);j<remnantparts1.size();++j) {
          res1+=remnantparts1[j].first.e();
        }
        for (size_t j(0);j<remnantparts2.size();++j) {
          res2+=remnantparts2[j].first.e();
        }
        ATOOLS::Vec4D hardparton[2];
        for (size_t j(0);j<2;++j) {
          if (sub->m_flip) hardparton[j]=ATOOLS::Vec4D(mom[j][0],Vec3D(-mom[j]));
          else             hardparton[j]=ATOOLS::Vec4D(mom[j][0],Vec3D( mom[j]));
        }
        // incoming partons might need to be flipped due to particle sorting
        bool flip(hardparton[0][3]<0);
        double u1((beamparts[0].first.e()-hardparton[flip?1:0][0])/res1);
        double u2((beamparts[1].first.e()-hardparton[flip?0:1][0])/res2);
        // filling
        for (size_t j(0);j<remnantparts1.size();++j) {
          HepMC::FourVector momentum(u1*remnantparts1[j].first.px(),
                                     u1*remnantparts1[j].first.py(),
                                     u1*remnantparts1[j].first.pz(),
                                     u1*remnantparts1[j].first.e());
          HepMC::GenParticle* outpart
              = new HepMC::GenParticle(momentum,remnantparts1[j].second,1);
          subvertex->add_particle_out(outpart);
        }
        for (size_t j(0);j<remnantparts2.size();++j) {
          HepMC::FourVector momentum(u2*remnantparts2[j].first.px(),
                                     u2*remnantparts2[j].first.py(),
                                     u2*remnantparts2[j].first.pz(),
                                     u2*remnantparts2[j].first.e());
          HepMC::GenParticle* outpart
              = new HepMC::GenParticle(momentum,remnantparts2[j].second,1);
          subvertex->add_particle_out(outpart);
        }
        if (beamparticles.size()==2) {
          subevent->set_beam_particles(beam[0],beam[1]);
        }
      }
      subevent->add_vertex(subvertex);
      // not enough info in subevents to set PDFInfo properly,
      // so leave it blank
      HepMC::PdfInfo subpdfinfo;
      subevent->set_pdf_info(subpdfinfo);
      // add weight
      std::vector<double> subweight; subweight.push_back(sub->m_result);
      subevent->weights()=subweight;
      // set the event number (could be used to identify correlated events)
      subevent->set_event_number(ATOOLS::rpa->gen.NumberOfGeneratedEvents());
      m_subeventlist.push_back(subevent);
    }
  }
  return true;
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs,
                                    double weight)
{
  if (blobs->empty()) {
    msg_Error()<<"Error in "<<METHOD<<"."<<std::endl
               <<"   Empty list - nothing to translate into HepMC."<<std::endl
               <<"   Continue run ... ."<<std::endl;
    return true;
  }
  if (p_event!=NULL) delete p_event;
  for (size_t i=0; i<m_subeventlist.size();++i)
    delete m_subeventlist[i];
  m_subeventlist.clear();
  p_event = new HepMC::GenEvent();
  return Sherpa2HepMC(blobs, *p_event, weight);
}

bool HepMC2_Interface::Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs,
                                         double weight)
{
  if (blobs->empty()) {
    msg_Error()<<"Error in "<<METHOD<<"."<<std::endl
               <<"   Empty list - nothing to translate into HepMC."<<std::endl
               <<"   Continue run ... ."<<std::endl;
    return true;
  }
  if (p_event!=NULL) delete p_event;
  DeleteGenSubEventList();
  p_event = new HepMC::GenEvent();
  return Sherpa2ShortHepMC(blobs, *p_event, weight);
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob * blob, HepMC::GenVertex *& vertex)
{
  if (m_ignoreblobs.count(blob->Type())) return false;
  int count = m_blob2genvertex.count(blob);
  if (count>0) {
    vertex = m_blob2genvertex[blob];
    return true;
  }
  else {
    ATOOLS::Vec4D pos = blob->Position();
    HepMC::FourVector position(pos[1],pos[2],pos[3],pos[0]);
    vertex = new HepMC::GenVertex(position,blob->Id());
    vertex->weights().push_back(1.);
  }

  bool okay = 1;
  HepMC::GenParticle * _particle;
  for (int i=0;i<blob->NInP();i++) {
    if (Sherpa2HepMC(blob->InParticle(i),_particle)) {
      vertex->add_particle_in(_particle);
    }
    else okay = 0;
  }
  for (int i=0;i<blob->NOutP();i++) {
    if (Sherpa2HepMC(blob->OutParticle(i),_particle)) {
      vertex->add_particle_out(_particle);
    }
    else okay = 0;
  }
  m_blob2genvertex.insert(std::make_pair(blob,vertex));
  if (!okay) {
    msg_Error()<<"Error in HepMC2_Interface::Sherpa2HepMC(Blob,Vertex)."<<std::endl
		       <<"   Continue event generation with new event."<<std::endl;
  }
  if (msg_LevelIsDebugging()) {
    ATOOLS::Vec4D check = blob->CheckMomentumConservation();
    double test         = ATOOLS::Vec3D(check).Abs();
    if (ATOOLS::dabs(1.-vertex->check_momentum_conservation()/test)>1.e-5 &&
	ATOOLS::dabs(test)>1.e-5)
      {
	msg_Error()<<"ERROR in "<<METHOD<<std::endl
			   <<"   Momentum not conserved. Continue."<<std::endl
			   <<"ERROR in Blob -> Vertex : "<<vertex->check_momentum_conservation()
			   <<" <- "<<test<<" "<<check
			   <<std::endl<<(*blob)<<std::endl;
	vertex->print(std::cout);
	std::cout<<"--------------------------------------------------------"<<std::endl;
      }
  }
  return okay;
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Particle * parton,HepMC::GenParticle *& particle)
{
  int count = m_particle2genparticle.count(parton);
  if (count>0) {
    particle = m_particle2genparticle[parton];
    return true;
  }

  ATOOLS::Vec4D mom  = parton->Momentum();
  HepMC::FourVector momentum(mom[1],mom[2],mom[3],mom[0]);
  int status=11;
  if (parton->DecayBlob()==NULL ||
      m_ignoreblobs.count(parton->DecayBlob()->Type())!=0) {
    status=1;
  }
  else {
    if (parton->DecayBlob()->Type()==ATOOLS::btp::Hadron_Decay ||
        parton->DecayBlob()->Type()==ATOOLS::btp::Hadron_Mixing) {
      status=2;
    }
    else if (parton->DecayBlob()->Type()==ATOOLS::btp::Signal_Process ||
             (parton->ProductionBlob() &&
              parton->ProductionBlob()->Type()==ATOOLS::btp::Signal_Process)) {
      status=3;
    }
    else if (parton->DecayBlob()->Type()==ATOOLS::btp::Bunch) {
      status=4;
    }
    else {
      status=11;
    }
  }
  particle = new HepMC::GenParticle(momentum,parton->Flav().HepEvt(),status);
  for (int i=1;i<3;i++) {
    if (parton->GetFlow(i)>0) particle->set_flow(i,parton->GetFlow(i));
  }
  m_particle2genparticle.insert(std::make_pair(parton,particle));
  return true;
}

void HepMC2_Interface::DeleteGenSubEventList()
{
  for (size_t i=0; i<m_subeventlist.size();++i)
    delete m_subeventlist[i];
  m_subeventlist.clear();
}
