#include "SHERPA/Tools/HepMC2_Interface.H"

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Run_Parameter.H"

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
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs,
                                    HepMC::GenEvent& event, double weight)
{
#ifdef USING__HEPMC2__UNITS
  event.use_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
#endif
  event.set_event_number(ATOOLS::rpa.gen.NumberOfDicedEvents());
  m_blob2genvertex.clear();
  m_particle2genparticle.clear();
  HepMC::GenVertex * vertex;
  std::vector<HepMC::GenParticle*> beamparticles;
  for (ATOOLS::Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    if (Sherpa2HepMC(*(blit),vertex)) {
      event.add_vertex(vertex);
      if ((*blit)->Type()==ATOOLS::btp::Signal_Process) {
        event.set_signal_process_vertex(vertex);
        if((*blit)->NInP()==2) {
          kf_code fl1=(*blit)->InParticle(0)->Flav().HepEvt();
          kf_code fl2=(*blit)->InParticle(1)->Flav().HepEvt();
          double x1=(*blit)->InParticle(0)->Momentum()[0]/rpa.gen.PBeam(0)[0];
          double x2=(*blit)->InParticle(1)->Momentum()[0]/rpa.gen.PBeam(1)[0];
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
        std::vector<double> weights; weights.push_back(weight);
        event.weights()=weights;
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
  return true;
}

bool HepMC2_Interface::Sherpa2ShortHepMC(ATOOLS::Blob_List *const blobs,
                                         HepMC::GenEvent& event, double weight)
{
#ifdef USING__HEPMC2__UNITS
  event.use_units(HepMC::Units::GEV,
                  HepMC::Units::MM);
#endif
  event.set_event_number(ATOOLS::rpa.gen.NumberOfDicedEvents());
  HepMC::GenVertex * vertex=new HepMC::GenVertex();
  std::vector<HepMC::GenParticle*> beamparticles;
  for (ATOOLS::Blob_List::iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    Blob* blob=*blit;
    for (int i=0;i<blob->NInP();i++) {
      if (blob->InParticle(i)->ProductionBlob()==NULL) {
        Particle* parton=blob->InParticle(i);
        ATOOLS::Vec4D mom  = parton->Momentum();
        HepMC::FourVector momentum(mom[1],mom[2],mom[3],mom[0]);
        HepMC::GenParticle* inpart = new HepMC::GenParticle(momentum,parton->Flav().HepEvt(),2);
        vertex->add_particle_in(inpart);
        beamparticles.push_back(inpart);
      }
    }
    for (int i=0;i<blob->NOutP();i++) {
      if (blob->OutParticle(i)->DecayBlob()==NULL) {
        Particle* parton=blob->OutParticle(i);
        ATOOLS::Vec4D mom  = parton->Momentum();
        HepMC::FourVector momentum(mom[1],mom[2],mom[3],mom[0]);
        vertex->add_particle_out(new HepMC::GenParticle(momentum,parton->Flav().HepEvt(),1));
      }
    }

    if ((*blit)->Type()==ATOOLS::btp::Signal_Process) {
      if((*blit)->NInP()==2) {
        kf_code fl1=(*blit)->InParticle(0)->Flav().HepEvt();
        kf_code fl2=(*blit)->InParticle(1)->Flav().HepEvt();
        double x1=(*blit)->InParticle(0)->Momentum()[0]/rpa.gen.PBeam(0)[0];
        double x2=(*blit)->InParticle(1)->Momentum()[0]/rpa.gen.PBeam(1)[0];
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
      std::vector<double> weights; weights.push_back(weight);
      event.weights()=weights;
    }
  }
  event.add_vertex(vertex);
  if (beamparticles.size()==2) {
    event.set_beam_particles(beamparticles[0],beamparticles[1]);
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
  p_event = new HepMC::GenEvent();
  return Sherpa2HepMC(blobs, *p_event, weight);
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob * blob, HepMC::GenVertex *& vertex)
{
  if (blob->Type()==ATOOLS::btp::Bunch) return false;
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
  int stat = int(parton->Status());
  if (parton->DecayBlob()!=NULL &&
      m_ignoreblobs.count(parton->DecayBlob()->Type())==0) stat = 2;
  else if (stat==2 || stat==4) stat=1;
  if (stat==2) {
    if (parton->DecayBlob()->Type()==ATOOLS::btp::Signal_Process ||
	parton->ProductionBlob()->Type()==ATOOLS::btp::Signal_Process) stat = 3;
  }
  particle = new HepMC::GenParticle(momentum,parton->Flav().HepEvt(),stat);
  for (int i=1;i<3;i++) {
    if (parton->GetFlow(i)>0) particle->set_flow(i,parton->GetFlow(i));
  }
  m_particle2genparticle.insert(std::make_pair(parton,particle));
  return true;
}

