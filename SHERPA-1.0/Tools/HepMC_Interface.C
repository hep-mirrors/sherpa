#include "HepMC_Interface.H"
#include "Vector.H"
#include "Message.H"
#include "Run_Parameter.H"

#ifdef CLHEP_SUPPORT
#include "CLHEP/Vector/LorentzVector.h"
#endif

namespace CLHEP { class Dummy {}; }

using namespace CLHEP;

using namespace SHERPA;
using namespace HepMC;
using namespace HepPDT;

HepMC_Interface::HepMC_Interface():
#ifdef CLHEP_SUPPORT
  p_event(new HepMC::GenEvent())
#else
  p_event(NULL)
#endif
{
  InitTheMap();
}

HepMC_Interface::~HepMC_Interface()
{
#ifdef CLHEP_SUPPORT
  delete p_event;
  delete p_particledatatable;
#endif
}

void HepMC_Interface::InitTheMap() 
{
#ifdef CLHEP_SUPPORT
  p_particledatatable = new DefaultConfig::ParticleDataTable();
  ATOOLS::Fl_Iter fli;
  HepPDT::TableBuilder *build = new HepPDT::TableBuilder(*p_particledatatable);
  for (ATOOLS::Flavour flav=fli.first();flav!=ATOOLS::Flavour(ATOOLS::kf::none);flav=fli.next()) {
    if (flav.IsOn() && flav.Size()==1) {
      HepPDT::TempParticleData& tpd1 = build->getParticleData(HepPDT::ParticleID(flav.HepEvt()));
      tpd1.tempParticleName=flav.TexName();
      tpd1.tempCharge=flav.Charge();
      tpd1.tempMass=HepPDT::Measurement(flav.Mass(),0.);
      tpd1.tempSpin=HepPDT::SpinState(flav.Spin(),0.,0.);
      tpd1.tempWidth=HepPDT::Measurement(flav.Width(),0.);
      build->addParticle(tpd1);
      flav=flav.Bar();
      HepPDT::TempParticleData& tpd2 = build->getParticleData(HepPDT::ParticleID(flav.HepEvt()));
      tpd2.tempParticleName=flav.TexName();
      tpd2.tempCharge=flav.Charge();
      tpd2.tempMass=HepPDT::Measurement(flav.Mass(),0.);
      tpd2.tempSpin=HepPDT::SpinState(flav.Spin(),0.,0.);
      tpd2.tempWidth=HepPDT::Measurement(flav.Width(),0.);
      build->addParticle(tpd2);
    }
  }
  delete build;
  p_particledatatable->writeParticleData(ATOOLS::msg.Debugging());
  ATOOLS::msg.Debugging()<<std::endl;
#endif
}

bool HepMC_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs)
{
#ifdef CLHEP_SUPPORT
  if (blobs->empty()) {
    ATOOLS::msg.Error()<<"Error in HepMC_Interface::Sherpa2HepMC(Blob_List)."<<std::endl
		       <<"   Empty list - nothing to translate into HepMC standard."<<std::endl
		       <<"   Continue run ... ."<<std::endl;
    return true;
  }
  if (p_event!=NULL) delete p_event;
  p_event = new HepMC::GenEvent();
  m_blob2vertex.clear();
  m_parton2particle.clear();
  GenVertex * vertex;
  std::string type;
  for (ATOOLS::Blob_Iterator blit=blobs->begin();blit!=blobs->end();++blit) {
    if (Sherpa2HepMC(*(blit),vertex)) {
      p_event->add_vertex(vertex);
      type = (*blit)->Type();
      if (type.find(std::string("Signal Process"))!=std::string::npos) 
	p_event->set_signal_process_vertex(vertex);
    }
  }
  if (ATOOLS::msg.Debugging()) {
    int charge = 0;
    for ( GenEvent::particle_const_iterator p = p_event->particles_begin();
	  p != p_event->particles_end(); ++p ) {
      if ((*p)->status()==1) charge += (*p)->particleID().threeCharge();
    }
    if (charge!=(ATOOLS::rpa.gen.Beam1().IntCharge()+ATOOLS::rpa.gen.Beam2().IntCharge())) {
      ATOOLS::msg.Error()<<"ERROR in HepMC_Interface::Sherpa2HepMC(ATOOLS::Blob_List *):"<<std::endl
			 <<"   Charge not conserved. Continue."<<std::endl;
    }
  }
#endif
  return true;
}

bool HepMC_Interface::Sherpa2HepMC(ATOOLS::Blob * blob,HepMC::GenVertex *& vertex) 
{
#ifdef CLHEP_SUPPORT
  if (blob->Type()==ATOOLS::btp::Bunch) return false;
  int count = m_blob2vertex.count(blob->Id());
  if (count>0) {
    vertex = m_blob2vertex[blob->Id()];
  }
  else {
    ATOOLS::Vec4D pos = blob->Position();
    HepLorentzVector position(pos[1],pos[2],pos[3],pos[0]);
    vertex = new GenVertex(position,blob->Id());
    vertex->weights().push_back(1.);
  }

  bool okay = 1;
  GenParticle * _particle;
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
  m_blob2vertex.insert(std::make_pair(blob->Id(),vertex));
  if (!okay) {
    ATOOLS::msg.Error()<<"Error in HepMC_Interface::Sherpa2HepMC(Blob,Vertex)."<<std::endl
		       <<"   Continue event generation with new event."<<std::endl;
  }
  if (ATOOLS::msg.Debugging()) {
    ATOOLS::Vec4D check = blob->CheckMomentumConservation();
    double test         = ATOOLS::Vec3D(check).Abs();
    if (ATOOLS::dabs(1.-vertex->check_momentum_conservation()/test)>1.e-6 &&
	ATOOLS::dabs(test)>1.e-6)
      {
	ATOOLS::msg.Error()<<"ERROR in HepMC_Interface::Sherpa2HepMC(ATOOLS::Blob_List *):"<<std::endl
			   <<"   Momentum not conserved. Continue."<<std::endl
			   <<"ERROR in Blob -> Vertex : "<<vertex->check_momentum_conservation()
			   <<" <- "<<test<<" "<<check
			   <<std::endl<<(*blob)<<std::endl;
	vertex->print(std::cout);
	std::cout<<"--------------------------------------------------------"<<std::endl;
      }
  }
  return okay;
#endif
  return true;
}

bool HepMC_Interface::Sherpa2HepMC(ATOOLS::Particle * parton,HepMC::GenParticle *& particle) 
{
#ifdef CLHEP_SUPPORT
  long int number = (long int)(parton), count = m_parton2particle.count(number);
  if (count>0) {
    particle = m_parton2particle[number];
    return true;
  }

  ATOOLS::Vec4D mom  = parton->Momentum();
  HepLorentzVector momentum(mom[1],mom[2],mom[3],mom[0]);
  int stat = parton->Status();
  if (parton->DecayBlob()!=NULL) stat = 2;
  else if (stat==2) stat=1;
  if (stat==2) {
    if (parton->DecayBlob()->Type()==ATOOLS::btp::Signal_Process ||
	parton->ProductionBlob()->Type()==ATOOLS::btp::Signal_Process) stat = 3;
  }
  particle = new GenParticle(momentum,parton->Flav(),stat);
  for (int i=1;i<3;i++) {
    if (parton->GetFlow(i)>0) particle->set_flow(i,parton->GetFlow(i));
  }
  m_parton2particle.insert(std::make_pair(number,particle));
#endif
  return true;
}

