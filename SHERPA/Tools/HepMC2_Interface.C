#include "HepMC2_Interface.H"

#include "Blob_List.H"
#include "Particle.H"
#include "Vector.H"
#include "Run_Parameter.H"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/PdfInfo.h"
#include "CXXFLAGS_PACKAGES.H"
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
                                    HepMC::GenEvent& event)
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
          if ((*blit)->InParticle(0)->Flav().IsAnti()) fl1=-fl1;
          if ((*blit)->InParticle(0)->Flav().Kfcode()==kf_gluon) fl1=0;
          kf_code fl2=(*blit)->InParticle(1)->Flav().HepEvt();
          if ((*blit)->InParticle(1)->Flav().IsAnti()) fl2=-fl2;
          if ((*blit)->InParticle(1)->Flav().Kfcode()==kf_gluon) fl2=0;
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
        Blob_Data_Base *info((**blit)["ME_Weight"]);
        double weight=1.0;
        if (info!=NULL) weight=info->Get<double>();
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

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob_List *const blobs)
{
  if (blobs->empty()) {
    msg_Error()<<"Error in "<<METHOD<<"."<<std::endl
               <<"   Empty list - nothing to translate into HepMC."<<std::endl
               <<"   Continue run ... ."<<std::endl;
    return true;
  }
  if (p_event!=NULL) delete p_event;
  p_event = new HepMC::GenEvent();
  return Sherpa2HepMC(blobs, *p_event);
}

bool HepMC2_Interface::Sherpa2HepMC(ATOOLS::Blob * blob, HepMC::GenVertex *& vertex)
{
  if (blob->Type()==ATOOLS::btp::Bunch) return false;
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
  if (parton->DecayBlob()!=NULL) stat = 2;
               else if (stat==2) stat=1;
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

