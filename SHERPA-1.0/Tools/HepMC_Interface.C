#include "HepMC_Interface.H"
#include "Vector.H"
#include "Message.H"

#include "CLHEP/Vector/LorentzVector.h"

using namespace SHERPA;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace HepMC;

HepMC_Interface::HepMC_Interface()
{
  InitTheMap();
}

void HepMC_Interface::InitTheMap() 
{
  p_particledatatable = new ParticleDataTable;
  ParticleData * pdata;
  Fl_Iter fli;
  double lifetime;
  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav=fli.next()) {
    if (flav.IsOn() && flav.Size()==1) {
      if (flav.IsStable()) lifetime = -1;
                      else lifetime = clifetime_from_width(flav.Width());
      if (flav.IsSusy()) msg.Error()<<"Potential error in HepMC_Interface::InitTheMap()."<<endl
				    <<"   Kf codes of PDG and Sherpa do not coincide for Susys !"<<endl
				    <<"   Continue and hope for the best."<<endl;
      pdata    = new ParticleData(flav.Name(),flav.Kfcode(),
				  flav.Charge(),flav.PSMass(),lifetime,flav.Spin());
      m_flavour2particle[flav] = pdata;
      if (!p_particledatatable->insert(pdata)) msg.Error()<<"Error in HepMC_Interface::InitTheMap()."<<endl
							  <<"   Did not succeed to fill in "<<flav.Name()<<endl
							  <<"   Continue and hope for the best."<<endl;
      flav     = flav.Bar();
      pdata    = new ParticleData(flav.Name(),-flav.Kfcode(),
				  flav.Charge(),flav.PSMass(),lifetime,flav.Spin());
      m_flavour2particle[flav] = pdata;
      if (!p_particledatatable->insert(pdata)) msg.Error()<<"Error in HepMC_Interface::InitTheMap()."<<endl
							  <<"   Did not succeed to fill in "<<flav.Name()<<endl
							  <<"   Continue and hope for the best."<<endl;
    }
  }
  msg.Debugging()<<"#######################################################"<<endl;
  p_particledatatable->print();
  msg.Debugging()<<"#######################################################"<<endl;
}

bool HepMC_Interface::Sherpa2HepMC(Blob_List * _blobs,GenEvent *& _event)
{
  if (_blobs->empty()) {
    msg.Error()<<"Error in HepMC_Interface::Sherpa2HepMC(Blob_List)."<<endl
	       <<"   Empty list - nothing to translate into HepMC standard."<<endl
	       <<"   Continue run ... ."<<endl;
    return 1;
  }
  if (_event==NULL) _event = new GenEvent();

  GenVertex * vertex;
  string type;
  for (Blob_Iterator blit=_blobs->begin();blit!=_blobs->end();++blit) {
    Sherpa2HepMC(*(blit),vertex);
    _event->add_vertex(vertex);
    type = (*blit)->Type();
    if (type.find(string("Signal Process"))>-1) _event->set_signal_process_vertex(vertex);
  }
}

bool HepMC_Interface::Sherpa2HepMC(Blob * _blob,GenVertex *& _vertex) 
{
  int count = m_blob2vertex.count(_blob->Id());
  if (count>0) {
    _vertex = m_blob2vertex[_blob->Id()];
  }
  else {
    Vec4D pos = _blob->Position();
    HepLorentzVector position(pos[1],pos[2],pos[3],pos[0]);
    _vertex = new GenVertex(position,_blob->Id());
    _vertex->weights().push_back(1.);
  }

  bool okay = 1;
  GenParticle * _particle;
  for (int i=0;i<_blob->NInP();i++) {
    if (Sherpa2HepMC(_blob->InParton(i),_particle)) {
      _vertex->add_particle_in(_particle);
    }
    else okay = 0;
  }
  for (int i=0;i<_blob->NOutP();i++) {
    if (Sherpa2HepMC(_blob->OutParton(i),_particle)) {
      _vertex->add_particle_out(_particle);
    }
    else okay = 0;
  }
  m_blob2vertex[_blob->Id()] = _vertex;
  if (!okay) {
    msg.Error()<<"Error in HepMC_Interface::Sherpa2HepMC(Blob,Vertex)."<<endl
	       <<"   Continue event generation with new event."<<endl;
  }
  return okay;
}

bool HepMC_Interface::Sherpa2HepMC(Parton * _parton,GenParticle *& _particle) 
{
  int count = m_parton2particle.count(_parton->Number());
  if (count>0) {
    _particle = m_parton2particle[_parton->Number()];
    return 1;
  }

  Vec4D mom  = _parton->Momentum();
  HepLorentzVector momentum(mom[1],mom[2],mom[3],mom[0]);
  int pdg_id = m_flavour2particle[_parton->Flav()]->pdg_id();
  _particle  = new GenParticle(momentum,pdg_id,_parton->Status());
  _particle->suggest_barcode(_parton->Number());
  for (int i=1;i<3;i++) {
    if (_parton->GetFlow(i)>0) _particle->set_flow(i,_parton->GetFlow(i));
  }
  m_parton2particle.insert(std::make_pair(_parton->Number(),_particle));
  return 1;
}

