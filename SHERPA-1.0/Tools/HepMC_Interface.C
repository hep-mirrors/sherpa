#ifdef _USE_HEPMC_


#include "HepMC_Interface.H"
#include "Vector.H"
#include "Message.H"

#include "CLHEP/Vector/LorentzVector.h"

using namespace SHERPA;
using namespace ATOOLS;
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
      if (flav.Bar()!=flav) {
	flav     = flav.Bar();
	pdata    = new ParticleData(flav.Name(),-flav.Kfcode(),
				    flav.Charge(),flav.PSMass(),lifetime,flav.Spin());
	m_flavour2particle[flav] = pdata;
	if (!p_particledatatable->insert(pdata)) msg.Error()<<"Error in HepMC_Interface::InitTheMap()."<<endl
							    <<"   Did not succeed to fill in "<<flav.Name()<<endl
							    <<"   Continue and hope for the best."<<endl;
      }
    }
  }
  msg.Debugging()<<"#######################################################"<<endl;
  p_particledatatable->print();
  msg.Debugging()<<"#######################################################"<<endl;
}

void HepMC_Interface::ResetTheLinks()
{
  m_blob2vertex.clear();
  m_particle2particle.clear();
}

bool HepMC_Interface::Sherpa2HepMC(Blob_List * _blobs,GenEvent *& _event)
{
  ResetTheLinks();
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
    vertex = NULL;
    if (Sherpa2HepMC(*(blit),vertex)) {
      _event->add_vertex(vertex);
      type = (*blit)->Type();
      if (type.find(string("Signal Process"))>-1) _event->set_signal_process_vertex(vertex);
    }
  }
}

bool HepMC_Interface::Sherpa2HepMC(Blob * _blob,GenVertex *& _vertex) 
{
  if (_blob->NInP()==1 && _blob->NOutP()==1) {
    if (_blob->InParticle(0)->Number()==_blob->OutParticle(0)->Number()) return 0;
  }

  int count = m_blob2vertex.count(_blob->Id());
  if (count>0) {
    _vertex = m_blob2vertex[_blob->Id()];
  }
  else {
    Vec4D pos = _blob->Position();
    HepLorentzVector position(pos[1],pos[2],pos[3],pos[0]);
    _vertex = new GenVertex(position,-_blob->Id()-1);
    _vertex->weights().push_back(1.);
  }

  bool take, okay = 1;
  GenParticle * _particle;
  Particle *    _SHparticle;

  for (int i=0;i<_blob->NInP();i++) {
    _SHparticle = _blob->InParticle(i);
    take        = 1;
    for (int j=0;j<_blob->NOutP();j++) {
      if (_SHparticle->Number()==_blob->OutParticle(j)->Number()) take = 0;
    }
    if (take) {
      if (Sherpa2HepMC(_SHparticle,_particle)) {
	_vertex->add_particle_in(_particle);
      }
      else okay = 0;
    }
  }
  for (int i=0;i<_blob->NOutP();i++) {
    _SHparticle = _blob->OutParticle(i);
    take        = 1;
    for (int j=0;j<_blob->NInP();j++) {
      if (_SHparticle->Number()==_blob->InParticle(j)->Number()) take = 0;
    }
    if (take) {
      if (Sherpa2HepMC(_blob->OutParticle(i),_particle)) {
	_vertex->add_particle_out(_particle);
      }
      else okay = 0;
    }
  }

  if (_vertex->particles_in_size()>0 || _vertex->particles_out_size()>0)
    m_blob2vertex.insert(std::make_pair(_blob->Id(),_vertex));
  else {
    delete _vertex;
    return 0;
  }

  if (!okay) {
    msg.Error()<<"Error in HepMC_Interface::Sherpa2HepMC(Blob,Vertex)."<<endl
	       <<"   Continue event generation with new event."<<endl;
  }

  return okay;
}

bool HepMC_Interface::Sherpa2HepMC(Particle * _SHparticle,GenParticle *& _particle) 
{
  int count = m_particle2particle.count(_SHparticle->Number());
  if (count>0) {
    _particle = m_particle2particle[_SHparticle->Number()];
    return 1;
  }

  Vec4D mom  = _SHparticle->Momentum();
  HepLorentzVector momentum(mom[1],mom[2],mom[3],mom[0]);
  int pdg_id = m_flavour2particle[_SHparticle->Flav()]->pdg_id();
  _particle  = new GenParticle(momentum,pdg_id,_SHparticle->Status());
  _particle->suggest_barcode(_SHparticle->Number());
  for (int i=1;i<3;i++) {
    if (_SHparticle->GetFlow(i)>0) _particle->set_flow(i,_SHparticle->GetFlow(i));
  }
  m_particle2particle.insert(std::make_pair(_SHparticle->Number(),_particle));
  return 1;
}

#endif
