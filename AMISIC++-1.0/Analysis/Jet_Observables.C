#include "Jet_Observables.H"

using namespace AMISIC;

Jet_Observables::Jet_Observables(int _m_typekt,double _m_minkt,double _m_maxkt,int _m_nbinskt,
				 int _m_typept,double _m_minpt,double _m_maxpt,int _m_nbinspt,
				 int _m_typey,double _m_miny,double _m_maxy,int _m_nbinsy,
				 int _m_jets,int _m_maxmult,double _m_ycut,double _m_ptcut):
  Primitive_Observable_Base(),
//   p_jetfinder(new ATOOLS::Jet_Algorithm(_m_ycut,4)),
  p_primkt(new Primordial_KT(_m_typekt,_m_minkt,_m_maxkt,_m_nbinskt)),
  p_jetmulti(new Jet_Multiplicity(0,0.0,(double)_m_maxmult,_m_maxmult,_m_ycut,_m_ptcut)), 
  p_jetpt(new Jet_PT(_m_typept,_m_minpt,_m_maxpt,_m_nbinspt,_m_jets)),
  p_jety(new Jet_Y(_m_typey,_m_miny,_m_maxy,_m_nbinsy,_m_jets))
{
//   p_jety->SetJetFinder(p_jetfinder);
//   p_jetpt->SetJetFinder(p_jetfinder);
//   p_jetmulti->SetJetFinder(p_jetfinder);
  name=std::string("jet_observables");
}

Jet_Observables::~Jet_Observables()
{
  delete p_primkt;
  delete p_jetmulti;
  delete p_jetpt;
  delete p_jety;
//   delete p_jetfinder;
}

void Jet_Observables::Evaluate(const ATOOLS::Particle_List &particles,double weight)
{
  p_primkt->Evaluate(particles,weight);
  p_jetmulti->Evaluate(particles,weight);
  p_jetpt->Evaluate(particles,weight);
  p_jety->Evaluate(particles,weight);
}

void Jet_Observables::Evaluate(const ATOOLS::Blob_List &blobs,double weight) 
{
  ATOOLS::Particle_List particles;
  for (ATOOLS::Blob_Const_Iterator bit=blobs.begin();bit!=blobs.end();++bit) { 
    for (int i=0;i<(*bit)->NOutP();++i) {
      ATOOLS::Particle *cur=(*bit)->OutParticle(i);
      if (cur->Status()==1) { 
	particles.push_back(cur);
      }
    }
  }
  Evaluate(particles,weight);
}

void Jet_Observables::EndEvaluation() 
{
  p_jety->EndEvaluation(); 
  p_jetpt->EndEvaluation(); 
  p_jetmulti->EndEvaluation();
  p_primkt->EndEvaluation();
}

void Jet_Observables::Output(std::string pathname) 
{
  p_jety->Output(pathname);
  p_jetpt->Output(pathname);
  p_jetmulti->Output(pathname);
  p_primkt->Output(pathname);
}
