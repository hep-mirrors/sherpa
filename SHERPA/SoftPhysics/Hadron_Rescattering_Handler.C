#include "SHERPA/SoftPhysics/Hadron_Rescattering_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <cmath>
#include <algorithm>

using namespace SHERPA;
using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

Hadron_Rescattering_Handler::Hadron_Rescattering_Handler() :
  m_name("Hadron_Rescattering"), m_on(true),
  m_rescattering(Hadron_Rescatterings(true)),
  m_deuterons(Deuteron_Coalescence( 3.,            //0.7,  b_coal in fm
				    0.197326/(3.), // p_coal in GeV from uncertainty principle : p = hbar/(2r)
				    CoalescenceModel::Exponential
				    // or HardDisk, choose as you wish
				    ))
{
  if (!m_on) return;
  m_rescattering.Initialize();
  m_particles.clear();
}

Hadron_Rescattering_Handler::~Hadron_Rescattering_Handler() {}

void Hadron_Rescattering_Handler::HarvestParticles(Blob * blob) {
  if (m_treatedblobs.find(blob)!=m_treatedblobs.end()) return;
  for (size_t i=0;i<blob->NOutP();i++)  {
    Particle * part1 = blob->OutParticle(i);
    if (m_particles.size()>0) {
      for (std::set<Particle * >::iterator pit=m_particles.begin();
	   pit!=m_particles.end();pit++) {
	if (part1==(*pit)) continue;
	Hadron_Collision * cc = new Hadron_Collision(part1,(*pit));
	if (!cc->On() || cc->B2()>10000.) { delete cc; continue; }
	m_candidates[cc->T_Lab()] = cc;
      }
    }
    m_particles.insert(part1);
  }
  msg_Out()<<"   * "<<blob->NOutP()<<" particles harvested from "<<blob->Id()<<", "
	   <<"now "<<m_candidates.size()<<" candidates and "
	   <<m_particles.size()<<" particles in list.\n";
  blob->UnsetStatus(blob_status::needs_hadronRescatter);
}

bool Hadron_Rescattering_Handler::operator()() {
  while (!m_candidates.empty()) {
    Hadron_Collision * candidate = m_candidates.begin()->second;
    // There is a subtlety here: we may need to have a competition between
    // coalescence and rescattering - maybe check coalescence first before
    // rescattering?  (this would be an implicit bias, but probably harmless).
    // We could make this switchable ...
    p_blob = m_rescattering(candidate);
    if (p_blob) {
      EliminateOptions(candidate);
      return true;
    }
    m_candidates.erase(m_candidates.begin());
  }
  return false;
}

void Hadron_Rescattering_Handler::EliminateOptions(Hadron_Collision * candidate) {
  msg_Out()<<METHOD<<" for "<<m_particles.size()<<" particles and "
	   <<m_candidates.size()<<" candidate scatters.\n";
  for (size_t i=0;i<2;i++) {
    Particle * part = candidate->InPart(i);
    set<Particle *>::iterator pit = m_particles.find(part);
    if (pit!=m_particles.end()) m_particles.erase(pit);
    map<double,Hadron_Collision *>::iterator cit=m_candidates.begin();
    while(cit!=m_candidates.end()) {
      if (cit->second->InPart(0)==part || cit->second->InPart(1)==part) {
	cit = m_candidates.erase(cit);
      }
      else cit++;
    }
  }
  msg_Out()<<METHOD<<" results in "<<m_particles.size()<<" particles and "
	   <<m_candidates.size()<<" candidate scatters.\n";
}

void Hadron_Rescattering_Handler::CleanUp(const size_t & mode) {
  m_treatedblobs.clear();
  m_particles.clear();
  for (auto &c : m_candidates) delete c.second;
  m_candidates.clear();
}

void Hadron_Rescattering_Handler::Output() {
  msg_Out() << "======================================================\n"
            << "Hadron_Rescattering_Handler::Output(): rescattering is "
            << (m_on ? "ON" : "OFF") << "\n"
            << "======================================================\n";
}

Blob* Hadron_Rescattering_Handler::Rescatter() {
  return nullptr;
}

Blob* Hadron_Rescattering_Handler::RescatterCoalescence() {
  return nullptr;
}

