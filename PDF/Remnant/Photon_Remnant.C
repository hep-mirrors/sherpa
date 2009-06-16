#include "PDF/Remnant/Photon_Remnant.H"

#include "ATOOLS/Org/Exception.H"

using namespace PDF;

Photon_Remnant::Photon_Remnant(const unsigned int _m_beam):
  Remnant_Base(rtp::photon,_m_beam) {}

bool Photon_Remnant::FillBlob(ATOOLS::Blob *beamblob,
			      ATOOLS::Particle_List *particlelist)
{
  if (p_partner==NULL) {
    THROW(critical_error,"Partner Remnant not set.");
  }
  for (size_t j=0;j<m_extracted.size();++j) {
    beamblob->AddToOutParticles(m_extracted[j]);
    if (particlelist!=NULL) {
      m_extracted[j]->SetNumber(particlelist->size());
      particlelist->push_back(m_extracted[j]);
    }
  }
  return true;
}

bool Photon_Remnant::AdjustKinematics()
{
  return true;
}
