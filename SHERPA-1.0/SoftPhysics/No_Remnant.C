#include "No_Remnant.H"

#include "Exception.H"

using namespace SHERPA;

No_Remnant::No_Remnant(const unsigned int _m_beam):
  Remnant_Base(rtp::intact,_m_beam) {}

bool No_Remnant::FillBlob(ATOOLS::Blob *beamblob,
			  ATOOLS::Particle_List *particlelist)
{
  if (p_partner==NULL) {
    throw(ATOOLS::Exception(ATOOLS::ex::critical_error,
			    "Partner Remnant not set.",
			    "No_Remnant","FillBlob"));
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

bool No_Remnant::AdjustKinematics()
{
  return true;
}
