#include "No_Remnant.H"

using namespace SHERPA;

No_Remnant::No_Remnant(const unsigned int _m_beam):
  Remnant_Base(No_Remnant::Intact,_m_beam) {}

bool No_Remnant::FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist)
{
  if (p_partner==NULL) {
    ATOOLS::msg.Error()<<"No_Remnant::FillBlob(..): "
		       <<"Partner Remnant not set! Abort."<<std::endl;
    exit(129);
  }
  // fill blob
  for (size_t j=0;j<m_parton[1].size();++j) {
    beamblob->AddToOutParticles(m_parton[1][j]);
    if (particlelist!=NULL) {
      m_parton[1][j]->SetNumber(particlelist->size());
      particlelist->push_back(m_parton[1][j]);
    }
  }
  return true;
}

bool No_Remnant::AdjustKinematics()
{
  return true;
}
