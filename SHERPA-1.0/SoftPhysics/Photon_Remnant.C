#include "Photon_Remnant.H"

using namespace SHERPA;

Photon_Remnant::Photon_Remnant(const unsigned int _m_beam):
  Remnant_Base(Photon_Remnant::Photon,_m_beam) {}

bool Photon_Remnant::FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist)
{
  if (p_partner==NULL) {
    ATOOLS::msg.Error()<<"Photon_Remnant::FillBlob(..): "
		       <<"Partner Remnant not set! Abort."<<std::endl;
    exit(129);
  }
  // fill blob
  ATOOLS::Vec4D ptot=m_pbeam;
  for (size_t j=0;j<m_parton[1].size();++j) {
    beamblob->AddToOutParticles(m_parton[1][j]);
    if (particlelist!=NULL) {
      m_parton[1][j]->SetNumber(particlelist->size());
      particlelist->push_back(m_parton[1][j]);
    }
  }
  return true;
}
