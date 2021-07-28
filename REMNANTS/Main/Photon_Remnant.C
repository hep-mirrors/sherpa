#include "REMNANTS/Main/Photon_Remnant.H"
#include "ATOOLS/Org/Exception.H"

using namespace REMNANTS;

Photon_Remnant::Photon_Remnant(const unsigned int _m_beam):
  Remnant_Base(rtp::photon,_m_beam) {}

bool Photon_Remnant::FillBlob(ATOOLS::Blob *beamblob,
			      ATOOLS::Particle_List *particlelist)
{
  msg_Out()<<METHOD<<" not fully implemented yet.  Will exit.\n";
  exit(1);
  if (p_partner==NULL) {
    THROW(critical_error,"Partner Remnant not set.");
  }
  for (ATOOLS::Part_Iterator pmit=m_extracted.begin();
       pmit!=m_extracted.end();pmit++) {
    beamblob->AddToOutParticles(*pmit);
    if (particlelist!=NULL) {
      (*pmit)->SetNumber(particlelist->size());
      particlelist->push_back(*pmit);
    }
  }
  return true;
}

bool Photon_Remnant::AdjustKinematics()
{
  return true;
}
