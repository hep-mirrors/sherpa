#include "Electron_Remnant.H"

#include "Exception.H"

using namespace SHERPA;

Electron_Remnant::Electron_Remnant(PDF::ISR_Handler *isrhandler,
				   const unsigned int _m_beam):
  Remnant_Base(rtp::electron,_m_beam)
{
  if (isrhandler==NULL) {
    THROW(fatal_error,"Cannot proceed without ISR and Beam Handler.");
  }
  p_pdfbase=isrhandler->PDF(m_beam);
}

bool Electron_Remnant::FillBlob(ATOOLS::Blob *beamblob,
				ATOOLS::Particle_List *particlelist)
{
  if (p_partner==NULL) {
    THROW(critical_error,"Partner Remnant not set.");
  }
  p_beamblob=beamblob;
  m_pbeam=beamblob->InParticle(0)->Momentum();
  ATOOLS::Vec4D ptot=m_pbeam;
  for (size_t j=0;j<m_extracted.size();++j) {
    beamblob->AddToOutParticles(m_extracted[j]);
    if (particlelist!=NULL) {
      m_extracted[j]->SetNumber(particlelist->size());
      particlelist->push_back(m_extracted[j]);
    }
    ptot=ptot-m_extracted[j]->Momentum();
  }
  if (!ATOOLS::IsZero(ptot[0])) {
    ATOOLS::Particle *rem = 
      new ATOOLS::Particle(-1,ATOOLS::Flavour(ATOOLS::kf::photon),ptot);
    rem->SetNumber((long int)rem);
    rem->SetInfo('F');
    beamblob->AddToOutParticles(rem);
    p_last[0]=rem;
    if (particlelist!=NULL) {
      rem->SetNumber(particlelist->size());
      particlelist->push_back(rem);
    }
  }
  return true;
}

bool Electron_Remnant::AdjustKinematics()
{
  // if (p_partner->Type()!=Hadron) return Remnant_Base::AdjustKinematics();
  return true;
}

ATOOLS::Flavour Electron_Remnant::
ConstituentType(const ATOOLS::Flavour &flavour) 
{
  if (flavour.Kfcode()==ATOOLS::kf::e) return flavour;
  return ATOOLS::kf::none;
}
