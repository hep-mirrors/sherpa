#include "Electron_Remnant.H"

using namespace SHERPA;

Electron_Remnant::Electron_Remnant(PDF::ISR_Handler *isrhandler,const unsigned int _m_beam,double _m_scale):
  Remnant_Base(Electron_Remnant::Electron,_m_beam),
  m_scale(-_m_scale)
{
  if (isrhandler==NULL) {
    ATOOLS::msg.Error()<<"Electron_Remnant::Electron_Remnant(NULL,"<<m_beam<<"): "
		       <<"Cannot proceed without ISR and Beam Handler! Abort."<<std::endl;
    exit(129);
  }
  p_pdfbase=isrhandler->PDF(m_beam);
}

bool Electron_Remnant::FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist)
{
  if (p_partner==NULL) {
    ATOOLS::msg.Error()<<"Electron_Remnant::FillBlob(..): "
		       <<"Partner Remnant not set! Abort."<<std::endl;
    exit(129);
  }
  // p_last remains NULL since in DIS the hadron accounts for momentum conservation
  p_beamblob=beamblob;
  m_pbeam=beamblob->InParticle(0)->Momentum();
  // fill blob
  ATOOLS::Vec4D ptot=m_pbeam;
  for (size_t j=0;j<m_parton[1].size();++j) {
    beamblob->AddToOutParticles(m_parton[1][j]);
    if (particlelist!=NULL) {
      m_parton[1][j]->SetNumber(particlelist->size());
      particlelist->push_back(m_parton[1][j]);
    }
    ptot=ptot-m_parton[1][j]->Momentum();
  }
  // eventually extract single photon
  if (!ATOOLS::IsZero(ptot[0])) {
    ATOOLS::Particle *rem = new ATOOLS::Particle(-1,ATOOLS::Flavour(ATOOLS::kf::photon),ptot);
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
