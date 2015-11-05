#include "SHERPA/Main/Filter.H"
#include "SHERPA/Tools/Output_Base.H"

using namespace SHERPA;
using namespace ATOOLS;


Filter::Filter() :
  m_N_charged_min(1), m_N_charged_max(15),
  m_etamax(2.5), m_ptmin(0.1)
{
  msg_Out()<<"In "<<METHOD<<".\n";
}

Filter::~Filter() {
  msg_Out()<<"In "<<METHOD<<".\n";
}

bool Filter::operator() (Blob_List * blobs) const {
  size_t N_charged(NumberOfChargedFSParticles(blobs));
  if (N_charged>=m_N_charged_min && N_charged<=m_N_charged_max) return true;
  return false;
}

size_t Filter::NumberOfChargedFSParticles(Blob_List * blobs) const {
  Particle_List particles(blobs->ExtractParticles(1,1));
  size_t Nparticles(0);
  for (Particle_List::iterator partiter=particles.begin();
       partiter!=particles.end();++partiter) {
    Particle * particle(*partiter);
    Vec4D      partmom(particle->Momentum());
    if (particle->DecayBlob())         continue;
    if (particle->Flav().Charge()==0)  continue;
    if (dabs(partmom.Eta())>m_etamax)  continue;
    if (dabs(partmom.PPerp())<m_ptmin) continue;
    msg_Out()<<(*particle)<<"\n";
    Nparticles++;
  }
  msg_Out()<<"  -------------------------------------------------\n";
  msg_Out()<<METHOD<<": "<<particles.size()<<"/"<<Nparticles<<"\n";
  msg_Out()<<"  -------------------------------------------------\n";
  return Nparticles;
}
