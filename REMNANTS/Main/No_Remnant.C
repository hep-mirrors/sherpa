#include "REMNANTS/Main/No_Remnant.H"
#include "ATOOLS/Org/Exception.H"

using namespace REMNANTS;
using namespace ATOOLS;

No_Remnant::No_Remnant(const unsigned int _m_beam):
  Remnant_Base(rtp::intact,_m_beam) { }

bool No_Remnant::FillBlob(ParticleMomMap *ktmap,const bool & copy) {
  if (m_extracted.size()==0) {
    THROW(critical_error,"No particles extracted from intact beam.");
  }
  else if (m_extracted.size()>1) {
    THROW(critical_error,"Too many particles extracted from intact beam.");
  }
  p_beamblob->AddToOutParticles(*m_extracted.begin());
  return true;
}

bool No_Remnant::TestExtract(const Flavour &flav,const Vec4D &mom) {
  if (mom[0]>p_beam->OutMomentum()[0])
    return false;
  return true;
}
