#include "REMNANTS/Main/No_Remnant.H"
#include "ATOOLS/Org/Exception.H"

using namespace REMNANTS;
using namespace ATOOLS;

No_Remnant::No_Remnant(const size_t & beam,const size_t & tag):
  Remnant_Base(Flavour(kf_none),beam,tag) { }

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
  if (flav!=p_beam->Bunch()) return false;
  for (size_t i=0;i<4;i++) {
    double diff = ((mom[i]-p_beam->OutMomentum(m_tag)[i])/
		   (mom[i]+p_beam->OutMomentum(m_tag)[i])); 
    if (diff>1.e-6) {
      msg_Error()<<"Error in "<<METHOD<<": difference in four-momenta = "
		 <<diff<<" = "<<mom<<" - "<<p_beam->OutMomentum()<<".\n";
      return false;
    }
  }
  return true;
}
