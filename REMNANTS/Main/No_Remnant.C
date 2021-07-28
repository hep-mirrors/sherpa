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
  msg_Out()<<METHOD<<"("<<m_beam<<"): p = "<<InMomentum()<<", "
	   <<"but x = "<<p_beam->X()<<"\n"
	   <<(*p_beamblob)<<"\n";
  return true;
}

bool No_Remnant::TestExtract(const Flavour &flav,const Vec4D &mom) {
  if (flav!=p_beam->Bunch()) return false;
  for (size_t i=0;i<4;i++) {
    double diff = ((mom[i]-p_beam->OutMomentum()[i])/
		   (mom[i]+p_beam->OutMomentum()[i])); 
    if (diff>1.e-6) {
      msg_Error()<<"Error in "<<METHOD<<": difference in four-momenta = "
		 <<diff<<" = "<<mom<<" - "<<p_beam->OutMomentum()<<".\n";
      return false;
    }
  }
  return true;
}
