#include "REMNANTS/Main/Remnant_Base.H"
#include "REMNANTS/Tools/Colour_Generator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Momentum_Shifter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "BEAM/Main/Beam_Base.H"
#include <algorithm>


using namespace REMNANTS;
using namespace ATOOLS;

std::ostream &REMNANTS::operator<<(std::ostream &ostr,const rtp::code code)
{
  switch (code) {
  case rtp::none:   return ostr<<"None";
  case rtp::intact: return ostr<<"Intact";
  case rtp::hadron: return ostr<<"Hadron";
  case rtp::photon: return ostr<<"Photon";
  case rtp::lepton: return ostr<<"Lepton";
  }
  return ostr;
}

Remnant_Base::Remnant_Base(const rtp::code type,const unsigned int beam):
  m_type(type), m_beam(beam),
  p_beam(NULL), p_partner(NULL), p_beamblob(NULL), p_colours(NULL),
  m_rescale(true), m_scale2(-1.)
{}

Remnant_Base::~Remnant_Base() {}

bool Remnant_Base::Extract(ATOOLS::Particle *parton) 
{
  if (TestExtract(parton->Flav(),parton->Momentum())) {
    if (std::find(m_extracted.begin(),m_extracted.end(),parton)==m_extracted.end()) {
      m_extracted.push_back(parton);
      // right now this works for Hadron_Remnant only, all other remnants do not have
      // the notion of a spectator
      MakeSpectator(parton);
      for (size_t index=0;index<2;index++) p_colours->AddColour(m_beam,index,parton);
    }
    return true;
  }
  msg_Error()<<METHOD<<": Cannot extract particle:\n"<<(*parton)<<"\n  from:\n"
	     <<p_beam->Bunch()<<"\n";
  return false;
}

bool Remnant_Base::TestExtract(ATOOLS::Particle *parton) 
{ 
  if (parton==NULL) {
    msg_Error()<<"Error in "<<METHOD<<"():\n"
	       <<"   Called with NULL pointer.\n";
    return false;
  }
  return TestExtract(parton->Flav(),parton->Momentum());
}

Blob * Remnant_Base::MakeBlob() {
  p_beamblob = new Blob();
  p_beamblob->SetType(btp::Beam);
  p_beamblob->SetId();
  p_beamblob->SetBeam(m_beam);
  p_beamblob->SetStatus(blob_status::needs_beams |
			blob_status::needs_softUE);
  Particle * part = new Particle(-1,p_beam->Bunch(),p_beam->OutMomentum());
  part->SetNumber(0);
  part->SetBeam(m_beam);
  part->SetStatus(part_status::decayed);
  part->SetFinalMass();
  p_beamblob->AddToInParticles(part);
  return p_beamblob;
}

void Remnant_Base::Reset()
{
  m_extracted.clear();
  m_pbeam    = p_beam->InMomentum();
  p_beamblob = NULL;
}
