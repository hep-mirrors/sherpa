#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Momentum_Shifter.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include "REMNANTS/Tools/Colour_Generator.H"
#include <algorithm>

using namespace REMNANTS;
using namespace ATOOLS;

std::ostream &REMNANTS::operator<<(std::ostream &ostr, const rtp::code code) {
  switch (code) {
  case rtp::none:
    return ostr << "None";
  case rtp::intact:
    return ostr << "Intact";
  case rtp::hadron:
    return ostr << "Hadron";
  case rtp::photon:
    return ostr << "Photon";
  case rtp::lepton:
    return ostr << "Lepton";
  }
  return ostr;
}

Remnant_Base::
Remnant_Base(const ATOOLS::Flavour & flav, const size_t & beam, const size_t & tag) :
  m_beamflav(flav), m_type(FixType(m_beamflav)),
  m_beam(beam), m_tag(tag), p_beam(nullptr), p_ff(nullptr),
  p_beamblob(nullptr), m_position(Vec4D(0.,0.,0.,0.)),
  p_colours(nullptr), m_scale2(-1.) {}

Remnant_Base::~Remnant_Base() {
  if (p_ff=nullptr) { delete p_ff; p_ff = NULL; }
}

rtp::code Remnant_Base::FixType(ATOOLS::Flavour & flav) {
  if (flav.IsLepton()) return rtp::lepton;
  if (flav.IsHadron()) return rtp::hadron;
  if (flav.IsPhoton()) return rtp::photon;
  return rtp::none;
}

bool Remnant_Base::Extract(ATOOLS::Particle *parton) {
  ////////////////////////////////////////////////////////////////////////////
  // Extracting a parton from a remnant (usually stemming from a shower blob) 
  // and, if necessary, create a spectator to compensate flavour.
  ////////////////////////////////////////////////////////////////////////////
  if (TestExtract(parton->Flav(), parton->Momentum())) {
    if (std::find(m_extracted.begin(), m_extracted.end(), parton)==m_extracted.end()) {
      //msg_Out()<<"   * "<<METHOD<<"(beam = "<<m_beam<<"): add "<<parton->Number()<<"\n";
      m_extracted.push_back(parton);
      MakeSpectator(parton);
      for (size_t index = 0; index < 2; index++) p_colours->AddColour(m_beam, index, parton);
      m_residualE -= parton->E();
    }
    return true;
  }
  msg_Error() <<"==> "<<METHOD << ": Cannot extract particle:\n"
              << (*parton) << "\n  from: " << p_beam->Bunch()
              << " with momentum " << p_beam->OutMomentum() << "\n";
  return false;
}

bool Remnant_Base::TestExtract(ATOOLS::Particle *parton) {
  if (parton == nullptr) {
    msg_Error() << "Error in " << METHOD << "():\n"
                << "   Called with NULL pointer.\n";
    return false;
  }
  ////////////////////////////////////////////////////////////////////////////
  // TODO: In Multiple_Interactions.C, TestExtract is called for the hard-interacting 
  // parton, after it has been extracted here. As it has been already been extracted,
  // it will fail now if E_parton > 0.5 * beam energy - this must be checked.
  ////////////////////////////////////////////////////////////////////////////
  if (std::find(m_extracted.begin(), m_extracted.end(), parton) != m_extracted.end())
    return true;
  return TestExtract(parton->Flav(), parton->Momentum());
}

Blob *Remnant_Base::MakeBlob() {
  p_beamblob = new Blob();
  p_beamblob->SetType(btp::Beam);
  p_beamblob->SetId();
  p_beamblob->SetBeam(m_beam);
  p_beamblob->SetStatus(blob_status::needs_beams | blob_status::needs_softUE);
  p_beamblob->SetPosition(m_position);
  Particle *part = new Particle(-1, p_beam->Bunch(m_tag), p_beam->OutMomentum(m_tag));
  part->SetNumber(0);
  part->SetBeam(m_beam);
  part->SetStatus(part_status::decayed);
  part->SetFinalMass();
  p_beamblob->AddToInParticles(part);
  //msg_Out()<<"   * "<<METHOD<<"(beam = "<<m_beam<<"):\n"<<(*p_beamblob)<<"\n";
  return p_beamblob;
}

Vec4D Remnant_Base::IncomingMomentum() { return p_beam->OutMomentum(m_tag); }

void Remnant_Base::Reset(const bool & resc,const bool &DIS) {
  //msg_Out()<<"   * "<<METHOD<<"(beam = "<<m_beam<<")\n";
  m_extracted.clear();
  p_beamblob = nullptr;
}
