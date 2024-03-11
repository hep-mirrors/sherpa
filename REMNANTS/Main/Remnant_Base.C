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
  case rtp::pomeron:
    return ostr << "Pomeron";
  case rtp::reggeon: return ostr << "Reggeon";
  }
  return ostr;
}

Remnant_Base::Remnant_Base(const ATOOLS::Flavour& flav, const size_t& beam, const size_t& tag)
    : m_beamflav(flav), m_type(FixType(m_beamflav)), m_beam(beam), m_tag(tag), p_beam(nullptr),
      p_ff(nullptr), p_beamblob(nullptr), m_position(Vec4D(0., 0., 0., 0.)), m_residualE(0.),
      m_scale2(-1.)
{}

Remnant_Base::~Remnant_Base() {
  if (p_ff!=nullptr) { delete p_ff; p_ff = nullptr; }
}

rtp::code Remnant_Base::FixType(ATOOLS::Flavour & flav) {
  if (flav.Kfcode() == kf_pomeron) return rtp::pomeron;
  if (flav.Kfcode() == kf_reggeon) return rtp::reggeon;
  if (flav.IsLepton()) return rtp::lepton;
  if (flav.IsHadron()) return rtp::hadron;
  if (flav.IsPhoton()) return rtp::photon;
  return rtp::none;
}

Particle* Remnant_Base::MakeParticle(const Flavour& flav)
{
  Particle* part = new Particle(-1, flav, Vec4D(0., 0., 0., 0.), 'B');
  part->SetNumber();
  part->SetBeam(m_beam);
  part->SetPosition(m_position + (*p_ff)());
  return part;
}

void Remnant_Base::CompensateColours(Colour_Generator* colours)
{
  while (!colours->Colours(m_beam, 0).empty() && !colours->Colours(m_beam, 1).empty() &&
         colours->Colours(m_beam, 0) != colours->Colours(m_beam, 1)) {
    Particle* gluon = MakeParticle(Flavour(kf_gluon));
    for (size_t i = 0; i < 2; i++) gluon->SetFlow(i + 1, colours->NextColour(m_beam, i));
    gluon->SetPosition(m_position + (*p_ff)());
    m_spectators.push_back(gluon);
  }
}

bool Remnant_Base::Extract(ATOOLS::Particle* parton, Colour_Generator* colours)
{
  // Extracting a parton from a remnant (usually stemming from a shower blob)
  // and, if necessary, create a spectator to compensate flavour.
  if (TestExtract(parton->Flav(), parton->Momentum())) {
    if (std::find(m_extracted.begin(), m_extracted.end(), parton)==m_extracted.end()) {
      m_extracted.push_back(parton);
      // Spectators compensate for flavour, i.e. they are only created for quarks.
      MakeSpectator(parton, colours);
      for (size_t index = 0; index < 2; index++) colours->AddColour(m_beam, index, parton);
      m_residualE -= parton->E();
    }
    return true;
  }
  msg_Error() << METHOD << ": Cannot extract particle:\n"
              << (*parton) << "\n  from: " << p_beam->Bunch()
              << " with momentum " << p_beam->OutMomentum()
              << ", difference = " << parton->Momentum()-p_beam->OutMomentum()
              << "\n";
  return false;
}

bool Remnant_Base::TestExtract(ATOOLS::Particle *parton) {
  if (parton == nullptr) {
    msg_Error() << "Error in " << METHOD << "():\n"
                << "   Called with NULL pointer.\n";
    return false;
  }
  // TODO: In Multiple_Interactions.C, TestExtract is called for the hard-interacting
  // parton, after it has been extracted here. As it has been already been extracted,
  // it will fail now if E_parton > 0.5 * beam energy - this must be checked.
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
  return p_beamblob;
}

Vec4D Remnant_Base::IncomingMomentum() { return p_beam->OutMomentum(m_tag); }

void Remnant_Base::Reset(const bool & resc,const bool &DIS) {
  m_extracted.clear();
  p_beamblob = nullptr;
}
