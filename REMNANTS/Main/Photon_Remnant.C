#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "BEAM/Main/Beam_Base.H"
#include "REMNANTS/Main/Photon_Remnant.H"
#include "REMNANTS/Tools/Colour_Generator.H"
#include <algorithm>

using namespace REMNANTS;
using namespace ATOOLS;

Photon_Remnant::Photon_Remnant(PDF::PDF_Base *pdf, const unsigned int beam)
    : Remnant_Base(rtp::photon, beam), p_pdf(pdf), p_partons(&(pdf->Partons())),
      m_beamflav(pdf->Bunch()), p_remnant(nullptr), p_recoiler(nullptr),
      p_spectator(nullptr), m_alpha(0.), m_gamma(1.), m_beta(-1.5),
      m_invb(1. / (m_beta + 1)), m_LambdaQCD(0.25) {
  m_scale2 = Max(4.0, pdf->Q2Min());
}

bool Photon_Remnant::FillBlob(ATOOLS::Blob *beamblob,
                              ATOOLS::Particle_List *particlelist) {
  if (p_partner == nullptr) {
    THROW(critical_error, "Partner Remnant not set.");
  }
  for (auto &pmit : m_extracted) {
    beamblob->AddToOutParticles(pmit);
    if (particlelist != nullptr) {
      pmit->SetNumber(particlelist->size());
      particlelist->push_back(pmit);
    }
  }
  return true;
}

Particle *Photon_Remnant::MakeParticle(const Flavour &flav) {
  Particle *part = new Particle(-1, flav, Vec4D(0., 0., 0., 0.), 'B');
  part->SetNumber();
  part->SetBeam(m_beam);
  return part;
}

bool Photon_Remnant::FillBlob(ParticleMomMap *ktmap, const bool &copy) {
  m_residualE = p_beam->OutMomentum()[0];
  // Add remnants, diquark and quark, if necessary.
  if (!p_remnant)
    MakeRemnants();
  // Possibly adjust final pending colours with extra gluons - in prinicple one
  // may have to check that they are not singlets ....
  CompensateColours();
  // Assume all remnant bases already produced a beam blob = p_beamblob
  MakeLongitudinalMomenta(ktmap, copy);
  bool colourconserved = p_beamblob->CheckColour(true);
  if (!colourconserved) {
    msg_Error() << "Error in " << METHOD << " for \n" << (*p_beamblob) << "\n";
    p_colours->Output();
    return false;
  }
  return true;
}

void Photon_Remnant::CompensateColours() {
  while (!p_colours->Colours(m_beam, 0).empty() &&
         !p_colours->Colours(m_beam, 1).empty()) {
    Particle *gluon = MakeParticle(Flavour(kf_gluon));
    for (size_t i = 0; i < 2; i++)
      gluon->SetFlow(i + 1, p_colours->NextColour(m_beam, i));
    m_spectators.push_back(gluon);
  }
}

void Photon_Remnant::Reset(const bool &DIS) {
  Remnant_Base::Reset();
  while (!m_spectators.empty()) {
    Particle *part = m_spectators.front();
    if (part->ProductionBlob())
      part->ProductionBlob()->RemoveOutParticle(part);
    if (part->DecayBlob())
      part->DecayBlob()->RemoveInParticle(part);
    delete part;
    m_spectators.pop_front();
  }
  m_spectators.clear();
  m_residualE = p_beam->OutMomentum()[0];
  p_remnant = p_recoiler = nullptr;
}

void Photon_Remnant::Output() {
  msg_Out() << METHOD << "(" << m_beam << ", " << m_beamflav << ").\n"
            << "   Partons are { ";
  for (const auto &p_parton : *p_partons) {
    msg_Out() << " " << p_parton;
  }
  msg_Out() << "}.\n";
}

Flavour Photon_Remnant::RemnantFlavour(const Flavour &flav) {
  return Flavour(-1. * flav.Kfcode());
}

bool Photon_Remnant::TestExtract(const Flavour &flav, const Vec4D &mom) {
  // Is flavour element of flavours allowed by PDF?
  if (p_partons->find(flav) == p_partons->end()) {
    msg_Error() << METHOD << ": flavour " << flav << " not found.\n";
    return false;
  }
  // Still enough energy?
  if (mom[0] > m_residualE) {
    msg_Error() << METHOD << ": too much momentum " << mom[0] << " "
                << "> E = " << m_residualE << ".\n";
    return false;
  }
  // Still enough energy?  And in range?
  m_x = mom[0] / (m_rescale ? m_residualE : p_beam->OutMomentum()[0]);
  if (m_x < p_pdf->XMin() || m_x > p_pdf->XMax()) {
    msg_Error() << METHOD << ": out of limits, x = " << m_x << ".\n";
    return false;
  }
  return true;
}

void Photon_Remnant::MakeLongitudinalMomenta(ParticleMomMap *ktmap,
                                             const bool &copy) {
  // Calculate the total momentum that so far has been extracted through
  // the shower initiators and use it to determine the still available
  // momentum; the latter will be successively reduced until the
  // rest is taken by the quark.
  // TODO: check if spectators are handled correctly
  Vec4D availMom;
  availMom = p_beam->OutMomentum();
  for (Part_Iterator pmit = m_extracted.begin(); pmit != m_extracted.end();
       pmit++) {
    availMom -= (*pmit)->Momentum();
    p_beamblob->AddToOutParticles(*pmit);
  }
  for (Part_Iterator pmit = m_spectators.begin(); pmit != m_spectators.end();
       pmit++) {
    Particle *part = (*pmit);
    part->SetMomentum(availMom);
    p_beamblob->AddToOutParticles(part);
  }
}

void Photon_Remnant::MakeSpectator(Particle *parton) {
  // If a shower initiator is a sea-quark or antiquark, a corresponding
  // antiflavour has to be added to the spectators.
  p_spectator = nullptr;
  Flavour flav = parton->Flav();
  if (flav.IsQuark()) {
    p_spectator = MakeParticle(flav.Bar());
    p_spectator->SetFlow((flav.Bar().IsAnti() ? 2 : 1), -1);
    p_colours->AddColour(m_beam, (flav.Bar().IsAnti() ? 1 : 0), p_spectator);
    m_spectators.push_front(p_spectator);
    msg_Out() << METHOD << ": Spectator created for beam " << m_beam
              << " with colour "
              << p_spectator->GetFlow((flav.Bar().IsAnti() ? 2 : 1)) << "\n";
    msg_Out() << METHOD << ": p_colours = " << p_colours << "\n";
  }
  if (flav.IsGluon()) {
    // TODO: Get quark-antiquark pair with colours corresponding to the gluon
  }
}

bool Photon_Remnant::MakeRemnants() {
  // If no valence quark has been extracted to date, a quark-diquark
  // pair must be constructed.  the idea is to pick one of the three flavours
  // at random for the quark, add it to the spectators, then construct the
  // "conjugate" diquark and add it as well to the spectators
  /*  Flavour flav;
    size_t index;
      int random = int(ran->Get() * (p_partons->size()-3));
      std::_Rb_tree_const_iterator<Flavour> flit = p_partons->begin();
      for (size_t i = 0; i < random; i++)
        flit++;
      flav = (*flit);
      Particle *part = MakeParticle(flav);
      index = ((flav.IsQuark() && !flav.IsAnti()) ||
               (flav.IsDiQuark() && flav.IsAnti()))
                  ? 0
                  : 1;
      part->SetFlow(index + 1, p_colours->NextColour(m_beam, index));
      m_spectators.push_back(part);

      p_remnant = p_recoiler = MakeParticle(RemnantFlavour(flav));
      p_remnant->SetFlow(2 - index, p_colours->NextColour(m_beam, 1 - index));
      m_spectators.push_front(p_recoiler);
      return true;*/
  size_t index;
  if (!m_extracted.empty()) {
    // TODO: this below is actually redundant at the moment, because we're
    // trying to implement for one remnant particle only
    for (auto &pmit : m_extracted) {
      Flavour flav = pmit->Flav();
      // Particle *part = MakeParticle(flav);
      index = (flav.IsQuark() && !flav.IsAnti()) ? 0 : 1;
      // part->SetFlow(index + 1, p_colours->NextColour(m_beam, index));
      // m_spectators.push_back(part);
      p_remnant = p_recoiler = MakeParticle(flav);
      p_remnant->SetFlow(2 - index, p_colours->NextColour(m_beam, 1 - index));
      // p_colours->AddColour(m_beam, index, p_remnant);
      msg_Out() << METHOD << ": Remnant created for beam " << m_beam
                << " with colour " << p_remnant->GetFlow(index + 1) << "\n";
      msg_Out() << METHOD << ": p_colours = " << p_colours << "\n";
      // m_spectators.push_front(p_recoiler);
    }
    return true;
  }
  msg_Out() << METHOD << ": No remnants have been extracted, please check. \n";
  return false;
}