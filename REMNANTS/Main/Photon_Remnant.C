#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "REMNANTS/Main/Photon_Remnant.H"
#include "REMNANTS/Tools/Colour_Generator.H"
#include <algorithm>

using namespace REMNANTS;
using namespace ATOOLS;

Photon_Remnant::Photon_Remnant(PDF::PDF_Base *pdf, const unsigned int beam)
    : Remnant_Base(rtp::photon, beam), p_pdf(pdf),
      p_partons(&(pdf->Partons())) {}

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
  // CompensateColours();
  // Assume all remnant bases already produced a beam blob = p_beamblob
  MakeLongitudinalMomenta(ktmap, copy);
  if (!p_beamblob->CheckColour(true)) {
    msg_Error() << "Error in " << METHOD << " for \n" << (*p_beamblob) << "\n";
    p_colours->Output();
    return false;
  }
  return true;
}

void Photon_Remnant::CompensateColours() {
  /* Should not be necessary at the moment, as colours are always compensated
   * within the remnant.
   * */
  while (!p_colours->Colours(m_beam, 0).empty() &&
         !p_colours->Colours(m_beam, 1).empty()) {
    Particle *gluon = MakeParticle(Flavour(kf_gluon));
    for (int i = 0; i < 2; i++)
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

bool Photon_Remnant::TestExtract(const Flavour &flav, const Vec4D &mom) {
  // Is flavour element of flavours allowed by PDF?
  if (p_partons->find(flav) == p_partons->end()) {
    msg_Error() << METHOD << ": flavour " << flav << " not found.\n";
    return false;
  }
  // Still enough energy?
  if (mom[0] > m_residualE) {
    msg_Error() << METHOD << ": too much momentum " << mom[0] << " "
                << "> E = " << m_residualE << " for beam " << m_beam << "\n";
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
  Vec4D availMom;
  availMom = p_beam->OutMomentum();
  for (auto pmit : m_extracted) {
    availMom -= pmit->Momentum();
    if (copy) {
      Particle *pcopy = new Particle(*pmit);
      pcopy->SetNumber();
      pcopy->SetBeam(m_beam);
      p_beamblob->AddToOutParticles(pcopy);
    } else
      p_beamblob->AddToOutParticles(pmit);
    (*ktmap)[pmit] = Vec4D();
  }
  for (auto pmit : m_spectators) {
    Particle *part = pmit;
    if (part == m_spectators.back())
      part->SetMomentum(availMom);
    else {
      part->SetMomentum(SelectZ(part->Flav()) * availMom);
      availMom -= part->Momentum();
    }
    if (copy) {
      Particle *pcopy = new Particle(*part);
      pcopy->SetNumber();
      pcopy->SetBeam(m_beam);
      p_beamblob->AddToOutParticles(pcopy);
    } else
      p_beamblob->AddToOutParticles(part);
    (*ktmap)[part] = Vec4D();
  }
}

double Photon_Remnant::SelectZ(const Flavour &flav) {
  double random = ran->Get();
  if (isnan(random) || random >= 1. || random <= 0.) {
    msg_Error() << METHOD
                << ": Something went wrong in the momentum distribution for "
                   "the photon remnants. \n"
                << " Will retry. \n";
    return SelectZ(flav);
  }
  return random;
}

void Photon_Remnant::MakeSpectator(Particle *parton) {
  /* As of now, we only use one interaction in the photon processes.
   * Treatment of Multiple Parton Interactions will need special considerations
   * and will be implemented later on.
   * */
}

bool Photon_Remnant::MakeRemnants() {
  /* The remnant is constructed from the extracted shower initiators.
   * If it is a quark, we only have to generate the corresponding antiparticle.
   * If it is a gluon, we have to construct a quark-antiquark pair with
   * corresponding colours.
   * */
  if (m_extracted.empty()) {
    msg_Out() << METHOD
              << ": No remnants have been extracted, please check. \n";
    return false;
  }
  if (m_extracted.size() > 1.) {
    msg_Error() << METHOD
                << ": Mulitple Interactions are not yet implemented in the "
                   "photon remnant!\n";
  }
  // TODO: the for-loop below is actually redundant at the moment, because
  // we're trying to implement for one interacting particle only
  for (auto pmit : m_extracted) {
    if (pmit->Flav().IsGluon()) {
      // TODO: implement the gluon treatment
      // For now, implement simple treatment by choosing between the light
      // quarks according to their squared charge
      int factor = 1;
      Flavour flav;
      double rand = ran->Get();
      if (rand < 4. / 6.)
        flav = kf_u;
      else if (rand < 5. / 6.)
        flav = kf_d;
      else
        flav = kf_s;
      for (int i = 1; i < 3; i++) {
        p_remnant = p_recoiler = MakeParticle(factor * flav);
        p_remnant->SetFlow(i, pmit->GetFlow(3 - i));
        m_spectators.push_front(p_recoiler);
        factor = -1;
      }
      return true;
    } else {
      p_remnant = p_recoiler = MakeParticle(pmit->Flav().Bar());
      int index = p_remnant->Flav().IsAnti() ? 1 : 0;
      p_remnant->SetFlow(index + 1, pmit->GetFlow(2 - index));
      m_spectators.push_front(p_recoiler);
      return true;
    }
  }
  msg_Error() << METHOD << ": Remnants could not be created. \n";
  return false;
}