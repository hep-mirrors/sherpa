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
  if (m_extracted.empty() || m_spectators.empty()) {
    msg_Error() << METHOD
              << ": No remnants have been extracted, please check. \n";
    return false;
  }
  m_energy = p_beam->OutMomentum()[0];
  // Assume all remnant bases already produced a beam blob = p_beamblob
  MakeLongitudinalMomenta(ktmap, copy);
  if (!p_beamblob->CheckColour(true)) {
    msg_Error() << "Error in " << METHOD << " for \n" << (*p_beamblob) << "\n";
    p_colours->Output();
    return false;
  }
  return true;
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
  m_energy = p_beam->OutMomentum()[0];
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
  double x = mom[0] / m_energy;
  // Still enough energy?
  if (x > 1.) {
    msg_Debugging() << METHOD << ": too much momentum " << mom[0] << " "
                    << "> E = " << m_energy << " for beam " << m_beam
                    << "\n";
    return false;
  }
  // Still in range?
  if (x < p_pdf->XMin() || x > p_pdf->XMax()) {
    msg_Error() << METHOD << ": out of limits, x = " << x << ".\n";
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
    if (availMom[0] < 0)
      msg_Error() << METHOD << ": Negative Energy in Remnants! \n";
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
  while (std::isnan(random) || random >= 1. || random <= 0.) {
    msg_Error() << METHOD
                << ": Something went wrong in the momentum distribution for "
                   "the photon remnants. \n"
                << " Will retry. \n";
    random = ran->Get();
  }
  return random;
}

void Photon_Remnant::MakeSpectator(Particle *parton) {
  /* The remnant is constructed from the extracted shower initiators.
   * If it is a quark, we only have to generate the corresponding antiparticle.
   * If it is a gluon, we have to construct a quark-antiquark pair with
   * corresponding colours.
   * */
  Particle * part = NULL;
  Flavour flavour = parton->Flav();
  if (flavour.IsQuark()) {
    part = MakeParticle(flavour.Bar());
    int i = (part->Flav().IsAnti()?2:1);
    part->SetFlow(i, p_colours->NextColour(m_beam,i-1));
    //part->SetFlow((flavour.Bar().IsAnti()?2:1),-1);
    p_colours->AddColour(m_beam,(part->Flav().IsAnti()?1:0),part);
    m_spectators.push_front(part);
  } else if (flavour.IsGluon()) {
    // TODO: implement the gluon treatment
    // For now, implement simple treatment by choosing between the light
    // quarks according to their squared charge
    int factor = 1;
    Flavour quark;
    double rand = ran->Get();
    if (rand < 4. / 6.)
      quark = kf_u;
    else if (rand < 5. / 6.)
      quark = kf_d;
    else
      quark = kf_s;
    for (int i = 1; i < 3; i++) {
      part = MakeParticle(factor * quark);
      part->SetFlow(i, p_colours->NextColour(m_beam,i-1));
      p_colours->AddColour(m_beam,(part->Flav().IsAnti()?1:0),part);
      m_spectators.push_front(part);
      factor *= -1;
    }
  }
}
