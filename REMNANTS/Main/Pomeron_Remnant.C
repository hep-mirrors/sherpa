#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "REMNANTS/Main/Pomeron_Remnant.H"
#include <algorithm>

using namespace REMNANTS;
using namespace ATOOLS;

Pomeron_Remnant::Pomeron_Remnant(PDF::PDF_Base *pdf, const size_t& beam, const size_t& tag)
    : Remnant_Base(Flavour(kf_pomeron), beam, tag), p_pdf(pdf),
      p_partons(&(pdf->Partons())), p_recoiler(nullptr),
      m_LambdaQCD(0.25) {
  p_ff = new Form_Factor(pdf->Bunch());
}

Particle *Pomeron_Remnant::MakeParticle(const ATOOLS::Flavour &flav) {
  Particle *part = new Particle(-1, flav, Vec4D(0., 0., 0., 0.), 'B');
  part->SetNumber();
  part->SetBeam(m_beam);
  return part;
}

bool Pomeron_Remnant::FillBlob(ParticleMomMap *ktmap, const bool &copy) {
  if (m_extracted.empty()) {
    msg_Error() << METHOD
                << ": No remnants have been extracted, please check. \n";
    return false;
  }
  if (m_extracted.size()>1) {
    msg_Error()<<METHOD<<": Too many remnants have been extracted, please check.\n";
    return false;
  }
  MakeRemnants();
  msg_Debugging() << METHOD << ": Filling blob with remnants, extracted = "
                  << m_extracted << ", \n and spectators = " << m_spectators
                  << "\n";
  // Assume all remnant bases already produced a beam blob = p_beamblob
  MakeLongitudinalMomenta(ktmap, copy);
  return true;
}

void Pomeron_Remnant::Reset(const bool& resc,const bool &DIS) {
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
  m_residualE = p_beam->OutMomentum(m_tag)[0];
  p_recoiler = nullptr;
}

void Pomeron_Remnant::Output() const {
  msg_Out() << METHOD << "(" << m_beam << ", " << m_beamflav << ").\n"
            << "   Partons are { "<< Flavour(kf_gluon) << " }.\n";
}

bool Pomeron_Remnant::TestExtract(const Flavour &flav, const Vec4D &mom) {
  // Is flavour element of flavours allowed by PDF?
  if (p_partons->find(flav) == p_partons->end()) {
    msg_Error() << METHOD << ": flavour " << flav << " not found.\n";
    return false;
  }
  if (mom[0] < flav.HadMass()) {
    msg_Debugging() << METHOD << ": parton too soft, mass = " << flav.HadMass()
                    << " and energy = " << mom[0] << "\n";
    return false;
  }
  double required_energy =
      flav.IsGluon() ? m_LambdaQCD + mom[0]
                     : Max(flav.HadMass(),m_LambdaQCD) + mom[0] + m_LambdaQCD;
  if (m_residualE < required_energy) {
    msg_Debugging() << METHOD << ": not enough energy to accomodate particle mass. \n";
    return false;
  }
  double x = mom[0] / m_residualE;
  // Still in range?
  if (x < p_pdf->XMin() || x > p_pdf->XMax()) {
    msg_Error() << METHOD << ": out of limits, x = " << x << ".\n";
    return false;
  }
  return true;
}

void Pomeron_Remnant::MakeLongitudinalMomenta(ParticleMomMap *ktmap,
                                             const bool &copy) {
  // Calculate the total momentum that so far has been extracted through
  // the shower initiators and use it to determine the still available
  // momentum; the latter will be successively reduced until the
  // rest is taken by the quark.
  Vec4D availMom;
  auto part_extr = m_extracted.front();
  availMom = p_beam->OutMomentum(m_tag) - part_extr->Momentum();
  if (copy) {
    Particle *pcopy = new Particle(*part_extr);
    pcopy->SetNumber();
    pcopy->SetBeam(m_beam);
    p_beamblob->AddToOutParticles(pcopy);
  } else
    p_beamblob->AddToOutParticles(part_extr);
  (*ktmap)[part_extr] = Vec4D();
  msg_Debugging() << METHOD << ": Longitudinal momentum left for remnants = " << availMom
                  << "\n";
  double remnant_masses = 0.;
  for (Particle  const * pit : m_spectators) {
    remnant_masses += Max(pit->Flav().HadMass(), m_LambdaQCD);
  }
  if (remnant_masses > m_residualE)
    msg_Error() << METHOD << ": Warning, HadMasses of remnants = "
                << remnant_masses << " vs. residual energy = " << m_residualE << "\n";
  for (auto part : m_spectators) {
    if (availMom[0] < 0)
      msg_Error() << METHOD << ": Negative Energy in Remnants! \n";
    if (part == m_spectators.back()) {
      part->SetMomentum(availMom);
    } else {
      part->SetMomentum(SelectZ(part->Flav(), availMom[0], remnant_masses) * availMom);
      availMom -= part->Momentum();
      remnant_masses -= Max(part->Flav().HadMass(), m_LambdaQCD);
    }
  msg_Debugging() << METHOD << ": set momentum for "<<part->Flav()<<" to "
                  << part->Momentum() << "\n";
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

void Pomeron_Remnant::MakeRemnants() {
  Particle * part;
  if (m_extracted.front()->Flav()==Flavour(kf_gluon)){
    part = p_recoiler = MakeParticle(Flavour(kf_gluon));
    for (int i = 1; i < 3; ++i) {
      part->SetFlow(i, m_extracted.front()->GetFlow(3-i));
    }
    part->SetPosition(m_position+(*p_ff)());
    m_spectators.push_front(part);
  } else {
    part = p_recoiler = MakeParticle(m_extracted.front()->Flav().Bar());
    part->SetFlow(part->Flav().IsAnti() ? 2 : 1, Flow::Counter());
    Particle * g = MakeParticle(Flavour(kf_gluon));
    for (int i = 1; i < 3; ++i) {
      int c = part->GetFlow(i) != 0 ? part->GetFlow(i) : m_extracted.front()->GetFlow(i);
      g->SetFlow(3-i, c);
    }
    part->SetPosition(m_position+(*p_ff)());
    g->SetPosition(m_position+(*p_ff)());
    m_spectators.push_front(part);
    m_spectators.push_front(g);
  }
}
double Pomeron_Remnant::SelectZ(const Flavour &flav, double restmom,
                                double remnant_masses) const {
  // Give a random number to distribute longitudinal momenta, but this number
  // must respect the mass constraints
  double zmin = Max(flav.HadMass(), m_LambdaQCD) / restmom;
  double zmax = zmin + (restmom - remnant_masses) / restmom;
  // Taken from Hadron_Remnant, adapted the exponents for photon PDFs
  if (zmax < zmin) {
    msg_Debugging() << METHOD << ": Error, zmin, zmax = " << zmin <<", "<<zmax << "\n";
    return 0;
  }
  double z;
  double beta = -1.5; // Assumed to be the same as in hadron remnants
  double invb = 1. / (beta + 1);
  double rand = ran->Get();
  z = pow(rand*pow(zmax,beta+1.)+(1.-rand)*pow(zmin,beta+1.),invb);
  return z;
}
