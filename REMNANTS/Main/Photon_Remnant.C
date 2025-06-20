#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "REMNANTS/Main/Photon_Remnant.H"
#include "REMNANTS/Tools/Colour_Generator.H"
#include <algorithm>

using namespace REMNANTS;
using namespace ATOOLS;

Photon_Remnant::
Photon_Remnant(PDF::PDF_Base *pdf, const size_t & beam, const size_t & tag) :
  Remnant_Base(pdf->Bunch(), beam, tag),
  p_pdf(pdf), p_partons(&(pdf->Partons())),
  m_LambdaQCD(0.25), m_beta_quark(-1.), m_beta_gluon(-1.2),
  m_valence(false), p_spectator(nullptr),
  p_recoiler(nullptr) {
  p_ff = new Form_Factor(pdf->Bunch());
}

bool Photon_Remnant::FillBlob(Colour_Generator* colours, ParticleMomMap* ktmap, const bool& copy)
{
  if (m_extracted.empty()) {
    msg_Error() << METHOD
              << ": No remnants have been extracted, please check. \n";
    return false;
  }
  // In the photon, there needs to be at least one quark-antiquark pair,
  // this is tracked with the m_valence flag. Of these two, the antiquark will
  // be used as the recoiler later-on.
  MakeRemnants(colours);
  msg_Debugging() << METHOD << ": Filling blob with remnants, extracted = "
                  << m_extracted << ", \n and spectators = " << m_spectators
                  << "\n";
  FindRecoiler();
  // Possibly adjust final pending colours with extra gluons -
  // in prinicple one may have to check that they are not singlets ....
  CompensateColours(colours);
  // Assume all remnant bases already produced a beam blob = p_beamblob
  if (!MakeLongitudinalMomenta(ktmap, copy)) {
    msg_Debugging() << METHOD << ": Cannot put all particles on mass-shell, returning false.\n";
    return false;
  }
  if (!p_beamblob->CheckColour(true)) {
    msg_Error()<<"   * Error in "<<METHOD<<" (illegal colour) for \n"
	       <<(*p_beamblob)<<"\n";
    colours->Output();
    return false;
  }
  return true;
}

void Photon_Remnant::Reset(const bool & resc,const bool &DIS) {
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
  m_valence   = false;
  p_recoiler  = nullptr;
}

void Photon_Remnant::Output() const {
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
  // Still in range?
  double x = mom[0] / m_residualE;
  if (x < p_pdf->XMin() || x > p_pdf->XMax()) {
    msg_Tracking() << METHOD << ": out of limits, x = " << x << ".\n";
    return false;
  }
  return true;
}

bool Photon_Remnant::MakeLongitudinalMomenta(ParticleMomMap *ktmap,
                                             const bool &copy) {
  // Calculate the total momentum that so far has been extracted through
  // the shower initiators and use it to determine the still available
  // momentum; the latter will be successively reduced until the
  // rest is taken by the quark.
  Vec4D availMom = p_beam->OutMomentum(m_tag);
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
  msg_Debugging() << METHOD << ": "
		  <<"Longitudinal momentum left for remnants = " << availMom
                  << "\n";
  /* The momentum that remains needs to be distributed over the remnants.
   * Each parton should have an energy greater than its (hadron) mass, additionally energy must be conserved.
   * This requires book-keeping of two variables:
   * 1. The masses that still need too be generated
   * 2. The energy that has been used so far.
   * The kinetic energy is distributed by a rough fit to the photon PDFs in the
   * low-x and low-Q^2 region in the function SelectZ
   */
  double remnant_masses = 0.;
  for (Particle  const * pit : m_spectators) {
    remnant_masses += Max(pit->Flav().HadMass(), m_LambdaQCD);
  }
  if (remnant_masses > availMom[0]) {
    msg_Debugging() << METHOD
                    << ": Warning, HadMasses of remnants = " << remnant_masses
                    << " vs. residual energy = " << availMom[0] << "\n";
    return false;
  }
  for (auto part : m_spectators) {
    if (part == m_spectators.back()) {
      part->SetMomentum(availMom);
    } else {
      part->SetMomentum(SelectZ(part->Flav(), availMom[0], remnant_masses) *
			availMom);
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
  return true;
}

double Photon_Remnant::SelectZ(const Flavour &flav, double restmom,
                               double remnant_masses) const {
  // Give a random number to distribute longitudinal momenta, but this number
  // must respect the mass constraints
  double zmin = Max(flav.HadMass(), m_LambdaQCD) / restmom;
  double zmax = zmin + (restmom - remnant_masses) / restmom;
  // Taken from Hadron_Remnant, adapted the exponents for photon PDFs
  if (zmax < zmin) {
    msg_Debugging() << METHOD << ": Error, zmin, zmax = "
		    << zmin <<", "<<zmax << "\n";
    return 0;
  }
  double z;
  double beta = flav.IsGluon() ? m_beta_gluon : m_beta_quark;
  double invb = 1. / (beta + 1);
  if (beta!=-1) {
    double rand = ran->Get();
    z = pow(rand*pow(zmax,beta+1.)+(1.-rand)*pow(zmin,beta+1.),invb);
  }
  else
    z = zmin * pow(zmax/zmin,ran->Get());
  return z;
}

void Photon_Remnant::MakeSpectator(Particle* parton, Colour_Generator* colours)
{
  /* The remnant is constructed from the extracted shower initiators.
   * If it is a quark, we only have to generate the corresponding antiparticle.
   * If it is a gluon, we do nothing for the moment, but will later make sure,
   * that there is at least one quark-antiquark pair in the remnants.
   * */
  p_spectator = nullptr;
  Flavour flavour = parton->Flav();
  if (!flavour.IsQuark()) return;
  p_spectator = MakeParticle(flavour.Bar());
  int i = (p_spectator->Flav().IsAnti()?2:1);
  p_spectator->SetFlow(i, -1);
  p_spectator->SetPosition(parton->XProd());
  colours->AddColour(m_beam, i - 1, p_spectator);
  m_spectators.push_front(p_spectator);
  if (!m_valence)
    m_valence = true;
}

void Photon_Remnant::MakeRemnants(Colour_Generator* colours)
{
  if (m_valence)
    return;
  Particle * part;
  Flavour quark;
  double rand = ran->Get();
  if (rand < 4. / 6.)
    quark = kf_u;
  else if (rand < 5. / 6.)
    quark = kf_d;
  else
    quark = kf_s;
  int factor = 1;
  for (int i = 1; i < 3; i++) {
    part = MakeParticle(factor * quark);
    part->SetFlow(i, colours->NextColour(m_beam, i - 1));
    part->SetPosition(m_position+(*p_ff)());
    m_spectators.push_front(part);
    factor *= -1;
  }
  m_valence = true;
}

void Photon_Remnant::FindRecoiler() {
  for (auto part : m_spectators) {
    if (part->Flav().IsQuark())
      p_recoiler = part;
  }
}
