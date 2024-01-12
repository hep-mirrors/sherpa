#include "REMNANTS/Tools/Beam_Decorrelator.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace REMNANTS;
using namespace ATOOLS;
using namespace std;


Beam_Decorrelator::Beam_Decorrelator() : m_on(false) {}

Beam_Decorrelator::~Beam_Decorrelator() {}

void Beam_Decorrelator::
Initialize(Remnant_Handler * const rhandler) {
  p_rhandler = rhandler;
  if (p_rhandler->Type()==strat::DIS1 || p_rhandler->Type()==strat::DIS2) {
    p_kperpGenerator = p_rhandler->GetKPerp();
    auto s = Settings::GetMainSettings()["REMNANTS"];
    m_expo    = s["SOFT_X_EXPONENT"].SetDefault(-2.0).Get<double>();
    m_xiP     = m_expo+1;
    m_invxiP  = 1./m_xiP;
    m_maxeta  = dabs(s["SOFT_ETA_RANGE"].SetDefault(7.5).Get<double>());
    m_mass2   = sqr(s["SOFT_MASS"].SetDefault(5.0).Get<double>());
    m_deltaM  = s["DELTA_MASS"].SetDefault(1.5).Get<double>();
    m_on      = (m_maxeta>0. && m_mass2>0.);
  }
}

bool Beam_Decorrelator::operator()(Blob * softblob) {
  if (!m_on) return true;
  p_softblob = softblob;
  // Check for pairs of partons; if they are colour-correlated and involve either
  // * one beam parton and one "proper" shower parton with |eta| < maxeta, or
  // * two beam partons from different beams
  // and if their mass is larger than the minimal mass m_mass2, they will emit a soft
  // gluon to be decorrelated in colour.
  for (size_t i=0;i<p_softblob->NOutP()-1;i++) {
    for (size_t j=i+1;j<p_softblob->NOutP();j++) {
      if (MustEmit(p_softblob->OutParticle(i),p_softblob->OutParticle(j))) SoftEmission();
    }
  }
  // Add the produced soft gluons to the blob.
  bool out=!m_softgluons.empty();
  while (!m_softgluons.empty()) {
    p_softblob->AddToOutParticles(m_softgluons.back());
    m_softgluons.pop_back();
  }
  //if (out) msg_Out()<<(*p_softblob)<<"\n";
  return true;
}

bool Beam_Decorrelator::MustEmit(Particle * pi, Particle * pj) {
  //return false;
  // Checks if the partons must emit a soft gluon, for conditions see above.
  // Ignore parton pairs from shower or from the same beam breakup
  if (pi->Info()=='I' || pj->Info()=='I' ||
      (pi->Beam()<=0 && pj->Beam()<=0) ||
      (pi->Beam()==pj->Beam() && pi->Beam()!=0)) return false;
  // Ignore parton pairs that are not colour-correlated
  if (!((pi->GetFlow(1)==pj->GetFlow(2) && pi->GetFlow(1)!=0) ||
	(pi->GetFlow(2)==pj->GetFlow(1) && pi->GetFlow(2)!=0))) return false;
  if (pi->Beam()>0 &&
      ((pj->Beam()>0 && pi->Momentum()[0]>pj->Momentum()[0]) || pj->Beam()<=0)) {
    p_beam = pi; p_spect = pj;
  }
  else {p_beam = pj; p_spect = pi; }
  m_pbeam  = p_beam->Momentum();
  m_pspect = p_spect->Momentum();
  m_Q2     = (m_pspect+m_pbeam).Abs2();
  // Spectator parton must be from final state and inside eta-range or from other beam,
  // and invariant mass of pair must be above threshold
  return ((p_spect->Beam()<0 || dabs(m_pspect.Eta())<m_maxeta) &&
	  m_Q2 > m_mass2);
}

bool Beam_Decorrelator::SoftEmission() {
  InitSoftEmission();
  return DefineKinematics();
}

void Beam_Decorrelator::InitSoftEmission() {
  // Construct invariants and c.m. system of emitter-spectator pair
  m_Q       = sqrt(m_Q2);
  m_mbeam   = p_beam->Flav().HadMass();  m_mbeam2  = sqr(m_mbeam);
  m_mspect  = p_spect->Flav().HadMass(); m_mspect2 = sqr(m_mspect);
  m_boost   = Poincare(m_pbeam+m_pspect);
  m_boost.Boost(m_pbeam);
  m_boost.Boost(m_pspect);
  m_rotat   = Poincare(m_pspect, m_Q*s_AxisM);
  m_rotat.Rotate(m_pbeam);
  m_rotat.Rotate(m_pspect);
}

bool Beam_Decorrelator::DefineKinematics() {
  // 1000 trials to produce a kinematics that works, with x from a simple monomial x^(m_xiP-1)
  // and the transverse momentum vector from the Primordial_KPerp
  m_minMbeam    = m_mbeam+m_deltaM;
  m_minMbeam2   = sqr(m_minMbeam);
  m_minMspect   = m_mspect+m_deltaM;
  m_minMspect2  = sqr(m_minMspect);
  double eps    = m_minMbeam2/m_Q2;
  double poweps = pow(eps,m_xiP);
  int    trials = 1000;
  do {
    m_x     = pow(ran->Get()*(1-poweps)+poweps,m_invxiP);
    m_ktvec = p_kperpGenerator->KT(p_beam);
  } while (!MakeKinematics() && (trials--)>0);
  if (trials<=0)
    msg_Tracking()<<METHOD<<": couldn't construct Primordial_KPerp kinematics.\n";
  return (trials>0);
}

bool Beam_Decorrelator::MakeKinematics() {
  // Constructing a kinematics in the c.m. frame
  m_kt2    = dabs(m_ktvec.Abs2());
  double x = m_x;
  double y = m_kt2/(m_Q2*x);
  double alpha, beta;
  if (m_mspect2<1.e-12) {
    beta  = 1.-y-(m_mbeam2+m_kt2)/((1.-x)*m_Q2);
    alpha = 1.;
  }
  else {
    double h1   = 1.-y+(m_mbeam2-m_mspect2+m_kt2)/((1.-x)*m_Q2);
    double h2   = m_mspect2/m_Q2*(1.-y);
    double disc = h1*h1-4.*h2;
    if (disc<0.) return false;
    beta  = (h1+sqrt(disc))/2.;
    alpha = 1.-m_mspect2/(m_Q2*beta);
  }
  // Simple kinematics checks on parameters
  if (y>x || y<0 || beta>(1.-y) || alpha<0. || alpha>1.) {
    return false;
  }
  double E     = m_Q/2.;
  Vec4D pi     = E * ((   alpha-x) * s_AxisP + (1.-beta-y) * s_AxisM) + m_ktvec;
  Vec4D pj     = E * (          x  * s_AxisP +          y  * s_AxisM) - m_ktvec;
  Vec4D pk     = E * ((1.-alpha)   * s_AxisP +       beta  * s_AxisM);
  // Simple kinematics checks on momenta: right mass shells, enogh mass in the beam-gluon pair
  // positive energies for all particles
  if (!IsZero(dabs(pi.Abs2()) - m_mbeam2) || !IsZero(dabs(pj.Abs2()) - 0.) ||
      !IsZero(dabs(pk.Abs2()) - m_mspect2) ||
      !IsZero(dabs((pi+pj+pk).Abs2()) - m_Q2) ||
      !(pi[0]>0.) || !(pj[0]>0.) || !(pk[0]>0.) ||
      (pi+pj).Abs2()<m_minMbeam2 || (pj+pk).Abs2()<m_minMspect2) {
    return false;
  }
  m_rotat.RotateBack(pi);
  m_rotat.RotateBack(pj);
  m_rotat.RotateBack(pk);
  m_boost.BoostBack(pi);
  m_boost.BoostBack(pj);
  m_boost.BoostBack(pk);
  // Updating the momenta of beam parton (the emitter) and spectator, adding a soft gluon
  p_beam->SetMomentum(pi);
  p_spect->SetMomentum(pk);
  Particle * gluon = new Particle(-1,Flavour(kf_gluon),pj,'B');
  m_softgluons.push_back(gluon);
  return true;
}

void Beam_Decorrelator::Reset() {
  m_softgluons.clear();
}


