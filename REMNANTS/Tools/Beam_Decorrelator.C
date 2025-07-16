#include "REMNANTS/Tools/Beam_Decorrelator.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace REMNANTS;
using namespace ATOOLS;

Beam_Decorrelator::Beam_Decorrelator() : m_on(false) {}

void Beam_Decorrelator::
Initialize(Remnant_Handler * const rhandler) {
  // Registering defaults
  auto s = Settings::GetMainSettings()["BEAM_DECORRELATOR"];
  s["SOFT_X_EXPONENT"].SetDefault(-2.0);
  s["SOFT_ETA_RANGE"].SetDefault(10.);
  s["SOFT_MASS"].SetDefault(7.5);
  s["DELTA_MASS"].SetDefault(0.1);

  p_rhandler = rhandler;
  m_on   = s["BEAM_DECORRELATOR"].SetDefault(0).Get<bool>();
  if (m_on &&
      (p_rhandler->Type()==strat::DIS1 || p_rhandler->Type()==strat::DIS2)) {
    // TODO: do the same for hh collisions?
    p_kperpGenerator = p_rhandler->GetKPerp();

    m_expo   = s["SOFT_X_EXPONENT"].Get<double>();
    m_xiP    = m_expo + 1;
    m_invxiP = 1. / m_xiP;
    m_maxeta = dabs(s["SOFT_ETA_RANGE"].Get<double>());
    m_mass2  = sqr(s["SOFT_MASS"].Get<double>());
    m_deltaM = s["DELTA_MASS"].Get<double>();
    m_on     = (m_maxeta > 0. && m_mass2 > 0.);
  }
}

bool Beam_Decorrelator::operator()(Blob * softblob) {
  if (!m_on) return true;
  p_softblob = softblob;
  // Check for pairs of partons; if they are colour-correlated and involve
  // * one beam and one FS parton, with |eta| < maxeta (TODO: not sure here?)
  // * two beam partons from different beams
  // and if their mass is larger than the minimal mass m_mass2, they will emit a
  // soft gluon to be decorrelated in colour.
  for (int i = 0; i < p_softblob->NOutP() - 1; i++) {
    for (int j = i + 1; j < p_softblob->NOutP(); j++) {
      if (MustEmit(p_softblob->OutParticle(i), p_softblob->OutParticle(j)))
        SoftEmission();
    }
  }
  // Add the produced soft gluons to the blob.
  while (!m_softgluons.empty()) {
    p_softblob->AddToOutParticles(m_softgluons.back());
    m_softgluons.pop_back();
  }
  return true;
}

bool Beam_Decorrelator::MustEmit(Particle * pi, Particle * pj) {
  // Checks if the partons must emit a soft gluon, for conditions see above.
  //msg_Out()<<METHOD<<": "<<pi->Number()<<" & "<<pj->Number()<<"\n";
  if (pi->Beam() < 0 && pj->Beam() < 0 ||
      pi->Info()=='I' || pj->Info()=='I')
    return false;
  // Ignore parton pairs that are not colour-correlated
  if (!((pi->GetFlow(1) == pj->GetFlow(2) && pi->GetFlow(1) != 0) ||
        (pi->GetFlow(2) == pj->GetFlow(1) && pi->GetFlow(2) != 0)))
    return false;
  if (pi->Info()=='B' || pj->Info()=='F' || pi->Momentum()[0] > pj->Momentum()[0]) {
    p_beam  = pi;
    p_spect = pj;
  } else {
    p_beam  = pj;
    p_spect = pi;
  }
  m_pbeam  = p_beam->Momentum();
  m_pspect = p_spect->Momentum();
  m_Q2     = (m_pspect+m_pbeam).Abs2();
  // Spectator parton must be from final state and inside eta-range or
  // from other beam, and invariant mass of pair must be above threshold
  return ( m_Q2 > m_mass2 &&
	   (p_spect->Beam()>=0 ||
	    (p_spect->Beam()<0 && dabs(m_pspect.Eta())<m_maxeta)));
}

bool Beam_Decorrelator::SoftEmission() {
  InitSoftEmission();
  Vec4D pi(0.,0.,0.,0.),pj(0.,0.,0.,0.),pk(0.,0.,0.,0.);
  if (!DefineKinematics(pi,pj,pk)) return false;
  return FillSoftEmission(pi,pj,pk);
}

void Beam_Decorrelator::InitSoftEmission() {
  // Construct invariants and c.m. system of emitter-spectator pair
  m_Q       = sqrt(m_Q2);
  m_mbeam2  = m_pbeam.Abs2();  m_mbeam  = sqrt(m_mbeam2);
  m_mspect2 = m_pspect.Abs2(); m_mspect = sqrt(m_mspect2);
  m_boost   = Poincare(m_pbeam+m_pspect);
  m_boost.Boost(m_pbeam);
  m_boost.Boost(m_pspect);
  m_rotat   = Poincare(m_pspect, m_Q*s_AxisM);
  m_rotat.Rotate(m_pbeam);
  m_rotat.Rotate(m_pspect);
}

bool Beam_Decorrelator::DefineKinematics(Vec4D & pi,Vec4D & pj,Vec4D & pk) {
  // 1000 trials to produce a kinematics that works, with x from a simple
  // monomial x^(m_xiP-1) and the transverse momentum vector from the
  // Primordial_KPerp
  m_minMbeam    = m_mbeam+m_deltaM;
  m_minMbeam2   = sqr(m_minMbeam);
  m_minMspect   = m_mspect+m_deltaM;
  m_minMspect2  = sqr(m_minMspect);
  double eps    = m_minMbeam2/m_Q2;
  double poweps = pow(eps,m_xiP);
  int    trials = 10;
  do {
    m_x     = pow(ran->Get()*(1-poweps)+poweps,m_invxiP);
    m_ktvec = p_kperpGenerator->KT(p_beam,m_Q2*m_x);
    trials--;
  } while (!MakeKinematics(pi,pj,pk) && trials > 0);
  if (trials<=0)
    msg_Tracking()<<METHOD<<": couldn't construct kinematics.\n";
  return (trials>0);
}

bool Beam_Decorrelator::MakeKinematics(Vec4D & pi,Vec4D & pj,Vec4D & pk) {
  // Constructing a kinematics in the c.m. frame
  m_kt2    = dabs(m_ktvec.Abs2());
  double x = m_x;
  double y = m_kt2/(m_Q2*x);
  double alpha, beta;
  if (m_mspect2<1.e-12) {
    beta  = y+(m_mbeam2+m_kt2)/((1.-x)*m_Q2);
    alpha = 0.;
  }
  else {
    double h1   = - (1.+y + (m_mbeam2+m_kt2-m_mspect2)/((1.-x)*m_Q2));
    double h2   =       y + (m_mbeam2+m_kt2-y*m_mspect2)/((1.-x)*m_Q2);
    double disc = h1*h1-4.*h2;
    if (disc<0.) return false;
    beta  = -(h1+sqrt(disc))/2.;
    alpha = m_mspect2/(m_Q2*(1.-beta));
  }
  if (dabs(beta-y)<1.e-12) beta=y;
  //msg_Out()<<METHOD<<"(x = "<<x<<", 1-y = "<<(1-y)<<", kt2 = "<<m_kt2<<"): "
  //	   <<"alpha, beta = "<<alpha<<", "<<beta<<"\n";
  // Simple kinematics checks on parameters
  if (y<0 || y>beta || alpha<0. || alpha>1.+1.e-6) {
    //msg_Out()<<"      ---> Limits killed it.\n";
    return false;
  }
  double E = m_Q/2.;
  pi       = E * ((1.-alpha-x) * s_AxisP + (   beta-y) * s_AxisM) + m_ktvec;
  pj       = E * (          x  * s_AxisP +          y  * s_AxisM) - m_ktvec;
  pk       = E * (    alpha    * s_AxisP + (1.-beta  ) * s_AxisM);
  // Simple kinematics checks on momenta: right mass shells, enough mass
  // in the beam-gluon pair positive energies for all particles
  if (!IsZero(dabs(pi.Abs2())-m_mbeam2) || !IsZero(dabs(pj.Abs2())) ||
      !IsZero(dabs(pk.Abs2())-m_mspect2) ||
      !IsZero(dabs((pi+pj+pk).Abs2()) - m_Q2) ||
      !(pi[0]>0.) || !(pj[0]>0.) || !(pk[0]>0.) ||
      (pi+pj).Abs2()<m_minMbeam2 || (pj+pk).Abs2()<m_minMspect2) {
    //msg_Out()<<"      ---> Vector checks:\n"
    //	     <<"      pi,j,k = "
    //	     <<dabs(pi.Abs2())-m_mbeam2<<", "
    //	     <<dabs(pj.Abs2())<<", "
    //	     <<dabs(pk.Abs2())-m_mspect2<<"\n"
    //	     <<"      Q2 = "<<IsZero(dabs((pi+pj+pk).Abs2()) - m_Q2)<<", "
    //	     <<"mins = "<<((pi+pj).Abs2()>m_minMbeam2)<<" (beam), "
    //	     <<((pj+pk).Abs2()>m_minMspect2)<<" (spect) : "
    //	     <<(pj+pk).Abs2()<<" > "<<m_minMspect2<<"\n";
    return false;
  }
  return true;
}

bool Beam_Decorrelator::
FillSoftEmission(ATOOLS::Vec4D & pi,ATOOLS::Vec4D & pj,ATOOLS::Vec4D & pk) {
  // Rotating and boosting back
  m_rotat.RotateBack(pi);
  m_rotat.RotateBack(pj);
  m_rotat.RotateBack(pk);
  m_boost.BoostBack(pi);
  m_boost.BoostBack(pj);
  m_boost.BoostBack(pk);
  // Updating the momenta of beam parton (the emitter) and spectator,
  // adjusting colours
  p_beam->SetMomentum(pi);
  p_spect->SetMomentum(pk);
  size_t index = (p_beam->GetFlow(1)>0 && p_spect->GetFlow(2)>0)?1:2;
  //msg_Out()<<METHOD<<":\n"<<(*p_beam)<<"\n"<<(*p_spect)<<"\n";
  p_beam->SetFlow(index,-1);
  // adding a soft gluon, adjusting its number and colours
  Particle * gluon = new Particle(-1,Flavour(kf_gluon),pj,'B');
  gluon->SetFlow(3-index,p_beam->GetFlow(index));
  gluon->SetFlow(index,p_spect->GetFlow(3-index));
  gluon->SetNumber();
  //msg_Out()<<"Now:\n"<<(*p_beam)<<"\n"<<(*p_spect)<<"\n"<<(*gluon)<<"\n";
  m_softgluons.push_back(gluon);
  return true;
}

void Beam_Decorrelator::Reset() {
  m_softgluons.clear();
}


