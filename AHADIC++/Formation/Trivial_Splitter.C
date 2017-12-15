#include "AHADIC++/Formation/Trivial_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Trivial_Splitter::Trivial_Splitter() {}

void Trivial_Splitter::Init() {
  p_constituents = hadpars->GetConstituents();
  // minmass is the mass of the lightest constituent (usually u/d quark)
  m_minmass      = p_constituents->MinMass();
  m_flavourselector.InitWeights();
  m_ktselector.Init(false);
  m_zselector.Init();
}

bool Trivial_Splitter::operator()(Proto_Particle * part1,
				  Proto_Particle * part2) {
  p_part1 = part1;
  p_part2 = part2;
  if (!InitKinematics(true)) return false;
  SelectFlavour();
  FixTransverseMomentum(true);
  ConstructRescueMomenta();
  
  p_part1->SetFlavour(m_newflav.Bar());
  p_part1->SetMomentum(m_q1mom);
  p_part2->SetFlavour(m_newflav);
  p_part2->SetMomentum(m_q2mom);
  return true;
}


bool Trivial_Splitter::operator()(Singlet * singlet) {
  p_singlet = singlet;
  list<Proto_Particle *>::iterator ppit1(p_singlet->begin()), ppit2(ppit1);
  ppit2++;
  p_part1    = (*ppit1);
  p_part2    = (*ppit2);
  m_spectmom = p_singlet->back()->Momentum();
  if (!InitKinematics(false)) return Rescue();
  do {
    SelectFlavour();
  } while (!FixTrialKinematics() || !CheckKinematics());

  p_part1->SetFlavour(m_newflav);
  p_part1->SetMomentum(m_q1mom);
  p_part2->SetMomentum(m_glumom);

  p_singlet->push_back(new Proto_Particle(m_newflav.Bar(),m_q2mom,'l'));
  return true;
}

bool Trivial_Splitter::InitKinematics(bool rescue) {
  m_Q2        = (p_part1->Momentum()+p_part2->Momentum()).Abs2();
  m_E         = sqrt(m_Q2)/2.;

  Vec4D mom1  = p_part1->Momentum(), mom2 = p_part2->Momentum();
  m_boost     = Poincare(mom1+mom2);
  m_boost.Boost(mom1);
  m_rotat     = Poincare(mom1,m_E*s_AxisP); 
  if (rescue) return (m_E>m_minmass);
  if (m_E<4.*m_minmass) return false;
  double R2   = 3.*sqr(m_minmass);
  double arg1 = 1.-4.*sqr(m_minmass)/m_Q2;
  double arg2 = sqr(1.-R2/m_Q2)-4.*sqr(m_minmass)/m_Q2;
  double zmin = Max(1.+R2/m_Q2-sqrt(arg2), 1.-sqrt(arg1))/2.;
  double zmax = Min(1.-R2/m_Q2+sqrt(arg2), 1.+sqrt(arg1))/2.;
  
  if (m_E<m_minmass || zmin>zmax || zmax<0.) return false;
  return true;
}

void Trivial_Splitter::SelectFlavour() {
  m_newflav      = m_flavourselector(m_E,true);
  m_popped_mass  = p_constituents->Mass(m_newflav);
  m_popped_mass2 = sqr(m_popped_mass);
}

void Trivial_Splitter::FixTransverseMomentum(bool rescue) {
  // for no transverse momentum replace m_ktmax = 0.
  m_ktmax = rescue? 0.: m_E-m_popped_mass-m_minmass/2.;
  m_kt    = m_ktmax>0.? m_ktselector(m_ktmax) : 0.;
  m_kt2   = m_kt*m_kt;
  m_phi   = 2.*M_PI*ran->Get();
  m_ktvec = m_kt * Vec4D(0.,cos(m_phi),sin(m_phi),0.);
}

bool Trivial_Splitter::FixBetaAndZ() {
  double R2   = (2.*m_popped_mass+m_minmass)*m_minmass;
  double mt2  = m_popped_mass2+m_kt2;
  double arg1 = 1.-4.*mt2/m_Q2;
  double arg2 = sqr(1.-R2/m_Q2)-4.*mt2/m_Q2;
  if (arg1<0. || arg2<0.) {
    msg_Error()<<"Error in "<<METHOD<<"("<<arg1<<", "<<arg2<<") "
	       <<"from "<<sqrt(m_Q2)<<" --> "
	       <<m_popped_mass<<" + "<<m_kt<<"\n";
    return false;
  }
  double zmin = Max(1+R2/m_Q2-sqrt(arg2), 1.-sqrt(arg1))/2.;
  double zmax = Min(1-R2/m_Q2+sqrt(arg2), 1.+sqrt(arg1))/2.;
  m_z    = m_zselector(zmin,zmax);
  m_beta = 1.-mt2/(m_Q2*m_z*(1.-m_z));
  return true;
}

bool Trivial_Splitter::ConstructMomenta() {
  double mt2_tilde = (m_popped_mass2+m_kt2)/m_Q2;
  m_q1mom  = m_E*(     m_z*s_AxisP +      mt2_tilde/m_z*s_AxisM)-m_ktvec;
  m_q2mom  = m_E*((1.-m_z)*s_AxisP + mt2_tilde/(1.-m_z)*s_AxisM)+m_ktvec;
  m_glumom = m_E*m_beta*s_AxisM;
  return true;
}  

bool Trivial_Splitter::FixTrialKinematics() {
  FixTransverseMomentum();
  if (!FixBetaAndZ() || !ConstructMomenta()) return false;
  m_rotat.RotateBack(m_q1mom);
  m_rotat.RotateBack(m_q2mom);
  m_rotat.RotateBack(m_glumom);
  m_boost.BoostBack(m_q1mom);
  m_boost.BoostBack(m_q2mom);
  m_boost.BoostBack(m_glumom);
  return true;
}
  
bool Trivial_Splitter::CheckKinematics() {
  // check if (last quark--gluon) pairing is heavy enough.
  return (sqrt((m_spectmom+m_q2mom).Abs()) > m_popped_mass + 2.*m_minmass);
}

bool Trivial_Splitter::Rescue() {
  // in this case, the invariant mass of the two gluons is below
  // the minimal mass of the lightest quark pair
  if (m_E<m_minmass) return false;
  SelectFlavour();
  FixTransverseMomentum(true);
  ConstructRescueMomenta();
  
  p_part1->SetFlavour(m_newflav.Bar());
  p_part1->SetMomentum(m_q1mom);
  p_part2->SetFlavour(m_newflav);
  p_part2->SetMomentum(m_q2mom);

  p_singlet->push_back(p_singlet->front());
  p_singlet->pop_front();
  return true;
}

void Trivial_Splitter::ConstructRescueMomenta() {
  double z = 0.5*(1.+sqrt(1.-4.*(m_kt2+m_popped_mass2)/m_Q2));
  m_q1mom  = m_E*(     z*s_AxisP + (1.-z)*s_AxisM)-m_ktvec;
  m_q2mom  = m_E*((1.-z)*s_AxisP +      z*s_AxisM)+m_ktvec;
  m_rotat.RotateBack(m_q1mom);
  m_rotat.RotateBack(m_q2mom);
  m_boost.BoostBack(m_q1mom);
  m_boost.BoostBack(m_q2mom);
}
