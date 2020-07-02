#include "CFPSHOWER++/Calculators/FF/SF_FF3_Coll.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace ATOOLS;

bool SF_FF3_Coll::operator()(Splitting & split,Configuration & config) {
  Init(split);
  return ConstructSystem(split);
}

bool SF_FF3_Coll::KinCheck(Splitting & split) {
  return (split.y()>=0.0 && split.y()<=1.0);
}

double SF_FF3_Coll::CalculateY(Splitting & split) {
  if (split.GetSplitter()->Flav().IsFermion())
    return split.t(0)/((1.-split.z(0))*split.Q2());
  if (split.GetSplitter()->Flav().IsVector()) {
    if (split.GetKernel()->GetFlavs()[0].IsVector())
      return split.t(0)/((1.-split.z(0))*split.z(0)*split.Q2());
    if (split.GetKernel()->GetFlavs()[0].IsFermion())
      return split.t(0)/split.Q2();
  }
  return split.t(0)/((1.-split.z(0))*split.z(0)*split.Q2());
}

bool SF_FF3_Coll::ConstructSystem(Splitting & split) {
  double mtot = 0.;
  for (size_t i=0;i<3;i++) mtot += sqrt(split.m2(i));
  if (sqrt(split.Q2())<mtot+sqrt(split.mspect2())) {
    return false;
  }
  double Q2     = split.Q2(), z1 = split.z(0), z2 = split.z(1), xi = z1/z2;
  double sai    = (mode&1 && split.IsEndPoint())? 0.: split.sai(); 
  double maij2  = split.msplit2(), mk2 = split.mspect2(), mj2 = split.m2(m_tags[2]);
  double Q2red  = Q2 - maij2 - mk2;
  double Q2p    = Q2 - sai - mj2 - mk2;
  double y      = split.t(0)/ (xi*Q2p);
  double xtilde = (xi*Q2red) / ((1.-y)*Q2p);
  split.Set_y(y);
  PHASIC::Kin_Args kin_args1(y,xtilde,split.phi(0));
  if (PHASIC::ConstructFFDipole(sai,mj2,
				maij2,mk2,
				split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
				kin_args1) < 0) return false;
  Vec4D  pai =  kin_args1.m_pi, pk = kin_args1.m_pk;
  double ma2 = split.m2(m_tags[0]), mi2 = split.m2(m_tags[1]);
  double y2  = ((sai<1.e-12) ? 0. : 1./(1.+(2.*pai*pk) / (sai-ma2-mi2)));
  PHASIC::Kin_Args kin_args2(y2,z2,split.phi(1));
  if (PHASIC::ConstructFFDipole(ma2,mi2,
				sai,mk2,
				pai,pk,
				kin_args2) < 0) return false;
  split.Set_ztilde(0,split.z(0));
  split.Set_ztilde(1,1.-split.z(0));
  split.SetMomentum(m_tags[0],kin_args2.m_pi);
  split.SetMomentum(m_tags[1],kin_args2.m_pj);
  split.SetMomentum(m_tags[2],kin_args1.m_pj);
  split.SetSpectatorMomentum(kin_args1.m_pk);
  return true;
}

void SF_FF3_Coll::CalculateWeight(Splitting & split) {
  double maij2 = split.msplit2(), mk2 = split.mspect2();
  double Q2red = split.Q2()-maij2-mk2;
  double Jac1  = Q2red/Lambda(split.Q2(),maij2,mk2);
  double sai   = split.sai(), xi = split.z(0)/split.z(1);
  double saik  = xi * Q2red + sai + mk2;
  double Jac2  = (saik-sai-mk2)/Lambda(saik,sai,mk2);
  m_weight = (Jac1 * Jac2) / (1. + xi*(sai+split.m2(m_tags[2])-maij2)/split.t(0));
}

bool SF_FF3_Coll::UpdateSystem(Splitting & split,Configuration & config) {
  split.GetSpectator()->SetMom(split.SpectatorMomentum());
  return true;
}


