#include "CFPSHOWER++/Calculators/FF/SF_FF13.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_FF13::SF_FF13(const Kernel_Info & info) : SF_Base(info) {}

double SF_FF13::Jacobean(const Splitting & split) const {
  double Q2    = split.Q2(), t = split.t();
  double maij2 = split.msplit2(), mk2 = split.mspect2(), mj2 = split.m2(2);
  double Jac1  = (Q2-maij2-mk2)/Lambda(Q2,maij2,mk2);
  double sai   = split.t2(), za = split.z(), xa = split.z2();
  double saik  = za/xa * (Q2-maij2-mk2) + sai + mk2;
  double Jac2  = (saik-sai-mk2)/Lambda(saik,sai,mk2);
  double total = (Jac1 * Jac2) / (1. + za/xa*(sai+mj2-maij2)/t);
}

bool SF_FF13::Construct(Splitting & split) {
  // Logic two consequtive splittings:
  // kin_args1  = (ij)_0 + k_0  ->  i_1 + j_1 + k_1 with i_1 = ij_2, massive
  // kin_args2  = i_1    + k_1  ->  i_2 + j_2 + k_1
  //              since i_1 is massive, with the right invariant mass, the spectator
  //              momentum in this splitting (k_1) will not change, and is only needed
  //              to fix the angles.
  double mtot = 0.;
  for (size_t i=0;i<3;i++) mtot += sqrt(split.m2(i));
  if (sqrt(split.Q2())<mtot+sqrt(split.mspect2())) return false;
  double Q2     = split.Q2(), sai   = split.t2(), z2 = split.z2(), z1 = split.z();
  double maij2  = split.msplit2(), mk2 = split.mspect2(), mj2 = split.m2(2);
  double Q2red  = Q2 - maij2 - mk2;
  double Q2p    = Q2 - sai - mj2 - mk2;
  double y      = (z1*split.t()) / (z2*Q2p);
  double xtilde = (z2*Q2) / (z1*(1.-y)*Q2p);
  split.Set_y(y);
  Kin_Args kin_args1(y,xtilde,split.phi());
  if (ConstructFFDipole(sai,mj2,
			maij2,mk2,
			split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
			kin_args1) < 0) return false;
  Vec4D pai    =  kin_args1.m_pi, pk = kin_args1.m_pk;
  double ma2   = split.m2(0), mi2 = split.m2(1);
  double y2    = ((sai<1.e-12) ? 0. : 1./(1.+(2.*pai*pk) / (sai-mi2)));
  Kin_Args kin_args2(y2,z2,split.phi2());
  if (ConstructFFDipole(ma2,mi2,
			sai,mk2,
			pai,pk,
			kin_args2)<0) return false;
  m_moms[0] = kin_args2.m_pi;
  m_moms[1] = kin_args2.m_pj;
  m_moms[2] = kin_args1.m_pj;
  m_specmom = kin_args1.m_pk;
  return true;
}
