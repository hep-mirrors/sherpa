#include "CFPSHOWER++/Calculators/FF/SF_FF13.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_FF13::SF_FF13(const Kernel_Info & info) : SF_Base(info) {
  m_subtract = subtract::none;
  if (m_flavs[m_tags[1]].IsGluon() ||
      m_flavs[m_tags[1]]==m_flavs[m_tags[0]].Bar()) {
    m_subtract = m_subtract | subtract::coll;
  }
  m_moms.resize(3);
  /*
  msg_Out()<<"     * "<<METHOD<<"[subtract = "<<m_subtract<<"] "
	   <<"for "<<m_split<<" -> ";
  for (size_t i=0;i<m_flavs.size();i++) msg_Out()<<m_flavs[i]<<" ";
  msg_Out()<<" [ ";
  for (size_t i=0;i<m_flavs.size();i++) msg_Out()<<m_tags[i]<<" ";
  msg_Out()<<"], (compared "<<m_split<<" with "<<m_flavs[m_tags[0]]
	   <<" ["<<m_tags[0]<<"])\n";
  */
}

double SF_FF13::Jacobean(const Splitting & split) const {
  double maij2 = split.msplit2(), mk2 = split.mspect2();
  double Q2red = split.Q2()-maij2-mk2;
  double Jac1  = Q2red/Lambda(split.Q2(),maij2,mk2);
  double sai   = split.sai(), xi = split.z(0)/split.z(1);
  double saik  = xi * Q2red + sai + mk2;
  double Jac2  = (saik-sai-mk2)/Lambda(saik,sai,mk2);
  double total = (Jac1 * Jac2) / (1. + xi*(sai+split.m2(m_tags[2])-maij2)/split.t(0));
  return total;
}

bool SF_FF13::Construct(Splitting & split,const int & mode) {
  // Logic two consequtive splittings:
  // kin_args1  = (ij)_0 + k_0  ->  i_1 + j_1 + k_1 with i_1 = ij_2, massive
  // kin_args2  = i_1    + k_1  ->  i_2 + j_2 + k_1
  //              since i_1 is massive, with the right invariant mass, the
  //              spectator momentum in this splitting (k_1) will not change,
  //              and is only needed to fix the angles.
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
  Kin_Args kin_args1(y,xtilde,split.phi(0));
  // msg_Out()<<METHOD<<": "<<m_split<<" --> "
  // 	   <<m_flavs[0]<<" "<<m_flavs[1]<<" "<<m_flavs[2]<<" "
  // 	   <<"["<<m_tags[0]<<" "<<m_tags[1]<<" "<<m_tags[2]<<"]\n"
  // 	   <<"   *** for y = "<<y<<", x = "<<xtilde<<", "
  // 	   <<"phi = "<<split.phi()<<", z = {"<<split.z()<<", "<<split.z2()<<"}, "
  // 	   <<"sai = "<<split.t2()<<", "<<split.Mode()<<", "
  // 	   <<"for t = "<<split.t()<<".\n";
  if (ConstructFFDipole(sai,mj2,
			maij2,mk2,
			split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
			kin_args1) < 0) return false;
  Vec4D  pai =  kin_args1.m_pi, pk = kin_args1.m_pk;
  double ma2 = split.m2(m_tags[0]), mi2 = split.m2(m_tags[1]);
  double y2  = ((sai<1.e-12) ? 0. : 1./(1.+(2.*pai*pk) / (sai-ma2-mi2)));
  //msg_Out()<<"   *** "<<METHOD<<" for sai = "<<split.t2()<<" -> "<<sai<<", "
  //	   <<"y2 = "<<y2<<", z2 = "<<split.z2()<<", phi = "<<split.phi2()<<".\n";
  Kin_Args kin_args2(y2,z2,split.phi(1));
  if (ConstructFFDipole(ma2,mi2,
			sai,mk2,
			pai,pk,
			kin_args2) < 0) return false;
  m_specmom = kin_args1.m_pk; // spectator momentum
  m_moms[m_tags[2]] = kin_args1.m_pj; // momentum of parton j
  m_moms[m_tags[0]] = kin_args2.m_pi; // momentum of "tagged" parton a
  m_moms[m_tags[1]] = kin_args2.m_pj; // momentum of parton i
  // msg_Out()<<"   *** p[a] = "<<m_moms[m_tags[0]]<<"\n"
  // 	   <<"   *** p[i] = "<<m_moms[m_tags[1]]<<"\n"
  // 	   <<"   *** p[j] = "<<m_moms[m_tags[2]]<<"\n"
  // 	   <<"   *** pk   = "<<m_specmom<<"\n";
  return true;
}
