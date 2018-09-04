#include "CFPSHOWER++/Calculators/II/SF_II.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_II::SF_II(const Kernel_Info & info) : SF_Base(info) {}

double SF_II::Jacobean(const Splitting & split) const {
  double J = 1.;
  return J;
}

bool SF_II::InitKinematics(Splitting & split) const {
  split.SetY(split.T() / (split.Q2() * (1.-split.Z())));
  split.SetX((split.Z()-split.Y()) / (1.-split.Y()));
  return true;
}

int SF_II::Construct(Splitting & split) const {
  //msg_Out()<<" *** "<<METHOD<<"(t = "<<sqrt(split.T())<<", z = "<<split.Z()<<") "
  //	   <<"--> y, x ="<<split.Y()<<", "<<split.X()<<" and phi = "<<split.phi()
  //	   <<" for "<<m_name<<"\n"
  //	   <<"     masses = "<<split.mi2()<<", "<<split.mj2()<<", "
  //	   <<split.mij2()<<", "<<split.mk2()<<" "   	   
  //	   <<"and Q2 from momenta = "
  //	   <<(split.GetSplitter()->Mom()+split.GetSpectator()->Mom()).Abs2()<<"\n";
  Kin_Args kin_ii(split.Y(),split.X(),split.phi());
  if (ConstructIIDipole(split.mi2(),split.mj2(),split.mij2(),split.mk2(),
			split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
			kin_ii) < 0) return -1;
  split.SetMom(0, kin_ii.m_pi);
  split.SetMom(1, kin_ii.m_pj);
  split.SetSpecMom(kin_ii.m_pk);
  return 1;
}
