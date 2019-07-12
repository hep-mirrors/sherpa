#include "CFPSHOWER++/Calculators/FF/SF_FF.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_FF::SF_FF(const Kernel_Info & info) : SF_Base(info) {}

double SF_FF::Jacobean(const Splitting & split) const {
  double Q2 = split.Q2(), Q2y = Q2*split.Y();
  return ( Q2  / Lambda(split.sijk(),split.mij2(),split.mk2()) *
	   Q2y / (Q2y + split.mi2()+split.mj2()-split.mij2()) );
}

bool SF_FF::InitKinematics(Splitting & split) const {
  if (sqrt(split.sijk())<sqrt(split.mi2())+sqrt(split.mj2())+sqrt(split.mk2())) return false;
  split.SetY(split.T() / (split.Q2() * (1.-split.Z())));
  split.SetX((split.Z()-split.Y()) / (1.-split.Y()));
  return true;
}

int SF_FF::Construct(Splitting & split) const {
  Kin_Args kin_ff(split.Y(),split.X(),split.phi());
  if (ConstructFFDipole(split.mi2(),split.mj2(),split.mij2(),split.mk2(),
			split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
			kin_ff) < 0) return -1;
  split.SetMom(0, kin_ff.m_pi);
  split.SetMom(1, kin_ff.m_pj);
  split.SetSpecMom(kin_ff.m_pk);
  return 1;
}
