#include "CFPSHOWER++/Calculators/FF/SF_FF12.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_FF12::SF_FF12(const Kernel_Info & info) : SF_Base(info) {}

double SF_FF12::Jacobean(const Splitting & split) const {
  double Q2 = split.Q2(), Q2y = Q2*split.y();
  return ( Q2  / Lambda(split.s(),split.msplit2(),split.mspect2()) *
	   Q2y / (Q2y + split.m2(0)+split.m2(1)-split.msplit2()) );
}

bool SF_FF12::InitKinematics(Splitting & split) const {
  if (sqrt(split.s())<sqrt(split.m2(0))+sqrt(split.m2(1))+sqrt(split.mspect2())) return false;
  split.Set_y(split.t() / (split.Q2() * (1.-split.z(0))));
  split.Set_x((split.z(0)-split.y()) / (1.-split.y()));
  return true;
}

int SF_FF12::Construct(Splitting & split) const {
  Kin_Args kin_ff(split.y(),split.x(),split.phi());
  if (ConstructFFDipole(split.m2(0),split.m2(1),
			split.msplit2(),split.mspect2(),
			split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
			kin_ff) < 0) return -1;
  split.SetMom(0, kin_ff.m_pi);
  split.SetMom(1, kin_ff.m_pj);
  split.SetSpecMom(kin_ff.m_pk);
  return 1;
}
