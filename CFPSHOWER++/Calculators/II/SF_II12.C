#include "CFPSHOWER++/Calculators/II/SF_II12.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_II12::SF_II12(const Kernel_Info & info) : SF_Base(info) {}

double SF_II12::Jacobean(const Splitting & split) const {
  double J = 1.;
  return J;
}

bool SF_II12::InitKinematics(Splitting & split) const {
  split.Set_y(split.t() / (split.Q2() * (1.-split.z(0))));
  split.Set_x((split.z(0)-split.y()) / (1.-split.y()));
  return true;
}

int SF_II12::Construct(Splitting & split) const {
  Kin_Args kin_ii(split.y(),split.x(),split.phi());
  if (ConstructIIDipole(split.m2(0),split.m2(1),
			split.msplit2(),split.mspect2(),
			split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
			kin_ii) < 0) return -1;
  split.SetMom(0, kin_ii.m_pi);
  split.SetMom(1, kin_ii.m_pj);
  split.SetSpecMom(kin_ii.m_pk);
  return 1;
}
