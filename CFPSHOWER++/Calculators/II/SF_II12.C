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

bool SF_II12::Construct(Splitting & split,const int & mode) {
  double z = split.z(0), y = split.t(0) / (split.Q2red() * (1.-z));
  double xtilde = (z-y) / (1.-y);
  Kin_Args kin_args(y,xtilde,split.phi(0));
  if (ConstructIIDipole(split.m2(0),split.m2(1),
			split.msplit2(),split.mspect2(),
			split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
			kin_args) < 0) return false;
  split.Set_y(y);
  m_moms[0] = kin_args.m_pi;
  m_moms[1] = kin_args.m_pj;
  m_specmom = kin_args.m_pk;
  return true;
}
