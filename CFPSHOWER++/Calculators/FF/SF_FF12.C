#include "CFPSHOWER++/Calculators/FF/SF_FF12.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_FF12::SF_FF12(const Kernel_Info & info) : SF_Base(info) {}

double SF_FF12::Jacobean(const Splitting & split) const {
  double Q2red = split.Q2red(), Q2red_y = Q2red * split.y(); 
  return ( Q2red / Lambda(split.Q2(),split.msplit2(),split.mspect2()) *
	   1. / (Q2red_y+split.m2(0)+split.m2(1)-split.msplit2())/
		       Q2red_y );
}

bool SF_FF12::Construct(Splitting & split) {
  if (sqrt(split.Q2()) <
      sqrt(split.m2(0))+sqrt(split.m2(1))+sqrt(split.mspect2())) return false;
  m_z = split.z();
  m_y = split.t()/ (split.Q2red() * (1.-m_z));
  m_ztilde = (m_z-m_y) / (1.-m_y);
  Kin_Args kinargs(m_y,m_ztilde,split.phi());
  if (ConstructFFDipole(split.m2(0),split.m2(1),
			split.msplit2(),split.mspect2(),
			split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
			kinargs) < 0) return false;
  split.Set_y(m_y);
  m_moms.resize(2);
  m_moms[0] = kinargs.m_pi;
  m_moms[1] = kinargs.m_pj;
  m_specmom = kinargs.m_pk;
  return true;
}
