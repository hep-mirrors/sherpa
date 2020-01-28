#include "CFPSHOWER++/Calculators/FF/SF_FF12.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_FF12::SF_FF12(const Kernel_Info & info) : SF_Base(info) {
  m_moms.resize(2);
}

double SF_FF12::Jacobean(const Splitting & split) const {
  double Q2red = split.Q2red(), Q2red_y = Q2red * split.y(); 
  return ( Q2red / Lambda(split.Q2(),split.msplit2(),split.mspect2()) *
	   Q2red_y /(Q2red_y+split.m2(0)+split.m2(1)-split.msplit2()) );
}

bool SF_FF12::Construct(Splitting & split,const int & mode) {
  if (sqrt(split.Q2()) <
      sqrt(split.m2(0))+sqrt(split.m2(1))+sqrt(split.mspect2())) return false;
  double z = split.z(), y = split.t()/(split.Q2red()*(1.-z)), ztilde = (z-y)/(1.-y);
  Kin_Args kinargs(y,ztilde,split.phi());
  if (ConstructFFDipole(split.m2(0),split.m2(1),
			split.msplit2(),split.mspect2(),
			split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
			kinargs) < 0) return false;
  split.Set_y(y);
  m_moms[0] = kinargs.m_pi;
  m_moms[1] = kinargs.m_pj;
  m_specmom = kinargs.m_pk;
  return true;
}
