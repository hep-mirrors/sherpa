#include "CFPSHOWER++/Calculators/FI/SF_FI12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_FI12::SF_FI12(const Kernel_Info & info) : SF_Base(info) {}

double SF_FI12::Jacobean(const Splitting & split) const {
  double xB    = split.GetSpectator()->XB();
  double ymass = split.y()* (1. + (split.msplit2()-split.m2(0)-split.m2(1))/split.Q2() );
  // added a security instance for the ratio of xB and modified splitting parameter
  if (xB/ymass > 1.) return 0.;
  // remember: IS particles are "barred"
  const Flavour & flav = split.GetSpectator()->Flav().Bar();
  // pdf estimates for initial state spectator before and after splitting
  size_t beam    = split.GetSpectator()->Beam()-1;
  double pdf_old = split.GetKernel()->GetXPDF(      xB,split.t(),flav,beam);
  double pdf_new = split.GetKernel()->GetXPDF(xB/ymass,split.t(),flav,beam);
  // return a non-zero Jacobean only if splitting is kinematically allowed
  if (dabs(pdf_old) > (split.GetKernel()->PDFMinValue() *
		       log(1.-xB)/log(1.-split.GetKernel()->PDFXMin())) ) {
    return (1.-split.y())/(1.-ymass) * pdf_new/pdf_old;
  }
  return 0.;
}

bool SF_FI12::Construct(Splitting & split,const int & mode) {
  double y = 1. - split.t()/(split.Q2red()*(1.-split.z()));
  Kin_Args kin_args(1.-y,split.z(),split.phi(),1|8);
  if (ConstructFIDipole(split.m2(0),split.m2(1),
			split.msplit2(),split.mspect2(),
			split.GetSplitter()->Mom(),-split.GetSpectator()->Mom(),
			kin_args) < 0) return false;
  split.Set_y(y);
  m_moms[0] = kin_args.m_pi;
  m_moms[1] = kin_args.m_pj;
  m_specmom = kin_args.m_pk;
  return true;
}
