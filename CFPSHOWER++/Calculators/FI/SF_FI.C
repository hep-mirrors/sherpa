#include "CFPSHOWER++/Calculators/FI/SF_FI.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_FI::SF_FI(const Kernel_Info & info) : SF_Base(info) {}

double SF_FI::Jacobean(const Splitting & split) const {
  double xB    = split.GetSpectator()->XB();
  double ymass = split.Y()* (1. + (split.mij2()-split.mi2()-split.mj2())/split.Q2() );
  // added a security instance for the ratio of xB and modified splitting parameter
  if (xB/ymass > 1.) return 0.;
  // remember: IS particles are "barred"
  const Flavour & flav = split.GetSpectator()->Flav().Bar();
  // pdf estimates for initial state spectator before and after splitting
  size_t beam    = split.GetSpectator()->Beam()-1;
  double pdf_old = split.GetKernel()->GetXPDF(      xB,split.T(),flav,beam);
  double pdf_new = split.GetKernel()->GetXPDF(xB/ymass,split.T(),flav,beam);
  // return a non-zero Jacobean only if splitting is kinematically allowed
  if (dabs(pdf_old) > (split.GetKernel()->PDFMinValue() *
		       log(1.-xB)/log(1.-split.GetKernel()->PDFXMin())) ) {
    return (1.-split.Y())/(1.-ymass) * pdf_new/pdf_old;
  }
  return 0.;
}

bool SF_FI::InitKinematics(Splitting & split) const {
  split.SetY(1. - split.T()/(split.Q2()*(1.-split.Z())));
  split.SetX(split.Z());
  return true;
}

int SF_FI::Construct(Splitting & split) const {
  Kin_Args kin_fi(1.-split.Y(),split.X(),split.phi(),1|8);
  if (ConstructFIDipole(split.mi2(),split.mj2(),split.mij2(),split.mk2(),
			split.GetSplitter()->Mom(),-split.GetSpectator()->Mom(),
			kin_fi) < 0) return -1;
  split.SetMom(0, kin_fi.m_pi);
  split.SetMom(1, kin_fi.m_pj);
  split.SetSpecMom(-kin_fi.m_pk);
  return 1;
}
