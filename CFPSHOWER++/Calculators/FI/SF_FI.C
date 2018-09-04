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
  double J     = 0.;
  double eta   = split.GetSpectator()->XB();
  double ymass = split.Y()* (1.0 + (split.mij2()-split.mi2()-split.mj2())/split.Q2() );
  // remember: IS particles are "barred"
  const Flavour & flav = split.GetSpectator()->Flav().Bar();
  size_t beam          = split.GetSpectator()->Beam()-1;
  double f_old = split.GetKernel()->GetXPDF(eta,split.T(),flav,beam);
  double f_new = split.GetKernel()->GetXPDF(eta/ymass,split.T(),flav,beam);
  if (dabs(f_old) > (split.GetKernel()->PDFMinValue() *
		     log(1.-eta)/log(1.-split.GetKernel()->PDFXMin())) ) {
    J = (1.-split.Y())/(1.-ymass) * f_new/f_old;
  }
  else {
  }
  return J;
}

bool SF_FI::InitKinematics(Splitting & split) const {
  split.SetY(1. - split.T()/(split.Q2()*(1.-split.Z())));
  split.SetX(split.Z());
  return true;
}

int SF_FI::Construct(Splitting & split) const {
  //msg_Out()<<" *** "<<METHOD<<"(t = "<<sqrt(split.T())<<", z = "<<split.Z()<<") "
  //	   <<"--> y, x ="<<split.Y()<<", "<<split.X()<<" and phi = "<<split.phi()
  //	   <<" for "<<m_name<<"\n"
  //	   <<"     masses = "<<split.mi2()<<", "<<split.mj2()<<", "
  //	   <<split.mij2()<<", "<<split.mk2()<<" "   	   
  //	   <<"and Q2 from momenta = "
  //	   <<(split.GetSplitter()->Mom()+split.GetSpectator()->Mom()).Abs2()<<"\n";
  Kin_Args kin_fi(1.-split.Y(),split.X(),split.phi(),1|8);
  if (ConstructFIDipole(split.mi2(),split.mj2(),split.mij2(),split.mk2(),
			split.GetSplitter()->Mom(),-split.GetSpectator()->Mom(),
			kin_fi) < 0) return -1;
  split.SetMom(0, kin_fi.m_pi);
  split.SetMom(1, kin_fi.m_pj);
  split.SetSpecMom(-kin_fi.m_pk);
  return 1;
}
