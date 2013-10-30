#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class TTBar_Core_Scale: public Core_Scale_Setter {
  public:

    TTBar_Core_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::CParam TTBar_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  double s(2.0*ampl->Leg(0)->Mom()*ampl->Leg(1)->Mom());
  double t(2.0*ampl->Leg(0)->Mom()*ampl->Leg(2)->Mom());
  double u(2.0*ampl->Leg(0)->Mom()*ampl->Leg(3)->Mom());
  if (ampl->Legs().size()>4) {
    msg_Tracking()<<METHOD<<"(): 2->"<<ampl->Legs().size()-2
		  <<" process. Returning \\hat{s}."<<std::endl;
    return PDF::CParam(s,s,0.0,s,-1);
  }
  double muf2(-1.0);
  Flavour f[4]={ampl->Leg(0)->Flav(),ampl->Leg(1)->Flav(),
		ampl->Leg(2)->Flav(),ampl->Leg(3)->Flav()};
  if (f[0].IsGluon() && f[1].IsGluon() &&
      f[2].IsQuark() && f[3]==f[2].Bar()) {
    double m2(sqr(f[2].Mass()));
    double Mt(dabs(1.0/6.0*(t*u-m2*(4.0*(m2+t)+m2*t/s))/(t*t)));
    double Mu(dabs(1.0/6.0*(u*t-m2*(4.0*(m2+u)+m2*u/s))/(u*u)));
    double disc((Mt+Mu)*ran->Get());
    muf2=dabs(disc>(f[3].IsAnti()?Mt:Mu)?u:t);
  }
  else if (f[0].IsQuark() && f[1]==f[0].Bar() &&
	   f[2].IsQuark() && f[3]==f[2].Bar() &&
	   f[0].Kfcode()!=f[2].Kfcode()) {
    muf2=dabs((f[0].IsAnti()!=f[2].IsAnti())?t:u);
  }
  else {
    THROW(fatal_error,"Invalid call");
  }
  double mur2(muf2), muq2(muf2);
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(muf2)<<"\n"
		 <<"  \\mu_r = "<<sqrt(mur2)<<"\n"
		 <<"  \\mu_q = "<<sqrt(muq2)<<"\n";
  msg_Debugging()<<"}\n";
  return PDF::CParam(muf2,muq2,0.0,mur2,-1);
}

DECLARE_ND_GETTER(TTBar_Core_Scale,"TTBar",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,TTBar_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new TTBar_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    TTBar_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"tt~ core scale"; 
}
