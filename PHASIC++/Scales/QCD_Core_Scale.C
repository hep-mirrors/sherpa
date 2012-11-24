#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class QCD_Core_Scale: public Core_Scale_Setter,
			     public ATOOLS::Tag_Replacer {
  public:

    QCD_Core_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::CParam QCD_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  double s(2.0*ampl->Leg(0)->Mom()*ampl->Leg(1)->Mom());
  double t(2.0*ampl->Leg(0)->Mom()*ampl->Leg(2)->Mom());
  double u(2.0*ampl->Leg(0)->Mom()*ampl->Leg(3)->Mom());
  double muf2(-1.0/(1.0/s+1.0/t+1.0/u)), mur2(muf2);
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(muf2)<<"\n"
		 <<"  \\mu_r = "<<sqrt(mur2)<<"\n";
  msg_Debugging()<<"}\n";
  return PDF::CParam(muf2,muf2,0.0,mur2,-1);
}

namespace PHASIC {

  DECLARE_ND_GETTER(QCD_Core_Scale_Getter,"QCD",
		    Core_Scale_Setter,Core_Scale_Arguments,true);

  Core_Scale_Setter *QCD_Core_Scale_Getter::operator()
    (const Core_Scale_Arguments &args) const
  {
    return new QCD_Core_Scale(args);
  }

  void QCD_Core_Scale_Getter::PrintInfo
  (std::ostream &str,const size_t width) const
  { 
    str<<"QCD core scale"; 
  }

}// end of namespace PHASIC
