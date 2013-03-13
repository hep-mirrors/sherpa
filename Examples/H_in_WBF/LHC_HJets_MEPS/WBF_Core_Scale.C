#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace MYSTUFF {

  class WBF_Core_Scale: public PHASIC::Core_Scale_Setter {
  public:

    WBF_Core_Scale(const PHASIC::Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class WBF_Core_Scale

}// end of namespace MYSTUFF


using namespace MYSTUFF;
using namespace PHASIC;
using namespace ATOOLS;

PDF::CParam WBF_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  double mh2(ampl->Leg(2)->Mom().Abs2());
  double muf2(mh2), mur2(muf2), q2(muf2);
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(muf2)<<"\n"
		 <<"  \\mu_r = "<<sqrt(mur2)<<"\n"
		 <<"  \\mu_q = "<<sqrt(q2)<<"\n";
  msg_Debugging()<<"}\n";
  return PDF::CParam(muf2,q2,0.0,mur2,-1);
}

DECLARE_ND_GETTER(WBF_Core_Scale,"WBF_Test",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,WBF_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new WBF_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    WBF_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"WBF core scale"; 
}
