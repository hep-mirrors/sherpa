#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class QCD_Core2_Scale: public Core_Scale_Setter {
  public:

    QCD_Core2_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::CParam QCD_Core2_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  double s(2.0*ampl->Leg(0)->Mom()*ampl->Leg(1)->Mom());
  double t(2.0*ampl->Leg(0)->Mom()*ampl->Leg(2)->Mom());
  double u(2.0*ampl->Leg(0)->Mom()*ampl->Leg(3)->Mom());


  double pT2 = ampl->Leg(3)->Mom().PPerp2();
  double deltaRap = fabs(ampl->Leg(2)->Mom().DY(ampl->Leg(3)->Mom()));
  double scale = pT2*exp(0.3*deltaRap);
  
  double muf2(scale), mur2(scale), q2(scale);
  msg_Debugging()<<METHOD<<"(): Set QCD2 scale {\n"
		 <<"  \\mu_f = "<<sqrt(muf2)<<"\n"
		 <<"  \\mu_r = "<<sqrt(mur2)<<"\n"
		 <<"  \\mu_q = "<<sqrt(q2)<<"\n";
  msg_Debugging()<<"}\n";
  return PDF::CParam(muf2,q2,0.0,mur2,-1);
}

DECLARE_ND_GETTER(QCD_Core2_Scale,"QCD2",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,QCD_Core2_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new QCD_Core2_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    QCD_Core2_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"QCD2 core scale"; 
}
