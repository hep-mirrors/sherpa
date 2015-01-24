#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class SingleTop_Core_Scale: public Core_Scale_Setter {
  public:

    SingleTop_Core_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::CParam SingleTop_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  double shat(2.0*ampl->Leg(0)->Mom()*ampl->Leg(1)->Mom());
  double that(2.0*ampl->Leg(0)->Mom()*ampl->Leg(2)->Mom());
  double uhat(2.0*ampl->Leg(0)->Mom()*ampl->Leg(3)->Mom());
  if (ampl->Legs().size()>4) {
    //    msg_Out()<<"Not enough clustering, so we are left with n legs= "<<ampl->Legs().size()<<std::endl;
    return PDF::CParam(shat,shat,0.0,shat,-1);
  }
  double muf2(-1.),mur2(-1.),muq2(-1.);
  Flavour f[4]={ampl->Leg(0)->Flav(),ampl->Leg(1)->Flav(),
		ampl->Leg(2)->Flav(),ampl->Leg(3)->Flav()};
   // tW production: muF = mT(top)
  if ((f[2].Kfcode()==24 && f[3].Kfcode()==6) ||
      (f[2].Kfcode()==6 && f[3].Kfcode()==24)) {
    muf2 = mur2 = muq2 = sqr(Flavour(kf_t).Mass()) +
      (f[2].Kfcode()==6?
       ampl->Leg(2)->Mom().PPerp2():ampl->Leg(3)->Mom().PPerp2());
  }
  // s-channel top: muF = \hat{s} 
  else if ((f[2].Kfcode()==5 && f[3].Kfcode()==6) ||
	   (f[2].Kfcode()==6 && f[3].Kfcode()==5)) {
    muf2 = mur2 = muq2 = shat;
  }
  // s-channel top, clustered into V^*+jet
  else if ((f[2].Kfcode()==24 && f[3].Strong() && f[3].QuarkFamily()!=3) ||
	   (f[2].Strong() && f[2].QuarkFamily()!=3 && f[3].Kfcode()==24)) {
    muf2 = muq2 = 
      (f[2].Kfcode()==23?
       ampl->Leg(2)->Mom().MPerp2():ampl->Leg(3)->Mom().MPerp2());
    mur2 =       
      (f[2].Kfcode()!=23?
       ampl->Leg(2)->Mom().PPerp2():ampl->Leg(3)->Mom().PPerp2());
  }
  // t-channel top: muF = \hat{t} 
  else if ((f[2].IsQuark() && f[2].QuarkFamily()!=3 && f[3].Kfcode()==6) ||
	   (f[2].Kfcode()==6 && f[3].IsQuark() && f[3].QuarkFamily()!=3)) {
    muf2 = mur2 = muq2 = dabs(that);
  }
  // top in initial state
  else if ((f[0].Kfcode()==6 && f[1].Strong() && f[1].Kfcode()!=6) ||
	   f[0].Strong() && f[0].Kfcode()!=6 && f[1].Kfcode()==6) {
    unsigned short int nottop=1;
    if (f[1].Kfcode()==6){ 
      if (f[0].Kfcode()==6) std::cout<<"tt initial state\n";
      nottop=0;
    }
    muf2 = muq2 = mur2 = (f[nottop].Kfcode()==5?shat:dabs(that));
  }
  // W in initial state
  else if ((f[0].Kfcode()==24 && f[1].Strong()) ||
	   (f[0].Strong() && f[1].Kfcode()==24)) {
    muf2 = mur2 = muq2 = dabs(that*uhat)/shat;
  }
  // pure QCD process
  else if (f[0].Strong() && f[0].QuarkFamily()!=3 &&
	   f[1].Strong() && f[1].QuarkFamily()!=3 && 
	   f[2].Strong() && f[3].Strong()) {
    muf2 = mur2 = muq2 = -1.0/(1.0/shat+1.0/that+1.0/uhat);    
  }
  else if ((f[0].Strong() && f[1].Strong()) && 
          (f[0].Kfcode()==5 || f[1].Kfcode()==5)
	   && f[2].Strong() && f[3].Strong())
    muf2 = mur2 = muq2 = -1.0/(1.0/shat+1.0/that+1.0/uhat);    
  if (muf2<0.) {
    msg_Out()<<METHOD<<": found something unexpected: "
	     <<f[0]<<" "<<f[1]<<" --> "<<f[2]<<" "<<f[3]<<",\n"
	     <<"   pt of jet = "<<ampl->Leg(3)->Mom().PPerp()<<" vs. "
	     <<"mass = "<<ampl->Leg(2)->Mom().Abs()<<".\n";
    exit(1);
  }
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(muf2)<<"\n"
		 <<"  \\mu_r = "<<sqrt(mur2)<<"\n"
		 <<"  \\mu_q = "<<sqrt(muq2)<<"\n";
  msg_Debugging()<<"}\n";
  return PDF::CParam(muf2,muq2,0.0,mur2,-1);
}

DECLARE_ND_GETTER(SingleTop_Core_Scale,"SingleTop",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,SingleTop_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  // check for oew == 2 or so, to force single top!
  return new SingleTop_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    SingleTop_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"single top core scale"; 
}
