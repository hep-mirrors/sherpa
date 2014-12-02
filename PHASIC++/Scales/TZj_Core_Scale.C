#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class TZj_Core_Scale: public Core_Scale_Setter {
  public:

    TZj_Core_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::CParam TZj_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  double s(2.0*ampl->Leg(0)->Mom()*ampl->Leg(1)->Mom());
  double t(2.0*ampl->Leg(0)->Mom()*ampl->Leg(2)->Mom());
  double u(2.0*ampl->Leg(0)->Mom()*ampl->Leg(3)->Mom());
  if (ampl->Legs().size()>4) {
    return PDF::CParam(s,s,0.0,s,-1);
  }
  double muf2(-1.0);
  Flavour f[4]={ampl->Leg(0)->Flav(),ampl->Leg(1)->Flav(),
		ampl->Leg(2)->Flav(),ampl->Leg(3)->Flav()};
  msg_Out()<<"Flavours are: "<<f[0]<<"  "<<f[1]<<"  "<<f[2]<<"  "<<f[3]<<'\n';
  if (f[0].Kfcode()==24 || f[1].Kfcode()==24){
    if (f[2].Kfcode()==6){
      //     muf2=sqr(f[2].Mass()+f[3].Mass()+ampl->Leg(2)->Mom().PPerp());
      muf2=sqr(Flavour(kf_t).Mass());
    }
    else if (f[3].Kfcode()==6){
      //      muf2=sqr(f[3].Mass()+f[2].Mass()+ampl->Leg(3)->Mom().PPerp());
      muf2=sqr(Flavour(kf_t).Mass());
    }
    else{
      if (ampl->Leg(2)->Flav()!=Flavour(kf_Z)) {
	//muf2=sqr(ampl->Leg(2)->Mom().PPerp()+f[3].Mass());
        muf2=sqr(Flavour(kf_t).Mass());
      }
      else {
	//muf2=sqr(ampl->Leg(3)->Mom().PPerp()+f[2].Mass());
	muf2=sqr(Flavour(kf_t).Mass());
      }
    }
  }
  else if (f[2].Kfcode()==24 || f[3].Kfcode()==24){
    if ((f[2].Kfcode()==6) + (f[3].Kfcode()==6)==1){
      //muf2=sqr(f[2].Mass()+f[3].Mass()+(f[2].Kfcode()==6?ampl->Leg(2)->Mom().PPerp():ampl->Leg(3)->Mom().PPerp()));
      muf2=sqr(Flavour(kf_t).Mass());
    }
    else if ((f[2].Kfcode()==23) + (f[3].Kfcode()==23)==1){
      //muf2=sqr(f[2].Mass()+f[3].Mass());
      muf2=sqr(Flavour(kf_t).Mass());
    }
    else{
      if (f[2].Kfcode()==24){
	//muf2=sqr(f[2].Mass()+ampl->Leg(3)->Mom().PPerp());
	muf2=sqr(Flavour(kf_t).Mass());
      }
      else {
	//	muf2=sqr(f[3].Mass()+ampl->Leg(2)->Mom().PPerp());
	muf2=sqr(Flavour(kf_t).Mass());
      }
    }
  }
  else if ((f[2].Kfcode()==6) + (f[3].Kfcode()==6)==1){
    Flavour top(NULL), jet(NULL);
    msg_Out()<<"single top found\n";
    for (int i=2; i<4; i++){
      if (f[i].Kfcode()==6) top=f[i];
      else jet=f[i];
    }
    if (jet.Kfcode()==23) muf2=sqr(175);
    else if (f[0].IsQuark() && f[1].IsQuark()) {
      if (jet.Kfcode()==5){
	muf2=sqr(Flavour(kf_t).Mass());
      }
      else{
	muf2=sqr(Flavour(kf_t).Mass());
      }
    }
    else if ((f[0].IsGluon() && f[1].IsQuark()) || (f[1].IsGluon() && f[0].IsQuark())) {
       muf2=sqr(Flavour(kf_t).Mass()); 
    } 
    else {
      THROW(fatal_error,"Invalid call");
    }
    //}
  }
  else if ((f[2].Kfcode()==23)+(f[3].Kfcode()==23)==1) muf2=sqr(175);
  else if ((f[0].Kfcode()==6) + (f[1].Kfcode()==6)==1){
    muf2=sqr(175);
  }
  else if (f[0].IsQuark() && f[1].IsQuark() && f[2].IsQuark() && f[3].IsQuark()){
    muf2=sqr(175);	
  }
  else if ((f[0].IsQuark() && f[1].IsGluon()) || (f[0].IsGluon() && f[1].IsQuark()))
    muf2=sqr(175);
  else if (f[2].IsGluon() && f[3].IsGluon() && f[0].IsQuark() && f[1].IsQuark())
    muf2=sqr(175);
  else if (f[0].IsGluon() && f[1].IsGluon() && f[2].IsQuark() && f[3].IsQuark())
    muf2=sqr(175);
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

DECLARE_ND_GETTER(TZj_Core_Scale,"TZj",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,TZj_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new TZj_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    TZj_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"tVj~ core scale"; 
}
