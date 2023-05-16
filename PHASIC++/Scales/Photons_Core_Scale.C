#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  class Photons_Core_Scale: public Core_Scale_Setter {
  public:

    Photons_Core_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::Cluster_Param Calculate(ATOOLS::Cluster_Amplitude *const ampl);

  };// end of class Scale_Setter_Base

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

PDF::Cluster_Param Photons_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
{
  std::vector<Vec4D> photon_moms;
  for (auto& leg: ampl->Legs()) {
    if (leg->Flav().Kfcode()==kf_photon) photon_moms.push_back(leg->Mom());
  }

  double muf2;
  if (photon_moms.size()==0) {
    double s(2.0*ampl->Leg(0)->Mom()*ampl->Leg(1)->Mom());
    double t(2.0*ampl->Leg(0)->Mom()*ampl->Leg(2)->Mom());
    double u(2.0*ampl->Leg(0)->Mom()*ampl->Leg(3)->Mom());
    muf2 = -1.0/(1.0/s+1.0/t+1.0/u)/sqrt(2.0);
  }
  else if (photon_moms.size()==1) {
    muf2 = photon_moms[0].PPerp2();
  }
  else {
    Vec4D photonsmom;
    for (auto& mom: photon_moms) photonsmom += mom;
    muf2 = photonsmom.Abs2();
  }

  double mur2(muf2), q2(muf2);
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(muf2)<<"\n"
		 <<"  \\mu_r = "<<sqrt(mur2)<<"\n"
		 <<"  \\mu_q = "<<sqrt(q2)<<"\n";
  msg_Debugging()<<"}\n";
  return PDF::Cluster_Param(NULL,q2,muf2,mur2,-1);
}

DECLARE_ND_GETTER(Photons_Core_Scale,"Photons",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,Photons_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new Photons_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
		    Photons_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Photons core scale"; 
}
