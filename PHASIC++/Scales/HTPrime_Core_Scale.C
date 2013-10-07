#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace PHASIC;

namespace PHASIC {

  class HTPrime_Core_Scale: public Core_Scale_Setter {

  public:

    HTPrime_Core_Scale(const Core_Scale_Arguments &args);

    ~HTPrime_Core_Scale();

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl);

  };

  HTPrime_Core_Scale::HTPrime_Core_Scale
  (const Core_Scale_Arguments &args): Core_Scale_Setter(args)
  {
  }

  HTPrime_Core_Scale::~HTPrime_Core_Scale()
  {
  }

  PDF::CParam HTPrime_Core_Scale::Calculate(Cluster_Amplitude *const ampl)
  {
    size_t idx_l1(0), idx_l2(0);
    double ht_hat_prime(0.0);
    for (size_t i=2; i<ampl->Legs().size(); ++i) {
      if (ampl->Leg(i)->Flav().IsLepton()) {
        if (ampl->Leg(i)->Flav().IsAnti()) {
          if (idx_l2>0) THROW(fatal_error, "Too many anti-leptons.");
          idx_l2=i;
        }
        else {
          if (idx_l1>0) THROW(fatal_error, "Too many leptons.");
          idx_l1=i;
        }
      }
      else if (ampl->Leg(i)->Flav().IsQuark() || ampl->Leg(i)->Flav().IsGluon()) {
        ht_hat_prime += ampl->Leg(i)->Mom().PPerp();
      }
      else {
        THROW(fatal_error, "Encountered non-lepton/parton.");
      }
    }
    if (idx_l1==0 || idx_l2==0) {
      THROW(fatal_error, "Did not find two leptons.");
    }

    Vec4D wmom(ampl->Leg(idx_l1)->Mom()+ampl->Leg(idx_l2)->Mom());
    ht_hat_prime += wmom.MPerp();
    double mu2 = sqr(ht_hat_prime/2.0);
    double muq2 = wmom.Abs2();

    msg_Debugging()<<METHOD<<"(): Set {\n"
                   <<"  \\mu_f = "<<sqrt(mu2)<<"\n"
                   <<"  \\mu_r = "<<sqrt(mu2)<<"\n"
                   <<"  \\mu_q = "<<sqrt(muq2)<<"\n";
    msg_Debugging()<<"}\n";
    return PDF::CParam(mu2,muq2,0.0,mu2,-1);
  }

}

DECLARE_ND_GETTER(HTPrime_Core_Scale,"HTPrime",
                  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,HTPrime_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new HTPrime_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,
                    HTPrime_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"HTPrime core scale"; 
}
