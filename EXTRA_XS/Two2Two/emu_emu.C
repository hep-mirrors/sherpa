#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"


using namespace EXTRAXS;
using namespace ATOOLS;
using namespace MODEL;
using namespace PHASIC;
using namespace std;


namespace EXTRAXS {

  class emu_emu : public ME2_Base {
  private:

  public:

    emu_emu(const External_ME_Args& args);

    double operator()(const ATOOLS::Vec4D_Vector& mom);

    double m_alpha, m_alphas;
    double MM2, ME2;
    Flavour_Vector m_flvs;
    bool m_ss; // same sign
  };

  emu_emu::emu_emu(const External_ME_Args& args)
    : ME2_Base(args)
  {
    m_oew=0;
    m_oqcd=0;
    m_alpha  = MODEL::s_model->ScalarConstant("alpha_QED");
    m_alphas  = MODEL::s_model->ScalarConstant("strong_cpl");
    m_flvs = args.Flavours();
    if(m_flvs[0].Charge()==m_flvs[1].Charge()) m_ss=true;
    else m_ss = false;
    ME2 = pow(Flavour(kf_e).Mass(),2);
    MM2 = pow(Flavour(kf_mu).Mass(),2);
  }

  double emu_emu::operator()(const ATOOLS::Vec4D_Vector& momenta)
  {
    double s((momenta[0]+momenta[1]).Abs2());
    double t((momenta[0]-momenta[2]).Abs2());
    double u((momenta[0]-momenta[3]).Abs2());
    double amp = (8.*M_PI*M_PI*(8.*ME2*MM2 + pow(ME2 + MM2 - s,2) 
                  + 2.*MM2*(-2*ME2 + t) + 2*ME2*(-2*MM2 + t) 
                  + pow(ME2 + MM2 - u,2))/t/t)/16.;
    return  2.*pow(M_PI,3)*sqr((*aqed)(s))*amp; // flux = 0.5
  }
}

DECLARE_TREEME2_GETTER(EXTRAXS::emu_emu,"emu_emu")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::emu_emu>::
operator()(const External_ME_Args &args) const
{
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=4) return NULL;
  if ( ( fl[0].Kfcode() == Flavour(kf_e) && fl[1].Kfcode() == Flavour(kf_mu)  ) || 
       ( fl[0].Kfcode() == Flavour(kf_mu) && fl[1].Kfcode() == Flavour(kf_e)  )
    )
  { 
    return new emu_emu(args);
  }
  return NULL;
}
