#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "ATOOLS/Phys/Pion_FormFactor.H"


using namespace EXTRAXS;
using namespace ATOOLS;
using namespace MODEL;
using namespace PHASIC;
using namespace std;


namespace EXTRAXS {

  class ee_PiPi : public ME2_Base {
  private:

  public:

    ee_PiPi(const External_ME_Args& args);

    double operator()(const ATOOLS::Vec4D_Vector& mom);

    double m_alpha, m_alphas;
    Flavour m_flv;
    std::unique_ptr<Pion_FormFactor> p_formfactor;
  };

  ee_PiPi::ee_PiPi(const External_ME_Args& args)
    : ME2_Base(args)
  {
    m_oew=0;
    m_oqcd=0;
    m_alpha  = MODEL::s_model->ScalarConstant("alpha_QED");
    m_alphas  = MODEL::s_model->ScalarConstant("strong_cpl");

    m_flv =  Flavour(kf_pi_plus);
    p_formfactor = std::unique_ptr<Pion_FormFactor>(new Pion_FormFactor());
  }

  double ee_PiPi::operator()(const ATOOLS::Vec4D_Vector& momenta)
  {
    // Eq 90 in 0912.0749
    double m1 = momenta[0].Mass();
    double s((momenta[2]+momenta[3]).Abs2());
    double t((momenta[0]-momenta[2]).Abs2());
    double masspi = m_flv.Mass();
    double amp = -t*(s+t) - m1*m1*(2*masspi*masspi-s-2*t) -sqr(m1*m1) - sqr(masspi*masspi) + 2*masspi*masspi*t;
    return  32*M_PI*M_PI*sqr((*aqed)(s))*amp/s/s;//*p_formfactor->Eval(s); // flux = 0.5
  }
}

DECLARE_TREEME2_GETTER(EXTRAXS::ee_PiPi,"ee_PiPi")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::ee_PiPi>::
operator()(const External_ME_Args &args) const
{
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=4) return NULL;
  if (fl[0]==Flavour(kf_e) && fl[1]==fl[0].Bar() &&
      (fl[2].Kfcode()==kf_pi_plus || fl[2].Kfcode()==-kf_pi_plus) && fl[3]==fl[2].Bar())
  {
    return new ee_PiPi(args);
  }
  return NULL;
}
