#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MyComplex.H"
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
    double s((momenta[2]+momenta[3]).Abs2());
    double theta = momenta[1].Theta(momenta[2]);
    double masspi = m_flv.Mass();
    double beta = sqrt(1.-4*sqr(masspi)/s);
    double betaE = sqrt(1.-4*sqr(momenta[0].Mass())/s);
    double costh2 = sqr(cos(theta));
    return sqr(4*M_PI*(*aqed)(s)) / 2. * (1.-sqr(betaE)*costh2) / sqr(betaE) * pow(beta,3) * p_formfactor->Eval(s);
    // return 4*sqr((*aqed)(s)) * (1.-sqr(betaE)*costh2)/ s/ betaE * pow(beta,3) * p_formfactor->Eval(s);
    // return 2./3*32 * M_PI * s * sqr((*aqed)(s)) / 8. / s * sqr(betaE) * sinth2 * pow(beta,3) * p_formfactor->Eval(s);
    // return 32 * s * sqr((*aqed)(s)) / 4. / s * sinth2 * pow(beta,3) * p_formfactor->Eval(s);
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
