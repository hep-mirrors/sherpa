#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "ATOOLS/Phys/FormFactor.H"


using namespace EXTRAXS;
using namespace ATOOLS;
using namespace MODEL;
using namespace PHASIC;
using namespace std;


namespace EXTRAXS {

  class ee_PiPi : public ME2_Base {
  private:

  public:

    ee_PiPi(const External_ME_Args& args,const finalstate::code &fs);

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    double Coulomb();
    double m_alpha, m_alphas, m_s;
    Flavour m_flv;
    std::unique_ptr<FormFactor> p_formfactor;
    finalstate::code m_fs;
  };

  ee_PiPi::ee_PiPi(const External_ME_Args& args, const finalstate::code &fs)
    : ME2_Base(args), m_fs(fs)
  {
    m_oew=0;
    m_oqcd=0;
    m_alpha = MODEL::s_model->ScalarConstant("alpha_QED");
    m_alphas = MODEL::s_model->ScalarConstant("strong_cpl");
    switch(m_fs){
      case finalstate::pion:
        m_flv =  Flavour(kf_pi_plus);
        break;
      case finalstate::kplus:
        m_flv = Flavour(kf_K_plus);
        break;
      case finalstate::off:
        m_flv = Flavour(args.m_outflavs[0]);
      default:
        THROW(fatal_error, "Unknown Final state");
        break;
    }
    p_formfactor = std::unique_ptr<FormFactor>(new FormFactor());
  }

  double ee_PiPi::Coulomb(){
    double m2 = 4.*m_flv.Mass()*m_flv.Mass();
    double v = 2.*sqrt((m_s-4.*m2)/m_s);
    v /= 1. + (m_s-m2)/m_s;
    double z = 2.*M_PI*(*aqed)(m_s)/v;
    return z/(1-exp(-z));
  }


  double ee_PiPi::operator()(const ATOOLS::Vec4D_Vector& momenta)
  {
    // Eq 90 in 0912.0749
    double m1 = momenta[0].Mass();
    m_s = (momenta[2]+momenta[3]).Abs2();
    double t((momenta[0]-momenta[2]).Abs2());
    double massFS = m_flv.Mass();
    double amp = -t*(m_s+t) - m1*m1*(2*massFS*massFS-m_s-2*t) -sqr(m1*m1) - sqr(massFS*massFS) + 2*massFS*massFS*t;
    switch(m_fs){
      case finalstate::pion:
        amp *= 1;
        break;
      case finalstate::kplus:
        amp *= Coulomb();
        break;
      case finalstate::off:
        break;
    }
     // = -t*(s+t) - m1*m1*(2*massFS*massFS-s-2*t) -sqr(m1*m1) - sqr(massFS*massFS) + 2*massFS*massFS*t;
    return  32*M_PI*M_PI*sqr((*aqed)(m_s))*amp/m_s/m_s;//*p_formfactor->Eval(s); // flux = 0.5
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

    return new ee_PiPi(args, finalstate::pion);
  }
  else if (fl[0]==Flavour(kf_e) && fl[1]==fl[0].Bar() &&
      (fl[2].Kfcode()==kf_K_plus || fl[2].Kfcode()==-kf_K_plus) && fl[3]==fl[2].Bar()){
      return new ee_PiPi(args, finalstate::kplus);
  }
  return NULL;
}
