#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/UFO/UFO_Model.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


namespace EXTRAXS {

  class ee_PiPi : public ME2_Base {
  private:

  public:

    ee_PiPi(const External_ME_Args& args);

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    double PionFormFactor(const double &s);

    double m_alpha, m_alphas;
    Flavour m_flv;
  };

  ee_PiPi::ee_PiPi(const External_ME_Args& args)
    : ME2_Base(args)
  {
    m_sintt=1;
    m_oew=0;
    m_oqcd=0;
    m_alpha  = MODEL::s_model->ScalarConstant("alpha_QED");
    m_alphas  = MODEL::s_model->ScalarConstant("strong_cpl");

    m_flv =  Flavour(kf_pi_plus);

  }

  double ee_PiPi::PionFormFactor(const double &s){
    // Todo add form factor
    return 1;
  }
  
  double ee_PiPi::operator()(const ATOOLS::Vec4D_Vector& momenta)
  {
    // Eq 90 in 0912.0749
    double s((momenta[2]+momenta[3]).Abs2());
    double t((momenta[0]-momenta[2]).Abs2());
    double masspi = m_flv.Mass();
    double beta = sqrt(1.-4*sqr(masspi)/s);
    double theta = momenta[0].Theta(momenta[2]);
    // cos(theta) = 1.+2.*t/s
    double sinth2 = 1.-sqr(1.+2.*t/s);
    // double sinth2 = sqr(sin(theta));
    if(sinth2 < 0) msg_Error()<<METHOD<<" Sin^2 less than zero: "<<sinth2<<std::endl;
    return sqr(4*M_PI*m_alpha)*2*M_PI/8. / s  * sinth2 * pow(beta,3) * PionFormFactor(s)/m_alphas;
  }
}

DECLARE_TREEME2_GETTER(ee_PiPi,"ee_PiPi")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,External_ME_Args,ee_PiPi>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=4) return NULL;
  if (fl[0]==Flavour(kf_e) && fl[1]==fl[0].Bar() &&
      (fl[2].Kfcode()==kf_pi_plus || fl[2].Kfcode()==-kf_pi_plus) && fl[3]==fl[2].Bar())
  {
    return new ee_PiPi(args);
  }
  return NULL;
}
