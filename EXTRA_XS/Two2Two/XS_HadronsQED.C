#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/UFO/UFO_Model.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


namespace EXTRAXS {

  class yy_bobo : public ME2_Base {
  private:
    double m_mass2, m_pref;
    int    m_charge;
  public:
    yy_bobo(const External_ME_Args& args);
    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };

  yy_bobo::yy_bobo(const External_ME_Args& args)
    : ME2_Base(args)
  {
    m_sintt = 1;
    m_oew   = 0;
    m_oqcd  = 0;
    const Flavour_Vector fl = args.Flavours();
    double charge = fl[2].Charge();
    if (charge == 0)
      THROW(fatal_error, "Cannot initialise yy -> BB for neutral bosons B.");
    m_mass2 = sqr(fl[2].HadMass());
    m_pref  = M_PI * sqr(MODEL::s_model->ScalarConstant("alpha_QED") *
			 sqr(charge));
    // doSelfTest(); // TODO trigger calculation based on message level
  }

  double yy_bobo::operator()(const ATOOLS::Vec4D_Vector& momenta)
  {
    // This matrix element to simulate the production of charged
    // boson pairs in EPA with scalar QED: no form factors etc.
    double shat((momenta[0]+momenta[1]).Abs2());
    if (4.*m_mass2>shat) return 0.;
    double mhat2   = m_mass2/shat;
    double vboson  = sqrt(1.-4.*mhat2);
    double logterm = log((sqrt(1./mhat2) + sqrt(1./mhat2-4.))/2.);
    // Calculation according to the paper
    double xs_perp = ( vboson * (1.+2.*mhat2) -
		       8.*mhat2 * (1. - mhat2) * logterm);
    double xs_para = ( vboson * (1.+6.*mhat2) -
		       8.*mhat2 * (1. - 3.*mhat2) * logterm);
    return m_pref/shat * (xs_perp + xs_para);
  }
}

DECLARE_TREEME2_GETTER(EXTRAXS::yy_bobo,"yy_bobo")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::yy_bobo>::
operator()(const External_ME_Args &args) const
{
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=4) return NULL;
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  if (fl[0]==Flavour(kf_photon) && fl[1]==Flavour(kf_photon) &&
      fl[2]==fl[3].Bar() && fl[2].IntSpin()==0 &&
      fl[2].Charge()!=0) {
    //    if (args.m_orders[0]==0 && args.m_orders[1]==2)
    //return new yy_bobo(args);
  }
  return NULL;
}
