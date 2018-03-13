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

  class ee_Y4S_BB : public ME2_Base {
  private:
    double m_mY2, m_mY2GY2;
  public:

    ee_Y4S_BB(const External_ME_Args& args);

    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };

  ee_Y4S_BB::ee_Y4S_BB(const External_ME_Args& args)
    : ME2_Base(args)
  {
    m_sintt=1;
    m_oew=0;
    m_oqcd=0;

    Flavour flY(kf_Upsilon_4S);
    m_mY2    = sqr(flY.HadMass());
    m_mY2GY2 = m_mY2*sqr(flY.Width());

  }
  
  double ee_Y4S_BB::operator()(const ATOOLS::Vec4D_Vector& momenta)
  {
    // This is a dummy matrix element to allow simulating ee -> Y(4S) -> B Bbar
    // events.
    // Nobody cares about the production process of B Bbar, so distribute it
    // flat according to phase space, only with a factor to get the
    // cross section cited in the BaBar physics book, page 75 (3.39nb)
//     return 1.29;
    // better approximation with propagator
    // still norm to xsec
    double p2((momenta[0]+momenta[1]).Abs2());
    return 0.007955/(sqr(p2-m_mY2)+m_mY2GY2);
  }
}

DECLARE_TREEME2_GETTER(ee_Y4S_BB,"ee_Y4S_BB")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,External_ME_Args,ee_Y4S_BB>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=4) return NULL;
  if (fl[0]==Flavour(kf_e) && fl[1]==fl[0].Bar() &&
      (fl[2].Kfcode()==kf_B || fl[2].Kfcode()==kf_B_plus) && fl[3]==fl[2].Bar())
  {
    return new ee_Y4S_BB(args);
  }
  return NULL;
}
