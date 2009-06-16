#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

#include "EXTRA_XS/Main/ME2_Base.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


namespace EXTRAXS {

  class ee_Y4S_BB : public ME2_Base {
  public:

    ee_Y4S_BB(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };

  ee_Y4S_BB::ee_Y4S_BB(const Process_Info& pi, const Flavour_Vector& fl)
    : ME2_Base(pi, fl)
  {
    m_sintt=1;
    m_oew=0;
    m_oqcd=0;
  }
  
  double ee_Y4S_BB::operator()(const ATOOLS::Vec4D_Vector& momenta)
  {
    // This is a dummy matrix element to allow simulating ee -> Y(4S) -> B Bbar
    // events.
    // Nobody cares about the production process of B Bbar, so distribute it
    // flat according to phase space, only with a factor to get the
    // cross section cited in the BaBar physics book, page 75 (3.39nb)
    return 1.29;
  }
}

DECLARE_ME2_GETTER(ee_Y4S_BB_Getter,"ee_Y4S_BB")
ME2_Base *ee_Y4S_BB_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_fi.NLOType()!=nlo_type::lo && pi.m_fi.NLOType()!=nlo_type::born)
    return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0]==Flavour(kf_e) && fl[1]==fl[0].Bar() &&
      (fl[2].Kfcode()==kf_B || fl[2].Kfcode()==kf_B_plus) && fl[3]==fl[2].Bar())
  {
    if (pi.m_fi.m_ps[0].m_fl.Kfcode()==kf_Upsilon_4S) {
      return new ee_Y4S_BB(pi, fl);
    }
  }
  return NULL;
}
