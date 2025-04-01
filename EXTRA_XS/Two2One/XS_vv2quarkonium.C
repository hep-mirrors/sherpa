#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flow.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "MODEL/UFO/UFO_Model.H"

#include "EXTRA_XS/Main/ME2_Base.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;


namespace EXTRAXS {
  class XS_gg_quarkonium : public ME2_Base {
  private:
    double m_mass, m_mass2, m_width; 
  public:
    XS_gg_quarkonium(const External_ME_Args& args);

    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);

  };
}

DECLARE_TREEME2_GETTER(EXTRAXS::XS_gg_quarkonium,"1XS_gg_quarkonium")
Tree_ME2_Base *ATOOLS::
Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,EXTRAXS::XS_gg_quarkonium>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=3) return NULL;
  if (fl[0].IsGluon() && fl[1].IsGluon()) {
    return new XS_gg_quarkonium(args);
  }
  return NULL;
}
