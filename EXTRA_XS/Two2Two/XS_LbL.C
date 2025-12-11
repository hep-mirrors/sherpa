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
  class yy_yy : public ME2_Base {
  private:
    int m_addRegge, m_addResonances, m_addMesons;
    Complex ReggeAmp(vector<int> & hels);
    Complex MesonAmp(vector<int> & hels);
    Complex ResonanceAmp(vector<int> & hels);
    Complex Amplitude(vector<int> & hels);
  public:
    yy_yy(const External_ME_Args& args);
    double operator()(const ATOOLS::Vec4D_Vector& mom);    
  };

  yy_yy::yy_yy(const External_ME_Args& args) {
  }
}

DECLARE_TREEME2_GETTER(EXTRAXS::yy_yy,"yy_yy")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
			      PHASIC::External_ME_Args,EXTRAXS::yy_yy>::
operator()(const External_ME_Args &args) const
{
  const Flavour_Vector fl = args.Flavours();
  if (fl.size()!=4) return NULL;
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  if (fl[0]==Flavour(kf_photon) && fl[1]==Flavour(kf_photon) &&
      fl[2]==Flavour(kf_photon) && fl[3]==Flavour(kf_photon) ) {
    return new yy_yy(args);
  }
  return NULL;
}
