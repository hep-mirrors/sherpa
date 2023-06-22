#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/UFO/UFO_Model.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "NEUTRINOS++/Current_Library/Current_ME.H"
#include "NEUTRINOS++/Current_Library/Scatter_Current_Base.H"
#include "NEUTRINOS++/Current_Library/Lepton_Lepton.H"
#include "NEUTRINOS++/Current_Library/Nucleon_Nucleon.H"
#include "NEUTRINOS++/Current_Library/Nucleon_Baryon_FFS.H"
#include "NEUTRINOS++/Current_Library/Nucleon_Baryon.H"

using namespace NEUTRINOS;
using namespace EXTRAXS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

namespace NEUTRINOS {

  class eP_eP : public ME2_Base {
  private:
    double m_mY2, m_mY2GY2;
    std::vector<int> m_indices;
    Current_ME           * p_ME;
    Scatter_Current_Base * p_JLL, * p_JNN;
  public:

    eP_eP(const External_ME_Args& args): ME2_Base(args)
    {
      m_indices.resize(2);
      m_indices[0] = 0;
      m_indices[1] = 2;
      p_JLL = new Lepton_Lepton(m_flavs,m_indices,"ee");
      m_indices[0] = 1;
      m_indices[1] = 3;
      p_JNN = new Nucleon_Nucleon(m_flavs,m_indices,"pp");
      m_indices.resize(4);
      m_indices[0] = 0;
      m_indices[1] = 2;
      m_indices[2] = 1;
      m_indices[3] = 3;
      p_ME  = new Current_ME(m_flavs, m_indices, "Current_ME");
      p_ME->SetCurrent1(p_JLL);
      p_ME->SetCurrent2(p_JNN);
      m_oew=2;
      m_oqcd=0;
    }
  
    double operator()(const ATOOLS::Vec4D_Vector& momenta)
    {
      return (*p_ME)(momenta);
    }

  };// end of class eP_eP
  
}

DECLARE_TREEME2_GETTER(eP_eP,"eP_eP")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,External_ME_Args,eP_eP>::
operator()(const External_ME_Args &args) const
{
  return NULL;
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  const Flavour_Vector fl=args.Flavours();
  if (fl.size()!=4) return NULL;
  if (fl[0]==Flavour(kf_e) && fl[2]==fl[0] &&
      fl[1].Kfcode()==kf_p_plus && fl[3].Kfcode()==kf_p_plus)
    return new eP_eP(args);
  return NULL;
}
