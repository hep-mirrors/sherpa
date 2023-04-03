#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/UFO/UFO_Model.H"
#include "MODEL/Main/Model_Base.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "NEUTRINOS++/Current_Library/Current_ME.H"
#include "NEUTRINOS++/Current_Library/Lepton_Lepton.H"
#include "NEUTRINOS++/Current_Library/Nucleon_Nucleon.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace NEUTRINOS;
using namespace std;


/*
   In all the differential cross sections the factor 1/16 Pi is cancelled
   by the factor 4 Pi for each alpha. Hence one Pi remains in the game.
*/

namespace EXTRAXS {
  class XS_lepton_nucleon : public ME2_Base {  // == XS_ffbar_ee but not XS_ffbar_f'fbar' !
  private:
    NEUTRINOS::Current_ME           * p_ME;
    NEUTRINOS::Scatter_Current_Base * p_JLL, * p_JNN;
    std::vector<int> m_indices;
  public:
    XS_lepton_nucleon(const External_ME_Args& args);
    ~XS_lepton_nucleon();
    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };
}

XS_lepton_nucleon::XS_lepton_nucleon(const External_ME_Args& args)
  : ME2_Base(args), p_ME(NULL), p_JLL(NULL), p_JNN(NULL)
{
  DEBUG_INFO("now entered EXTRAXS::XS_lepton_nucleon ...");
  m_indices.resize(4);
  m_indices[0] = 0; m_indices[1] = 2; m_indices[2] = 1; m_indices[3] = 3;
  p_JLL = new Lepton_Lepton(m_flavs,m_indices,"ee");
  p_JNN = new Nucleon_Nucleon(m_flavs,m_indices,"pp");
  p_ME  = new Current_ME(m_flavs, m_indices, "Current_ME");
  p_ME->SetCurrent1(p_JLL);
  p_ME->SetCurrent2(p_JNN);
}

XS_lepton_nucleon::~XS_lepton_nucleon() {
  if (p_ME)  delete p_ME;
  if (p_JLL) delete p_JLL;
  if (p_JNN) delete p_JNN;
}

double XS_lepton_nucleon::operator()(const ATOOLS::Vec4D_Vector& momenta) {
  msg_Out()<<METHOD<<":\n";
  for (size_t i=0;i<4;i++) msg_Out()<<"   "<<m_flavs[i]<<": "<<m_momenta[i]<<"\n";
  double result = (*p_ME)(momenta);
  msg_Out()<<"|M|^2 = "<<result<<".\n";
  exit(1);
}

DECLARE_TREEME2_GETTER(XS_lepton_nucleon,"XS_lepton_nucleon")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,External_ME_Args,XS_lepton_nucleon>::
operator()(const External_ME_Args& args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl=args.Flavours();
  msg_Out()<<METHOD<<" for "<<fl.size()<<" flavours: "
	   <<fl[0]<<" "<<fl[1]<<" --> "<<fl[2]<<" "<<fl[3]<<" "
	   <<"tags: "<<fl[1].IsBaryon()<<" && "<<(fl[1].IntSpin()==1)<<" "
	   <<fl[3].IsBaryon()<<" && "<<(fl[3].IntSpin()==1)<<" "<<"\n";
  if (fl.size()!=4) return NULL;
  if (fl[0].IsLepton() && fl[2].IsLepton() &&
      fl[1].IsBaryon() && fl[1].IntSpin()==1 &&
      fl[3].IsBaryon() && fl[3].IntSpin()==1) {
    return new XS_lepton_nucleon(args);
  }
  return NULL;
}
