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
#include "NEUTRINOS++/Current_Library/Nucleon_Baryon.H"

#include "NEUTRINOS++/Tools/Form_Factor_Parameter_Maps.H"

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

  //Check for GE, GM form factors
  cpl_info::code GE = cpl_info::GE;
  bool Bool_GMsGEs = ffs->ContainsFormFactorType(GE);

  //Check for fs, gs form factors
  cpl_info::code f1 = cpl_info::f1;
  bool Bool_fsgs = ffs->ContainsFormFactorType(f1);

  if (Bool_GMsGEs && Bool_fsgs) exit(1); //Can't mix declaration types TODO JW: Add exit statement
  if (!Bool_GMsGEs && !Bool_fsgs) exit(1); //Need some form factors! TODO JW: Add exit statement

  m_indices.resize(2);
  m_indices[0] = 2;
  m_indices[1] = 0;
  p_JLL  = new Lepton_Lepton(m_flavs,m_indices,"ee");
  m_indices.resize(2);
  m_indices[0] = 3;
  m_indices[1] = 1;
  if (Bool_GMsGEs) p_JNN  = new Nucleon_Nucleon(m_flavs,m_indices,"pp");
  if (Bool_fsgs) p_JNN  = new Nucleon_Baryon(m_flavs,m_indices,"pp");
  m_indices.resize(4);
  m_indices[0] = 2;
  m_indices[1] = 0;
  m_indices[2] = 3;
  m_indices[3] = 1;
  p_ME   = new Current_ME(m_flavs, m_indices, "Current_ME");
  msg_Info()<<METHOD<<": setting currents "<<p_JLL->Name()<<" + "<<p_JNN->Name()<<" into "<<p_ME->Name()<<".\n";
  p_ME->SetCurrent1(p_JLL);
  p_ME->SetCurrent2(p_JNN);
  m_oew  = 2;
  m_oqcd = 0;
}

XS_lepton_nucleon::~XS_lepton_nucleon() {
  if (p_ME)  delete p_ME;
  if (p_JLL) delete p_JLL;
  if (p_JNN) delete p_JNN;
}

double XS_lepton_nucleon::operator()(const ATOOLS::Vec4D_Vector& momenta) {
  bool printout = false;
  msg->SetPrecision(16);
  double result = (*p_ME)(momenta);

  if (printout) {
    double t = (momenta[0]-momenta[2]).Abs2(); 
    msg_Out()<<METHOD<<" |M|^2("<<sqrt(-t)<<" vs "<<momenta[2].PPerp()<<"), t = "<<t<<", "<<"result/t^2 = "<<(result/(t*t))<<".\n\n\n";

    for (size_t i=0;i<4;i++)msg_Out()<<"   "<<i<<"     "<<m_flavs[i]<<": " <<momenta[i]<<", "<<"m = "<<sqrt(Max(0.,momenta[i].Abs2()))<<".\n";
    msg_Out()<<"\n\n\n";

    msg_Out()<<"Joe Python Script"<<":\n";
    for (size_t i=0;i<4;i++) msg_Out()<<"p_"<<i+1<<" = " <<"Particle.fromFourMom("<<momenta[i][0]<<","<<momenta[i][1]<<","<<momenta[i][2]<<","<<momenta[i][3]<<", "<<"anti=False, flavour=None)"<<"\n";
    msg_Out()<<"Sherpa = "<<result<<"\n";
    msg_Out()<<"\n\n\n";
  }
  return result;
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
