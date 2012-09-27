#include "COMIX/Main/Single_Dipole_Term.H"

#include "COMIX/Main/Process_Group.H"
#include "PDF/Main/ISR_Handler.H"
#include "COMIX/Phasespace/PS_Generator.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/CXXFLAGS.H"

using namespace COMIX;
using namespace PHASIC;

Single_Dipole_Term::Single_Dipole_Term
(COMIX::Single_Process *const rs,
 NLO_subevt *const sub,NLO_subevt *const msub):
  p_bg(rs->GetAmplitude()), p_sub(sub), p_msub(msub)
{
  delete p_int;
  p_int=rs->Integrator();
  p_gen=rs->Generator();
  Init(rs->Info(),rs->Integrator()->Beam(),rs->Integrator()->ISR());
  m_name.erase(m_name.length()-3,1);
  m_name+="_RS"+ToString(sub->m_i)+"_"
    +ToString(sub->m_j)+"_"+ToString(sub->m_k);
}

Single_Dipole_Term::~Single_Dipole_Term()
{
  p_selector=NULL;
  p_scale=NULL;
  p_int=NULL;
}

double COMIX::Single_Dipole_Term::Differential
(const Cluster_Amplitude &ampl,int mode) 
{
  DEBUG_FUNC(Name());
  m_zero=false;
  p_int->ColorIntegrator()->SetPoint(&ampl);
  return PHASIC::Process_Base::Differential(ampl,mode);
}

double COMIX::Single_Dipole_Term::Partonic
(const Vec4D_Vector &p,const int mode) 
{
  Single_Dipole_Term *sp(this);
  if (mode==1) return m_lastxs;
  if (m_zero || !Selector()->Result()) return m_lastxs;
  for (size_t i(0);i<m_nin+m_nout;++i) {
    double psm(m_flavs[i].Mass());
    if (p[i][0]<psm) return m_lastxs;
  }
  sp->p_scale->CalculateScale(p);
  if (!p_bg->JetTrigger(Selector(),m_mcmode))
    return m_lastxs=0.0;
  if (m_mcmode==1) p_int->ColorIntegrator()->GeneratePoint();
  m_w=p_bg->KT2Trigger(p_sub,m_mcmode);
  if (m_w) sp->p_bg->Differential(p,p_sub);
  m_lastxs=-p_sub->m_me;
  m_w*=p_int->ColorIntegrator()->GlobalWeight();
  if (p_int->HelicityIntegrator()!=NULL) 
    m_w*=p_int->HelicityIntegrator()->Weight();
  m_w*=sp->KFactor();
  return m_lastxs*=m_w;
}

bool Single_Dipole_Term::Trigger(const ATOOLS::Vec4D_Vector &p)
{
  if (!Selector()->NoJetTrigger(p)) return false;
  return p_bg->SetMomenta(p);
}

bool Single_Dipole_Term::GeneratePoint()
{
}

bool Single_Dipole_Term::Combinable(const size_t &idi,const size_t &idj)
{
}

const Flavour_Vector &Single_Dipole_Term::CombinedFlavour(const size_t &idij)
{
}
