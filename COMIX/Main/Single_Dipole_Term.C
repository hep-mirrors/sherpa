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
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
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
  p_proc(rs), p_bg(rs->GetAmplitude()), p_sub(sub), p_msub(msub)
{
  p_gen=rs->Generator();
  Process_Info info(rs->Info());
  info.Combine(sub->m_i,sub->m_j,msub->p_fl[sub->m_ijt]);
  info.m_fi.m_nlotype&=~nlo_type::real;
  Init(info,rs->Integrator()->Beam(),rs->Integrator()->ISR(),1);
  p_rsint=rs->Integrator();
  m_name+="_RS"+ToString(sub->m_i)+"_"
    +ToString(sub->m_j)+"_"+ToString(sub->m_k);
  m_maxcpl=rs->MaxOrders();
  m_mincpl=rs->MinOrders();
  for (size_t i(0);i<m_maxcpl.size();++i) {
    m_maxcpl[i]=m_maxcpl[i]-info.m_fi.m_nlocpl[i];
    m_mincpl[i]=m_mincpl[i]-info.m_fi.m_nlocpl[i];
  }
  p_bg->FillCombinations(m_ccombs,m_cflavs,&m_brs,p_sub);
  if (rs->IsMapped())
    for (CFlavVector_Map::iterator fit(m_cflavs.begin());
	 fit!=m_cflavs.end();++fit)
      for (size_t i(0);i<fit->second.size();++i)
	fit->second[i]=rs->ReMap(fit->second[i],m_brs[fit->first]);
}

Single_Dipole_Term::~Single_Dipole_Term()
{
  p_scale=NULL;
}

Weights_Map COMIX::Single_Dipole_Term::Differential(
    const Cluster_Amplitude &ampl,
    Variations_Mode varmode,
    int mode)
{
  DEBUG_FUNC(Name());
  m_zero=false;
  p_rsint->ColorIntegrator()->SetPoint(&ampl);
  return PHASIC::Process_Base::Differential(ampl,varmode,mode);
}

double COMIX::Single_Dipole_Term::Partonic(const Vec4D_Vector &p, int mode)
{
  Single_Dipole_Term *sp(this);
  if (mode==1) return m_lastxs;
  if (m_zero || !Selector()->Result()) return m_lastxs;
  if (!p_bg->RSTrigger(Selector(),m_mcmode))
    return m_lastxs=0.0;
  sp->p_scale->CalculateScale(p);
  if (m_mcmode==1) p_rsint->ColorIntegrator()->GeneratePoint();
  m_w=p_bg->KT2Trigger(p_sub,m_mcmode);
  if (m_w) sp->p_bg->Differential(p_sub);
  m_lastxs=-p_sub->m_me;
  m_w*=p_rsint->ColorIntegrator()->GlobalWeight();
  if (p_rsint->HelicityIntegrator()!=NULL) 
    m_w*=p_rsint->HelicityIntegrator()->Weight();
  m_w*=sp->KFactor(2);
  return m_lastxs*=m_w;
}

bool Single_Dipole_Term::Trigger(const ATOOLS::Vec4D_Vector &p)
{
  Selector()->SetResult(1);
  return p_bg->SetMomenta(p);
}

bool Single_Dipole_Term::GeneratePoint()
{
  return false;
}

bool Single_Dipole_Term::Combinable(const size_t &idi,const size_t &idj)
{
  Combination_Set::const_iterator 
    cit(m_ccombs.find(std::pair<size_t,size_t>(idi,idj)));
  return cit!=m_ccombs.end();
}

const Flavour_Vector &Single_Dipole_Term::CombinedFlavour(const size_t &idij)
{
  CFlavVector_Map::const_iterator fit(m_cflavs.find(idij));
  if (fit==m_cflavs.end()) THROW(fatal_error,"Invalid request");
  return fit->second;
}

NLO_subevtlist *Single_Dipole_Term::GetRSSubevtList()
{
  return p_rsint->Process()->GetSubevtList();
}
