#include "EXTRA_XS/Main/Single_Process.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/FSR_Channel.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "MODEL/Main/Model_Base.H"

#include "EXTRA_XS/Main/ME2_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "METOOLS/Main/Spin_Structure.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using PHASIC::nlo_type;
using PHASIC::Process_Info;

Single_Process::Single_Process() :
  p_born_me2(NULL), p_virtual_me2(NULL), m_nlotype(nlo_type::lo)
{
}

Single_Process::~Single_Process()
{
  if (p_born_me2) delete p_born_me2;
  if (p_virtual_me2) delete p_virtual_me2;
}

bool Single_Process::Initialize()
{
  DEBUG_FUNC(&m_pinfo);
  DEBUG_VAR(m_pinfo);
  MODEL::s_model->GetCouplings(m_cpls);
  if (m_nin!=2) return false;
  
  // can't do resonant processes, with one exception: ee -> Y(4S) -> B Bbar
  if (m_pinfo.m_fi.m_ps.size()!=m_pinfo.m_fi.NExternal()) {
    if (m_pinfo.m_fi.m_ps[0].m_fl.Kfcode()!=kf_Upsilon_4S) {
      DEBUG_INFO("found decay process, which Internal can't handle.");
      return false;
    }
  }
  
  // can't do any BSM
  if (m_pinfo.m_mpiprocess==false && MODEL::s_model->Name()!="SM") {
    DEBUG_INFO("Requested BSM, Internal can't cope, it's too dumb...");
    return false;
  }

  m_nlotype=m_pinfo.m_fi.NLOType();
  
  if (m_nlotype==nlo_type::loop || m_nlotype==nlo_type::vsub) {
    DEBUG_INFO("searching loop process");
    p_virtual_me2=PHASIC::Virtual_ME2_Base::GetME2(m_pinfo);
    if (p_virtual_me2!=NULL) {
      DEBUG_INFO("found");
      return true;
    }
    else {
      DEBUG_INFO("not found ...");
      return false;
    }
  }
  else if (m_nlotype==nlo_type::lo || m_nlotype==nlo_type::born ||
           m_nlotype==nlo_type::real || m_nlotype==nlo_type::rsub) {
    DEBUG_INFO("searching tree process");
    p_born_me2=dynamic_cast<ME2_Base*>
      (PHASIC::Tree_ME2_Base::GetME2(m_pinfo));
    if (p_born_me2!=NULL) {
      DEBUG_INFO("found");
      p_born_me2->SetCouplings(m_cpls);
      m_oqcd=p_born_me2->OrderQCD();
      m_oew=p_born_me2->OrderEW();
      return true;
    }
    else {
      DEBUG_INFO("not found ...");
      return false;
    }
  }
  else {
    DEBUG_INFO("don't know about processes of type "<<m_nlotype);
    return false;
  }
}

double Single_Process::Partonic(const ATOOLS::Vec4D_Vector& momenta,
				const int mode) 
{
  if (mode==1) return m_lastxs;
  if (m_nlotype==nlo_type::lo && !Selector()->Result()) return m_lastxs=0.0;
  
  p_scale->CalculateScale(momenta);
  if (p_born_me2) {
    m_lastxs=(*p_born_me2)(momenta);
  }
  else if (p_virtual_me2) {
    p_virtual_me2->SetRenScale(p_scale->Scale(stp::ren));
    p_virtual_me2->Calc(momenta);
    m_lastxs=p_virtual_me2->Result().GetFinite();
  }
  return m_lastxs*=KFactor();
}

bool EXTRAXS::Single_Process::FillIntegrator
(PHASIC::Phase_Space_Handler *const psh)
{
  PHASIC::Multi_Channel *mc(psh->FSRIntegrator());
  mc->DropAllChannels();
  size_t sintt(7);
  if (GetME()) sintt=GetME()->SIntType();
  if (sintt&1)
    mc->Add(new PHASIC::S1Channel(m_nin,m_nout,(Flavour*)&Flavours().front()));
  if (sintt&2)
    mc->Add(new PHASIC::T1Channel(m_nin,m_nout,(Flavour*)&Flavours().front()));
  if (sintt&4)
    mc->Add(new PHASIC::U1Channel(m_nin,m_nout,(Flavour*)&Flavours().front()));
  return false;
}

bool Single_Process::Combinable(const size_t &idi,const size_t &idj)
{
  size_t sintt(7);
  if (GetME()) sintt=GetME()->SIntType();
  if ((idi==1 && idj==2) || (idi==4 && idj==8)) {
    return sintt&1;
  }
  else if ((idi==1 && idj==4) || (idi==2 && idj==8)) {
    return sintt&2;
  }
  else if ((idi==1 && idj==8) || (idi==2 && idj==4)) {
    return sintt&4;
  }
  else {
    return false;
  }
}

const Flavour_Vector &Single_Process::
CombinedFlavour(const size_t &idij)
{
  if (GetME()) return GetME()->CombinedFlavour(idij);
  static Flavour_Vector fls(1,kf_none);
  return fls;
}
