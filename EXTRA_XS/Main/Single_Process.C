#include "EXTRA_XS/Main/Single_Process.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/FSR_Channel.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "MODEL/Main/Model_Base.H"

#include "EXTRA_XS/Main/ME2_Base.H"
#include "EXTRA_XS/NLO/Virtual_ME2_Base.H"
#include "EXTRA_XS/Main/ME_Base.H"
#include "HELICITIES/Main/Spin_Structure.H"

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
  if (m_nin!=2 || (m_nout!=2 && m_nout!=3)) return false;
  
  // can't do resonant processes, with one exception: ee -> Y(4S) -> B Bbar
  if (m_pinfo.m_fi.m_ps.size()!=m_pinfo.m_fi.NExternal()) {
    if (m_pinfo.m_fi.m_ps[0].m_fl.Kfcode()!=kf_Upsilon_4S) {
      DEBUG_INFO("found decay process, which Internal can't handle.");
      return false;
    }
  }
  
  // can't do any BSM
  if (MODEL::s_model->Name()!="SM") {
    DEBUG_INFO("Requested BSM, Internal can't cope, it's too dumb...");
    return false;
  }

  m_nlotype=m_pinfo.m_fi.NLOType();
  
  if (m_nlotype==nlo_type::loop || m_nlotype==nlo_type::vsub) {
    DEBUG_INFO("searching loop process");
    p_virtual_me2=Virtual_ME2_Base::GetME2(m_pinfo);
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
    p_born_me2=ME2_Base::GetME2(m_pinfo);
    if (p_born_me2!=NULL) {
      DEBUG_INFO("found");
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

double Single_Process::Differential(const ATOOLS::Vec4D_Vector& momenta) 
{
  if (m_nlotype==nlo_type::lo && !Trigger()) return m_lastxs=m_last=0.0;
  
  p_scale->CalculateScale(p_int->PSHandler()->LabPoint());

  if (p_born_me2) {
    m_lastxs=(*p_born_me2)(momenta);
  }
  else if (p_virtual_me2) {
    p_virtual_me2->SetRenScale(p_scale->Scale(PHASIC::stp::ren));
    p_virtual_me2->Calc(momenta);
    m_lastxs=p_virtual_me2->Result().GetFinite();
  }

  if (p_int->ISR() && m_nin==2) {
    if (p_int->ISR()->On()) {
      p_int->ISR()->MtxLock();
      if (!p_int->ISR()->CalculateWeight(p_scale->Scale(PHASIC::stp::fac))) {
	p_int->ISR()->MtxUnLock();
	return m_last=m_lastlumi=0.0;
      }
      m_lastlumi=p_int->ISR()->Weight(&m_flavs.front()); 
      p_int->ISR()->MtxUnLock();
    }
    else m_lastlumi=1.;
  }
  else m_lastlumi=1.;
  return m_last=m_lastxs*m_lastlumi*KFactor();
}

double Single_Process::Differential2() 
{
  if (m_lastxs==0.0 || !Trigger()) return 0.0;
  if (p_int->ISR() && m_nin==2) {
    p_scale->CalculateScale2(p_int->PSHandler()->LabPoint());
    if (m_flavs[0]==m_flavs[1] || p_int->ISR()->On()==0)
      return 0.;
    p_int->ISR()->MtxLock();
    if (!p_int->ISR()->CalculateWeight2(p_scale->Scale(PHASIC::stp::fac))) {
      p_int->ISR()->MtxUnLock();
      return 0.0;
    }
    double tmp=m_lastxs*p_int->ISR()->Weight2(&m_flavs.front()); 
    p_int->ISR()->MtxUnLock();
    m_last+=tmp*=KFactor2();
    return tmp;
  }
  else {
    return 0.0;
  }
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
