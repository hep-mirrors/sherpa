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
using PHASIC::Process_Info;

Single_Process::Single_Process() :
  p_born_me2(NULL), p_virtual_me2(NULL), m_nlotype(nlo_type::lo), m_localFS(false),
  m_stashed(false)
{ }

Single_Process::~Single_Process()
{
  if (p_born_me2)    delete p_born_me2;
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
  if (/*m_pinfo.m_special!="MPI_Process" && */ MODEL::s_model->Name()!="SM") {
    DEBUG_INFO("Requested BSM, Internal can't cope, it's too dumb...");
    return false;
  }

  m_nlotype=m_pinfo.m_fi.NLOType();

  if (m_nlotype==nlo_type::loop) {
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
      m_maxcpl[0]=m_mincpl[0]=p_born_me2->OrderQCD();
      m_maxcpl[1]=m_mincpl[1]=p_born_me2->OrderEW();
      p_born_me2->FillCombinations(m_ccombs,m_cfls);
      m_sprimemin = p_born_me2->SPrimeMin()>0.?p_born_me2->SPrimeMin():-1.;
      m_sprimemax = p_born_me2->SPrimeMax()>0.?p_born_me2->SPrimeMax():-1.;
      if (m_flavs[2]==Flavour(kf_instanton)) StashOriginalFlavours(); 
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

void Single_Process::StashOriginalFlavours() {
  m_flavs.reserve(100);
  m_stashed = true;
  for (size_t i=0;i<m_flavs.size();++i) m_original_flavs.push_back(m_flavs[i]);
  m_original_nout = m_nout;
}

void Single_Process::RestoreOriginalFlavours() {
  m_flavs.clear();
  for (size_t i=0;i<m_original_flavs.size();++i) m_flavs.push_back(m_original_flavs[i]);
  p_int->SetMomenta(m_original_momenta);
  m_nout = m_original_nout;
}

void Single_Process::OverwriteOriginalWithLocalFlavoursAndMomenta() {
  m_flavs.clear();
  for (size_t i=0;i<p_born_me2->Flavours().size();i++) {
    m_flavs.push_back(p_born_me2->Flavours()[i]);
  }
  m_original_momenta = p_int->Momenta();
  p_int->SetMomenta(p_born_me2->Momenta());
  m_nout = m_flavs.size()-m_nin;
}

double Single_Process::Partonic(const ATOOLS::Vec4D_Vector& momenta,
                                Variations_Mode varmode,
                                int mode)
{
  if (mode==1) return m_mewgtinfo.m_B=m_lastbxs=m_lastxs;
  if (m_nlotype==nlo_type::lo && !Selector()->Result())
    return m_mewgtinfo.m_B=m_lastbxs=m_lastxs=0.0;
  if (!p_born_me2->FillFinalState(momenta)) return 0.;
  m_localFS = true;
  if (m_stashed) OverwriteOriginalWithLocalFlavoursAndMomenta();
  p_scale->CalculateScale(momenta);
  if (m_stashed) RestoreOriginalFlavours();
  m_localFS = false;
  if (p_born_me2) m_mewgtinfo.m_B=m_lastbxs=m_lastxs=(*p_born_me2)(momenta)*KFactor();
  else if (p_virtual_me2) {
    p_virtual_me2->SetRenScale(p_scale->Scale(stp::ren));
    p_virtual_me2->Calc(momenta);
    m_lastbxs=p_virtual_me2->ME_Born();
    m_mewgtinfo.m_VI=m_lastxs=p_virtual_me2->ME_Finite()*KFactor();
  }
  return m_lastxs;
}

bool EXTRAXS::Single_Process::FillIntegrator
(PHASIC::Phase_Space_Handler *const psh)
{
  PHASIC::Multi_Channel *mc(psh->FSRIntegrator());
  mc->DropAllChannels();
  if (m_nin==2 && m_nout==1 && m_flavs[2]==Flavour(kf_instanton)) {
    mc->Add(new PHASIC::NoChannel(m_nin,m_nout,(Flavour*)&Flavours().front()));
    return false;
  }
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
  if (m_ccombs.size()) {
    std::set<std::pair<size_t,size_t> >::const_iterator 
      cit(m_ccombs.find(std::pair<size_t,size_t>(idi,idj)));
    return cit!=m_ccombs.end();
  }
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

const Flavour_Vector &Single_Process::CombinedFlavour(const size_t &idij)
{
  if (m_cfls.size()) {
    std::map<size_t,ATOOLS::Flavour_Vector>::const_iterator fit(m_cfls.find(idij));
    if (fit==m_cfls.end()) THROW(fatal_error,"Invalid request");
    return fit->second;
  }
  if (GetME()) return GetME()->CombinedFlavour(idij);
  static Flavour_Vector fls(1,kf_none);
  return fls;
}

bool Single_Process::FillFinalState(const ATOOLS::Vec4D_Vector &p) {
  return true;
}

std::vector<std::vector<int> > * Single_Process::Colours() const {
  return m_localFS?(&p_born_me2->Colours()):NULL;
}

const bool Single_Process::HasInternalScale() const {
  return p_born_me2?p_born_me2->HasInternalScale():false;
}

const double Single_Process::InternalScale() const {
  return p_born_me2?p_born_me2->InternalScale():-1.;
}
