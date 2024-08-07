#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Scale_Setter_Base
#define PARAMETER_TYPE PHASIC::Scale_Setter_Arguments
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

Scale_Setter_Base::Scale_Setter_Base
(const Scale_Setter_Arguments &args):
  p_proc(args.p_proc),
  p_model(args.p_model), p_cpls(args.p_cpls), p_subs(NULL),
  m_scale(stp::size), m_coupling(args.m_coupling),
  m_nin(args.m_nin), m_nout(args.m_nout),
  m_l1(0), m_l2(0)
{
  Settings& s = Settings::GetMainSettings();
  s["MEPS"]["CORE_SCALE"].SetDefault("Default");
  for (size_t i(0);i<stp::size;++i) m_scale[i]=sqr(rpa->gen.Ecms());
  if (p_proc) {
    m_nin=p_proc->NIn();
    m_nout=p_proc->NOut();
  }
  size_t nl(0);
  if (p_proc) {
    for (size_t i(m_nin);i<p_proc->Flavours().size();++i) {
      if (p_proc->Flavours()[i].IsLepton()) {
        nl++;
        if      (nl==1) m_l1=i;
        else if (nl==2) m_l2=i;
        else           {m_l1=m_l2=0; break;}
      }
    }
  }
  m_p.resize(m_nin+m_nout);
}

bool Scale_Setter_Base::Initialize()
{
  return true;
}

void Scale_Setter_Base::SetCouplings()
{
  if (p_proc==NULL || p_proc->Integrator()->Beam()==NULL) return;
  DEBUG_FUNC(p_proc->Name());
  if (p_cpls==NULL) THROW(fatal_error,"No coupling information");
  p_subs=p_proc->GetSubevtList();
  Data_Reader read(" ",",","#",":");
  std::vector<std::vector<std::string> > helpsvv;
  read.SetString(m_coupling);
  read.MatrixFromString(helpsvv,"");
  for (size_t i(0);i<helpsvv.size();++i) {
    if (helpsvv[i].size()!=2) {
      if (helpsvv[i].size()==1 && helpsvv[i][0]=="None") break;
      THROW(fatal_error,"Invalid tag "+m_coupling+".");
    }
    Coupling_Map::iterator cit(p_cpls->lower_bound(helpsvv[i][0]));
    Coupling_Map::iterator eit(p_cpls->upper_bound(helpsvv[i][0]));
    if (cit!=eit) {
      int idx(ToType<int>(helpsvv[i][1]));
      if (idx<0) continue;
      if (idx>=m_scale.size())
	THROW(fatal_error,"Index too large for "+helpsvv[i][0]+".");
      for (;cit!=eit;++cit) {
	msg_Debugging()<<*cit->second<<" -> "<<helpsvv[i][1]<<"\n";
	if (cit->second->Sub()==NULL) cit->second->SetScale(&m_scale[idx]);
	else {
	  cit->second->Sub()->m_mu2.resize(m_scale.size());
	  cit->second->SetScale(&cit->second->Sub()->m_mu2[idx]);
	}
      }
    }
    else {
      msg_Error()<<METHOD<<"("<<p_proc->Name()<<"): Valid tags are\n ";
      for (Coupling_Map::const_iterator cit(p_cpls->begin());
	   cit!=p_cpls->end();++cit) msg_Error()<<" "<<cit->first;
      msg_Error()<<"\n";
      THROW(fatal_error,"Invalid coupling tag "+helpsvv[i][0]+".");
    }
  }
  m_fac.resize(2,1.0);
}

Scale_Setter_Base::~Scale_Setter_Base()
{
}

void Scale_Setter_Base::ShowSyntax(const size_t i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   // available scale choices\n\n";
  Scale_Getter_Function::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}

double Scale_Setter_Base::HTM() const
{
  double htm(0.0);
  for (size_t i(m_nin);i<m_p.size();++i) htm+=m_p[i].MPerp();
  return htm;
}

double Scale_Setter_Base::PTM() const
{
  //product of transverse masses of all massive particles
  double ptm(1.0);
  size_t n(0);
  for (size_t i(m_nin);i<m_p.size();++i){
      ATOOLS::Flavour flav =p_proc->Flavours()[i];
      if (flav.IsMassive()){
        ptm*=m_p[i].MPerp();
        n++;
     }
    }
  return pow(ptm,1./n);
}

double Scale_Setter_Base::HT() const
{
  double ht(0.0);
  for (size_t i(m_nin);i<m_p.size();++i) ht+=m_p[i].PPerp();
  return ht;
}

double Scale_Setter_Base::HTMprime() const
{
  if (m_l1==0 || m_l2==0) THROW(fatal_error,"Lepton indices not set.");
  double htmp((m_p[m_l1]+m_p[m_l2]).MPerp());
  for (size_t i(m_nin);i<m_p.size();++i)
    if (i!=m_l1 && i!=m_l2) htmp+=m_p[i].MPerp();
  return htmp;
}

double Scale_Setter_Base::HTprime() const
{
  if (m_l1==0 || m_l2==0) THROW(fatal_error,"Lepton indices not set.");
  double htp((m_p[m_l1]+m_p[m_l2]).MPerp());
  for (size_t i(m_nin);i<m_p.size();++i)
    if (i!=m_l1 && i!=m_l2) htp+=m_p[i].PPerp();
  return htp;
}

Vec4D Scale_Setter_Base::PSum() const
{
  Vec4D sum(0.0,0.0,0.0,0.0);
  for (size_t i(m_nin);i<m_p.size();++i) sum+=m_p[i];
  return sum;
}

double Scale_Setter_Base::hHT() const
{
  // hadronic H_T
  double htj(0.0);
  for (size_t i(m_nin);i<m_p.size();++i)
    if (p_proc->Flavours()[i].Strong())
      htj+=m_p[i].PPerp();
  return htj;
}

double Scale_Setter_Base::CalculateScale
(const ATOOLS::Vec4D_Vector &p,const size_t mode)
{
  DEBUG_FUNC((p_proc?p_proc->Name():""));
  if (!m_escale.empty()) {
    for (size_t i(0);i<m_escale.size();++i) m_scale[i]=m_escale[i];
    while (m_ampls.size()) {
      m_ampls.back()->Delete();
      m_ampls.pop_back();
    }
    if (p_subs) {
      for (size_t i(0);i<p_subs->size();++i) {
	NLO_subevt *sub((*p_subs)[i]);
	size_t ssz(Min(sub->m_mu2.size(),m_scale.size()));
	for (size_t j(0);j<ssz;++j) sub->m_mu2[j]=m_scale[j];
	if (sub->p_ampl) {
	  sub->p_ampl->Delete();
	  sub->p_ampl=NULL;
	}
      }
    }
    p_cpls->Calculate();
    return m_scale[stp::fac];
  }
  if (p_subs==NULL) {
    m_p.resize(p.size());
    for (size_t j(0);j<m_p.size();++j) m_p[j]=p[j];
    Calculate(p,mode);
  }
  else {
    msg_Debugging()<<"calculating scales for subevents"<<std::endl;
    for (int i(p_subs->size()-1);i>=0;--i) {
      NLO_subevt *sub((*p_subs)[i]);
      if (!sub->m_trig && !sub->IsReal()) {
        for (size_t j(0);j<sub->m_mu2.size();++j) sub->m_mu2[j]=-1.0;
	if (sub->p_ampl) {
	  sub->p_ampl->Delete();
	  sub->p_ampl=NULL;
	}
	continue;
      }
      m_p.resize(sub->m_n);
      for (size_t j(0);j<m_p.size();++j)
	m_p[j]=sub->p_mom[j][0] < 0.0 ?-sub->p_mom[j]:sub->p_mom[j];

      /* For EXTAMP::Dipole_Wrapper_Processes, the flavour config does
	 not match the flavour config of the corresponding subevent
	 (the former having real emission flavours, the latter born
	 flavours). Need to locally fix that here. */
      p_proc->SetCaller(static_cast<Process_Base*>(sub->p_proc));
      Flavour_Vector tmp = p_proc->Caller()->Flavours();
      p_proc->Caller()->SetFlavours(Flavour_Vector(sub->p_fl,&sub->p_fl[sub->m_n]));
      Calculate(Vec4D_Vector(m_p),mode);
      p_proc->Caller()->SetFlavours(tmp);

      if (i+1==p_subs->size()) m_escale=m_scale;
      size_t ssz(Min(sub->m_mu2.size(),m_scale.size()));
      for (size_t j(0);j<ssz;++j) sub->m_mu2[j]=m_scale[j];
      if (sub->p_ampl) {
	sub->p_ampl->Delete();
	sub->p_ampl=NULL;
      }
      if (m_ampls.size()) {
	sub->p_ampl=m_ampls.back();
	for (Cluster_Amplitude *ampl(sub->p_ampl);
	     ampl;ampl=ampl->Next()) ampl->SetProc(sub->p_proc);
	m_ampls.pop_back();
      }
    }
    m_scale=m_escale;
    m_escale.clear();
  }
  if (p_proc && p_proc->Integrator()->Beam()) p_cpls->Calculate();
  msg_Debugging()<<"\\mu_F = "<<sqrt(m_scale[stp::fac])<<std::endl;
  msg_Debugging()<<"\\mu_R = "<<sqrt(m_scale[stp::ren])<<std::endl;
  msg_Debugging()<<"\\mu_Q = "<<sqrt(m_scale[stp::res])<<std::endl;
  return m_scale[stp::fac];
}

bool Scale_Setter_Base::UpdateScale(const QCD_Variation_Params &var)
{
  return false;
}
