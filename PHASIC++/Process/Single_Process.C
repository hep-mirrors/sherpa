#include "PHASIC++/Process/Single_Process.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

Single_Process::Single_Process(): m_zero(false)
{
}

Single_Process::~Single_Process()
{
  for (Coupling_Map::const_iterator
	 cit(m_cpls.begin());cit!=m_cpls.end();++cit)
    delete cit->second;
}

size_t Single_Process::Size() const
{
  return 1;
}

Process_Base *Single_Process::operator[](const size_t &i)
{
  if (i==0) return this;
  return NULL;
}

void Single_Process::DeSelect()
{
  p_selected=NULL;
}

bool Single_Process::SelectOne()
{
  p_selected=this;
  return true;
}

Weight_Info *Single_Process::OneEvent() 
{
  SelectOne();
  return p_int->PSHandler()->OneEvent(this);
}

Weight_Info *Single_Process::WeightedEvent(const int mode) 
{
  SelectOne();
  return p_int->PSHandler()->WeightedEvent(this,mode);
}

double Single_Process::KFactor() const
{
  return p_kfactor->KFactor();
}

double Single_Process::BeamISRWeight
(const double& Q2,const int mode) const
{
  double wgt(1.0);
  if (m_nin!=2) return wgt;
  if (p_int->ISR() && p_int->ISR()->On()) {
    p_int->ISR()->MtxLock();
    int stat(mode?p_int->ISR()->CalculateWeight2(Q2):
	     p_int->ISR()->CalculateWeight(Q2));
    if (!stat) {
      p_int->ISR()->MtxUnLock();
      return 0.0;
    }
    wgt*=mode?p_int->ISR()->Weight2(&m_flavs.front()):
      p_int->ISR()->Weight(&m_flavs.front()); 
    p_int->ISR()->MtxUnLock();
  }
  if (p_int->Beam() && p_int->Beam()->On()) {
    p_int->Beam()->MtxLock();
    p_int->Beam()->CalculateWeight(Q2);
    wgt*=p_int->Beam()->Weight();
    p_int->Beam()->MtxUnLock();
  }
  return wgt;
}

void Single_Process::BeamISRWeight
(NLO_subevtlist *const subs,const int mode) const
{
  double muf2(subs->m_muf2);
  if (m_nin==2 && p_int->ISR() && p_int->ISR()->On()) {
    size_t nscales(0);
    for (size_t i(0);i<subs->size();++i) {
      NLO_subevt *sub((*subs)[i]);
      if (!IsEqual(sub->m_facscale,muf2)) {
	if (mode==0) (*sub)*=BeamISRWeight(sub->m_facscale,mode);
	else sub->m_result+=sub->m_me*BeamISRWeight(sub->m_facscale,mode);
	++nscales;
      }
    }
    if (nscales<subs->size()) {
      double lumi(BeamISRWeight(muf2,mode));
      for (size_t i(0);i<subs->size();++i) 
	if (IsEqual((*subs)[i]->m_facscale,muf2)) {
	  if (mode==0) (*(*subs)[i])*=lumi;
	  else (*subs)[i]->m_result+=(*subs)[i]->m_me*lumi;
	}
    }
  }
}

double Single_Process::Differential(const Vec4D_Vector &p)
{
  m_wgtinfo.m_w0 = m_last=0.0;
  p_int->SetMomenta(p);
  if (GetSubevtList()==NULL) {
    if (m_zero) return 0.0;
    Scale_Setter_Base *scs(ScaleSetter());
    if (IsMapped()) {
      p_mapproc->Integrator()->SetMomenta(p);
      scs=p_mapproc->ScaleSetter();
    }
    if (Partonic(p)==0.0) return 0.0;
    if (m_wgtinfo.m_nx==0) m_wgtinfo.m_w0 = m_lastxs;
    m_wgtinfo.m_renscale = scs->Scale(stp::ren);
    if (m_lastxs==0.0) return m_last=0.0; 
    return m_last=m_lastxs*BeamISRWeight(scs->Scale(stp::fac),0);
  }
  Partonic(p);
  NLO_subevtlist *subs(GetSubevtList());
  BeamISRWeight(subs,0);
  for (size_t i=0;i<subs->size();++i) m_last+=(*subs)[i]->m_result;
  return m_last;
}

double Single_Process::Differential2()
{
  if (m_nin!=2 || p_int->ISR()==NULL ||
      p_int->ISR()->On()==0) return 0.0;
  if (m_flavs[0]==m_flavs[1]) return 0.0;
  if (GetSubevtList()==NULL) {
    if (m_lastxs==0.0) return 0.0;
    Scale_Setter_Base *scs((IsMapped()?p_mapproc:this)->ScaleSetter());
    double tmp=Partonic2();
    if (tmp!=0.0) m_last+=tmp*=BeamISRWeight(scs->Scale(stp::fac),1);
    return tmp;
  }
  Partonic2();
  NLO_subevtlist *subs(GetSubevtList());
  BeamISRWeight(subs,1);
  double tmp(m_last);
  m_last=0.0;
  for (size_t i=0;i<subs->size();++i) m_last+=(*subs)[i]->m_result;
  return m_last-tmp;
}

double Single_Process::Partonic2()
{
  return m_lastxs;
}

bool Single_Process::CalculateTotalXSec(const std::string &resultpath,
					const bool create) 
{ 
  p_int->Reset();
  SP(Phase_Space_Handler) psh(p_int->PSHandler());
  if (p_int->ISR()) {
    if (m_nin==2) {
      if (m_flavs[0].Mass()!=p_int->ISR()->Flav(0).Mass() ||
          m_flavs[1].Mass()!=p_int->ISR()->Flav(1).Mass()) {
        p_int->ISR()->SetPartonMasses(&m_flavs.front());
      }
    }
  }
  psh->InitCuts();
  if (p_int->ISR())
    p_int->ISR()->SetSprimeMin(psh->Cuts()->Smin());
  psh->CreateIntegrators();
  p_int->SetResultPath(resultpath);
  p_int->ReadResults();
  exh->AddTerminatorObject(p_int);
  psh->InitIncoming();
  double var(p_int->TotalVar());
  msg_Info()<<METHOD<<"(): Calculate xs for '"
            <<m_name<<"' ("<<(p_gen?p_gen->Name():"")<<")"<<std::endl;
  double totalxs(psh->Integrate()/rpa.Picobarn());
  if (!IsEqual(totalxs,p_int->TotalResult())) {
    msg_Error()<<"Result of PS-Integrator and summation do not coincide!\n"
	       <<"  '"<<m_name<<"': "<<totalxs
	       <<" vs. "<<p_int->TotalResult()<<std::endl;
  }
  if (p_int->TotalXS()>0.0) {
    p_int->SetTotal();
    if (var==p_int->TotalVar()) {
      exh->RemoveTerminatorObject(p_int);
      return 1;
    }
    p_int->StoreResults();
    exh->RemoveTerminatorObject(p_int);
    return 1;
  }
  exh->RemoveTerminatorObject(p_int);
  return 0;
}

void Single_Process::SetScale(const Scale_Setter_Arguments &args)
{
  if (IsMapped()) return;
  Scale_Setter_Arguments cargs(args);
  cargs.p_proc=this;
  cargs.p_cpls=&m_cpls;
  p_scale = Scale_Setter_Base::Scale_Getter_Function::
    GetObject(m_pinfo.m_scale=cargs.m_scale,cargs);
  if (p_scale==NULL) THROW(fatal_error,"Invalid scale scheme");
}

void Single_Process::SetKFactor(const KFactor_Setter_Arguments &args)
{
  if (IsMapped()) return;
  KFactor_Setter_Arguments cargs(args);
  cargs.p_proc=this;
  p_kfactor = KFactor_Setter_Base::KFactor_Getter_Function::
    GetObject(m_pinfo.m_kfactor=cargs.m_kfac,cargs);
  if (p_kfactor==NULL) THROW(fatal_error,"Invalid kfactor scheme");
}

void Single_Process::SetLookUp(const bool lookup)
{
  m_lookup=lookup;
}

bool Single_Process::Combinable
(const size_t &idi,const size_t &idj)
{
  return true;
}

const Flavour_Vector &Single_Process::
CombinedFlavour(const size_t &idij)
{
  static Flavour_Vector fls(1,kf_none);
  return fls;
}
