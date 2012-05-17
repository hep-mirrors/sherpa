#include "PHASIC++/Process/Single_Process.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Process/POWHEG_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
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

Weight_Info *Single_Process::OneEvent(const int wmode,const int mode)
{
  p_selected=this;
  return p_int->PSHandler()->OneEvent(this,mode);
}

double Single_Process::KFactor() const
{
  if (p_kfactor) return p_kfactor->KFactor();
  return 1.0;
}

double Single_Process::BeamISRWeight
(const double& Q2,const int mode) const
{
  if (!m_use_biweight) return 1.;
  double wgt(1.0);
  if (m_nin!=2) return 0.5/p_int->Momenta()[0].Mass();
  if (p_int->ISR()) {
    wgt*=p_int->ISR()->Weight
      (mode,p_int->Momenta()[0],p_int->Momenta()[1],
       Q2,Q2,m_flavs[0],m_flavs[1]);
    double LQ2(Q2);
    ClusterAmplitude_Vector &ampls
      ((IsMapped()?p_mapproc:this)->ScaleSetter()->Amplitudes());
    if (ampls.size()) {
      DEBUG_FUNC(m_name<<", mode = "<<mode);
      Cluster_Amplitude *ampl(ampls.front());
      if (m_pinfo.Has(nlo_type::real)) ampl=ampl->Next();
      for (;ampl;ampl=ampl->Next()) {
	if (IsEqual(LQ2,ampl->KT2())) continue;
	if (ampl->Next()) {
	  if (ampl->Next()->Splitter()->Stat()==3) {
	    msg_Debugging()<<"Skip decay "<<
	      ID(ampl->Next()->Splitter()->Id())<<"\n";
	    continue;
	  }
	}
	Flavour f1(ampl->Leg(0)->Flav().Bar());
	Flavour f2(ampl->Leg(1)->Flav().Bar());
	if (MapProc() && LookUp()) {
	  f1=ReMap(f1,ampl->Leg(0)->Id());
	  f2=ReMap(f2,ampl->Leg(1)->Id());
	}
	msg_Debugging()<<"PDF ratio "<<f1<<"("<<ampl->Leg(0)->Flav().Bar()
		       <<"),"<<f2<<"("<<ampl->Leg(1)->Flav().Bar()
		       <<") at "<<sqrt(LQ2);
	double wd1=p_int->ISR()->Weight
	  (mode|2,-ampl->Leg(0)->Mom(),-ampl->Leg(1)->Mom(),LQ2,LQ2,f1,f2,0);
	double wd2=p_int->ISR()->Weight
	  (mode|4,-ampl->Leg(0)->Mom(),-ampl->Leg(1)->Mom(),LQ2,LQ2,f1,f2,0);
	LQ2=ampl->KT2();
	double wn1=p_int->ISR()->Weight
	  (mode|2,-ampl->Leg(0)->Mom(),-ampl->Leg(1)->Mom(),LQ2,LQ2,f1,f2,0);
	double wn2=p_int->ISR()->Weight
	  (mode|4,-ampl->Leg(0)->Mom(),-ampl->Leg(1)->Mom(),LQ2,LQ2,f1,f2,0);
	if (!IsZero(wn1) && !IsZero(wd1)) wgt*=wn1/wd1;
	if (!IsZero(wn2) && !IsZero(wd2)) wgt*=wn2/wd2;
	msg_Debugging()<<" / "<<sqrt(LQ2)<<" -> "
		       <<wn1/wd1<<" * "<<wn2/wd2<<" ( "<<wgt<<" )\n";
      }
    }
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
  double muf2(subs->back()->m_mu2[stp::fac]);
  if (m_nin==2 && p_int->ISR()) {
    size_t nscales(0);
    for (size_t i(0);i<subs->size();++i) {
      NLO_subevt *sub((*subs)[i]);
      if (!IsEqual(sub->m_mu2[stp::fac],muf2) && sub->m_me!=0.0) {
        sub->m_result+=sub->m_last[mode]=
          sub->m_me*BeamISRWeight(sub->m_mu2[stp::fac],mode);
	++nscales;
      }
    }
    if (nscales<subs->size()) {
      double lumi(BeamISRWeight(muf2,mode));
      for (size_t i(0);i<subs->size();++i) {
	if (IsEqual((*subs)[i]->m_mu2[stp::fac],muf2) &&
	    (*subs)[i]->m_me!=0.0) {
          (*subs)[i]->m_result+=(*subs)[i]->m_last[mode]=
            (*subs)[i]->m_me*lumi;
        }
      }
    }
  }
  else {
    for (size_t i(0);i<subs->size();++i) {
      (*subs)[i]->m_result+=(*subs)[i]->m_last[mode]=
        (*subs)[i]->m_me*BeamISRWeight((*subs)[i]->m_mu2[stp::fac],mode);
    }
  }
}

double Single_Process::Differential(const Vec4D_Vector &p)
{
  m_wgtinfo.m_w0=m_last[0]=0.0;
  p_int->SetMomenta(p);
  double flux=0.25/sqrt(sqr(p[0]*p[1])-p[0].Abs2()*p[1].Abs2());
  if (GetSubevtList()==NULL) {
    if (m_zero) return 0.0;
    Scale_Setter_Base *scs(ScaleSetter());
    if (IsMapped()) {
      p_mapproc->Integrator()->SetMomenta(p);
      scs=p_mapproc->ScaleSetter();
    }
    scs->SetCaller(this);
    if (Partonic(p,0)==0.0) return 0.0;
    if (m_wgtinfo.m_nx==0) m_wgtinfo.m_w0 = m_lastxs;
    m_wgtinfo*=flux;
    m_wgtinfo.m_mur2=scs->Scale(stp::ren);
    if (m_lastxs==0.0) return m_last[0]=0.0;
    return m_last[0]=m_lastxs*BeamISRWeight(scs->Scale(stp::fac),0);
  }
  Partonic(p,0);
  NLO_subevtlist *subs(GetSubevtList());
  BeamISRWeight(subs,0);
  for (size_t i=0;i<subs->size();++i) {
    m_last[0]+=(*subs)[i]->m_last[0];
    (*subs)[i]->m_mewgt*=flux;
  }
  return m_last[0];
}

double Single_Process::Differential2()
{
  m_last[1]=0.0;
  if (m_nin!=2 || p_int->ISR()==NULL ||
      p_int->ISR()->On()!=3) return 0.0;
  if (m_flavs[0]==m_flavs[1]) return 0.0;
  if (!p_int->ISR()->PDF(0)->Contains(m_flavs[1]) ||
      !p_int->ISR()->PDF(1)->Contains(m_flavs[0])) return 0.0;
  if (GetSubevtList()==NULL) {
    if (m_lastxs==0.0) return 0.0;
    Scale_Setter_Base *scs((IsMapped()?p_mapproc:this)->ScaleSetter());
    scs->SetCaller(this);
    m_last[1]=Partonic(p_int->Momenta(),1);
    if (m_last[1]!=0.0) m_last[1]*=BeamISRWeight(scs->Scale(stp::fac),1);
    return m_last[1];
  }
  Partonic(p_int->Momenta(),1);
  NLO_subevtlist *subs(GetSubevtList());
  BeamISRWeight(subs,1);
  for (size_t i=0;i<subs->size();++i) m_last[1]+=(*subs)[i]->m_last[1];
  return m_last[1];
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
        p_int->ISR()->SetPartonMasses(m_flavs);
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
  double totalxs(psh->Integrate()/rpa->Picobarn());
  if (!IsEqual(totalxs,p_int->TotalResult())) {
    msg_Error()<<"Result of PS-Integrator and summation do not coincide!\n"
	       <<"  '"<<m_name<<"': "<<totalxs
	       <<" vs. "<<p_int->TotalResult()<<std::endl;
  }
  if (p_int->Points()) {
    p_int->SetTotal();
    if (var==p_int->TotalVar()) {
      exh->RemoveTerminatorObject(p_int);
      return 1;
    }
    p_int->StoreResults();
    p_int->StoreBackupResults();
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
  m_pinfo.m_kfactor=cargs.m_kfac;
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

ATOOLS::ME_wgtinfo *Single_Process::GetMEwgtinfo()
{
  return &m_wgtinfo; 
}

ATOOLS::Flavour Single_Process::ReMap
(const ATOOLS::Flavour &fl,const size_t &id) const
{
  return fl;
}

Cluster_Amplitude *Single_Process::Cluster
(const size_t &mode,const double &kt2)
{
  MCatNLO_Process *mp(dynamic_cast<MCatNLO_Process*>(Parent()));
  if (mp) {
    Cluster_Amplitude *ampl(mp->GetAmplitude());
    if (ampl) return ampl;
  }
  POWHEG_Process *pp(dynamic_cast<POWHEG_Process*>(Parent()));
  if (pp) {
    Cluster_Amplitude *ampl(pp->GetAmplitude());
    if (ampl) return ampl;
  }
  if (!(mode&256)) {
    ClusterAmplitude_Vector &ampls
      ((IsMapped()?p_mapproc:this)->ScaleSetter()->Amplitudes());
    if (ampls.size()) {
      msg_Debugging()<<METHOD<<"(): Found "
		     <<ampls.size()<<" amplitude(s) ... ";
      if (p_int->InSwaped()) {
	msg_Debugging()<<"select 2nd.\n";
	Cluster_Amplitude *ampl(ampls.front()->CopyAll());
	for (Cluster_Amplitude *campl(ampl);
	     campl;campl=campl->Next()) {
	  std::swap<Cluster_Leg*>(campl->Legs()[0],campl->Legs()[1]);
	  for (size_t i(0);i<campl->Legs().size();++i) {
	    const Vec4D &p(campl->Leg(i)->Mom());
	    campl->Leg(i)->SetMom(Vec4D(p[0],-p[1],-p[2],-p[3]));
	  }
	}
	return ampl;
      }
      msg_Debugging()<<"select 1st.\n";
      return ampls.front()->CopyAll();
    }
  }
  PDF::Cluster_Definitions_Base* cd=p_shower->GetClusterDefinitions();
  int amode=cd->AMode(), cmode=mode;
  if (amode) cmode|=512;
  if (mode&512) cd->SetAMode(1);
  p_gen->SetClusterDefinitions(cd);
  Cluster_Amplitude* ampl(p_gen->ClusterConfiguration(this,cmode,kt2));
  if (ampl) ampl->Decays()=m_pinfo.m_fi.GetDecayInfos();
  cd->SetAMode(amode);
  return ampl;
}
