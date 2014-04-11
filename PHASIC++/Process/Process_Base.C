#include "PHASIC++/Process/Process_Base.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Channels/BBar_Multi_Channel.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Decay_Info.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;

Process_Base::Process_Base():
  p_parent(NULL), p_selected(this), p_mapproc(NULL), p_sproc(NULL),
  p_int(new Process_Integrator(this)), p_selector(NULL),
  p_cuts(NULL), p_gen(NULL), p_shower(NULL), p_mc(NULL),
  p_scale(NULL), p_kfactor(NULL),
  m_nin(0), m_nout(0), 
  m_oqcd(0), m_oew(0), m_tinfo(1), m_mcmode(0), m_cmode(0),
  m_lookup(false), m_use_biweight(true), p_apmap(NULL)
{
  m_last=m_lastb=0.0;
}

Process_Base::~Process_Base() 
{
  if (p_kfactor) delete p_kfactor;
  if (p_scale) delete p_scale;
  delete p_selector;
  delete p_int;
}

Process_Base *Process_Base::Selected()
{ 
  if (!p_selected) return NULL;
  if (p_selected!=this) return p_selected->Selected();
  return this; 
}

bool Process_Base::SetSelected(Process_Base *const proc)
{
  if (proc==this) {
    p_selected=this;
    return true;
  }
  if (IsGroup())
    for (size_t i(0);i<Size();++i)
      if ((*this)[i]->SetSelected(proc)) {
	p_selected=(*this)[i];
	return true;
      }
  return false;
}

Process_Base *Process_Base::Parent()
{ 
  if (p_parent && p_parent!=this) return p_parent->Parent();
  return this; 
}

bool Process_Base::GeneratePoint()
{
  return true;
}

void Process_Base::AddPoint(const double &value)
{
}

bool Process_Base::ReadIn(const std::string &pid)
{
  return true;
}

void Process_Base::WriteOut(const std::string &pid)
{
}

void Process_Base::EndOptimize()
{
}

void Process_Base::MPISync()
{
}

void Process_Base::SetFixedScale(const std::vector<double> &s)
{
  if (IsMapped()) p_mapproc->SetFixedScale(s);
  if (p_scale!=NULL) p_scale->SetFixedScale(s);
}

void Process_Base::SetSelectorOn(const bool on)
{
  Selector()->SetOn(on);
}

void Process_Base::SetUseBIWeight(bool on)
{
  m_use_biweight=on;
}

double Process_Base::Differential(const Cluster_Amplitude &ampl,int mode) 
{
  Vec4D_Vector p(ampl.Legs().size());
  for (size_t i(0);i<ampl.NIn();++i) p[i]=-ampl.Leg(i)->Mom();
  if (mode&16) THROW(not_implemented,"Invalid mode");
  for (size_t i(ampl.NIn());i<p.size();++i) p[i]=ampl.Leg(i)->Mom();
  if (mode&64) {
    if (mode&1) return 1.0;
    return Trigger(p);
  }
  bool selon(Selector()->On());
  if (!Trigger(p)) {
    if ((mode&1) && selon) {
      SetSelectorOn(false);
      Trigger(p);
    }
  }
  if (mode&2) {
    std::vector<double> s(ScaleSetter(1)->Scales().size(),0.0);
    s[stp::fac]=ampl.MuF2();
    s[stp::ren]=ampl.MuR2();
    s[stp::res]=ampl.MuQ2();
    if (s.size()>stp::size+stp::res)
      s[stp::size+stp::res]=ampl.KT2();
    SetFixedScale(s);
  }
  if (mode&4) SetUseBIWeight(false);
  if (mode&128) this->GeneratePoint(); 
  double res(this->Differential(p));
  if (mode&32) {
    SP(Phase_Space_Handler) psh(Parent()->Integrator()->PSHandler());
    res*=psh->Weight(p);
  }
  if (mode&4) SetUseBIWeight(true);
  if (mode&2) SetFixedScale(std::vector<double>());
  if (Selector()->On()!=selon) SetSelectorOn(selon);
  return res;
}

bool Process_Base::IsGroup() const
{
  return false;
}

bool Process_Base::FillIntegrator
(Phase_Space_Handler *const psh)
{
  return false;
}

bool Process_Base::InitIntegrator
(Phase_Space_Handler *const psh)
{
  if (!p_sproc) return true;
  DEBUG_FUNC(m_name);
  SetBBarMC(new BBar_Multi_Channel(this,p_sproc,psh));
  psh->SetFSRIntegrator(p_mc);
  return true;
}

void Process_Base::UpdateIntegrator
(Phase_Space_Handler *const psh)
{
}

class Order_KF {
public:
  bool operator()(const Subprocess_Info &a,const Subprocess_Info &b)
  { return a.m_fl.Kfcode()<b.m_fl.Kfcode(); }
  bool operator()(const Cluster_Leg *a,const Cluster_Leg *b)
  { return a->Flav().Kfcode()<b->Flav().Kfcode(); }
};// end of class Order_KF

class Order_IsoWeak {
public:
  bool operator()(const Subprocess_Info &a,const Subprocess_Info &b)
  { return a.m_fl.IsDowntype() && b.m_fl.IsUptype(); }
  bool operator()(const Cluster_Leg *a,const Cluster_Leg *b)
  { return a->Flav().IsDowntype() && b->Flav().IsUptype(); }
};// end of class Order_IsoWeak

class Order_Anti {
public:
  bool operator()(const Subprocess_Info &a,const Subprocess_Info &b)
  { return !a.m_fl.IsAnti() && b.m_fl.IsAnti(); }
  bool operator()(const Cluster_Leg *a,const Cluster_Leg *b)
  { return !a->Flav().IsAnti() && b->Flav().IsAnti(); }
};// end of class Order_Anti

class Order_SVFT {
public:
  bool operator()(const Subprocess_Info &a,const Subprocess_Info &b) 
  {
    if (a.m_fl.IsScalar() && !b.m_fl.IsScalar()) return true;
    if (a.m_fl.IsVector() && !b.m_fl.IsScalar() && 
	!b.m_fl.IsVector()) return true;
    if (a.m_fl.IsFermion() && !b.m_fl.IsFermion() && 
	!b.m_fl.IsScalar() && !b.m_fl.IsVector()) return true;
    return false;
  }
  bool operator()(const Cluster_Leg *a,const Cluster_Leg *b) 
  {
    if (a->Flav().IsScalar() && !b->Flav().IsScalar()) return true;
    if (a->Flav().IsVector() && !b->Flav().IsScalar() && 
	!b->Flav().IsVector()) return true;
    if (a->Flav().IsFermion() && !b->Flav().IsFermion() && 
	!b->Flav().IsScalar() && !b->Flav().IsVector()) return true;
    return false;
  }
};// end of class Order_SVFT

class Order_Mass {
public:
  int operator()(const Subprocess_Info &a,const Subprocess_Info &b) 
  { return a.m_fl.Mass()>b.m_fl.Mass(); }
  int operator()(const Cluster_Leg *a,const Cluster_Leg *b) 
  { return a->Flav().Mass()>b->Flav().Mass(); }
};// end of class Order_Mass

class Order_InvMass {
public:
  int operator()(const Subprocess_Info &a,const Subprocess_Info &b) 
  { return a.m_fl.Mass()<b.m_fl.Mass(); }
  int operator()(const Cluster_Leg *a,const Cluster_Leg *b) 
  { return a->Flav().Mass()<b->Flav().Mass(); }
};// end of class Order_InvMass

class Order_Coupling {
public:
  int operator()(const Subprocess_Info &a,const Subprocess_Info &b) 
  { return !a.m_fl.Strong() && b.m_fl.Strong(); }
  int operator()(const Cluster_Leg *a,const Cluster_Leg *b) 
  { return !a->Flav().Strong() && b->Flav().Strong(); }
};// end of class Order_Coupling

class Order_Priority {
public:
  int operator()(const Subprocess_Info &a,const Subprocess_Info &b) 
  { return a.m_fl.Priority() > b.m_fl.Priority(); }
  int operator()(const Cluster_Leg *a,const Cluster_Leg *b) 
  { return a->Flav().Priority() > b->Flav().Priority(); }
};// end of class Order_Priority

class Order_Multiplicity {
  FMMap* p_fmm;
public:
  Order_Multiplicity(FMMap* fmm) {p_fmm=fmm;}
  int operator()(const Subprocess_Info &a,const Subprocess_Info &b)
  {
    if ((*p_fmm)[int(a.m_fl.Kfcode())]==0 || 
	(*p_fmm)[int(b.m_fl.Kfcode())]==0) return 0;
    if ((*p_fmm)[int(a.m_fl.Kfcode())]>
	(*p_fmm)[int(b.m_fl.Kfcode())]) return 1;
    return 0;
  }
  int operator()(const Cluster_Leg *a,const Cluster_Leg *b)
  {
    if ((*p_fmm)[int(a->Flav().Kfcode())]==0 || 
	(*p_fmm)[int(b->Flav().Kfcode())]==0) return 0;
    if ((*p_fmm)[int(a->Flav().Kfcode())]>
	(*p_fmm)[int(b->Flav().Kfcode())]) return 1;
    return 0;
  }
};// end of class Order_Multiplicity

void Process_Base::SortFlavours(Subprocess_Info &info,FMMap *const fmm)
{
  if (info.m_ps.empty()) return;
  ATOOLS::Flavour heaviest(kf_photon);
  for (size_t i(0);i<info.m_ps.size();++i) {
    if (info.m_ps[i].m_fl.Mass()>heaviest.Mass()) heaviest=info.m_ps[i].m_fl;
    else if (info.m_ps[i].m_fl.Mass()==heaviest.Mass() &&
	     !info.m_ps[i].m_fl.IsAnti()) heaviest=info.m_ps[i].m_fl;
  }
  std::stable_sort(info.m_ps.begin(),info.m_ps.end(),Order_KF());
  std::stable_sort(info.m_ps.begin(),info.m_ps.end(),Order_IsoWeak());
  std::stable_sort(info.m_ps.begin(),info.m_ps.end(),Order_Anti());
  std::stable_sort(info.m_ps.begin(),info.m_ps.end(),Order_SVFT());
  if (fmm) 
    std::stable_sort(info.m_ps.begin(),info.m_ps.end(),Order_Multiplicity(fmm));
  if (heaviest.IsAnti())  
    std::stable_sort(info.m_ps.begin(),info.m_ps.end(),Order_InvMass());
  else std::stable_sort(info.m_ps.begin(),info.m_ps.end(),Order_Mass());
  std::stable_sort(info.m_ps.begin(),info.m_ps.end(),Order_Coupling());
  std::stable_sort(info.m_ps.begin(),info.m_ps.end(),Order_Priority());
  for (size_t i(0);i<info.m_ps.size();++i) SortFlavours(info.m_ps[i]);
}

void Process_Base::SortFlavours(Process_Info &pi,const int mode)
{
  FMMap fmm;
  for (size_t i(0);i<pi.m_ii.m_ps.size();++i) {
    const Flavour *hfl=&pi.m_ii.m_ps[i].m_fl;
    if (fmm.find(int(hfl->Kfcode()))==fmm.end()) 
      fmm[int(hfl->Kfcode())]=0;
    if (hfl->IsFermion()) {
      fmm[int(hfl->Kfcode())]+=10;
      if (!hfl->IsAnti()) fmm[int(hfl->Kfcode())]+=10;
    }
  }
  for (size_t i(0);i<pi.m_fi.m_ps.size();++i) {
    const Flavour *hfl=&pi.m_fi.m_ps[i].m_fl;
    if (fmm.find(int(hfl->Kfcode()))==fmm.end()) 
      fmm[int(hfl->Kfcode())]=0;
    if (hfl->IsFermion()) fmm[int(hfl->Kfcode())]++;
  }
  if (mode&1) SortFlavours(pi.m_ii,&fmm);
  SortFlavours(pi.m_fi,&fmm);
}

void Process_Base::Init(const Process_Info &pi,
			BEAM::Beam_Spectra_Handler *const beamhandler,
			PDF::ISR_Handler *const isrhandler,const int mode)
{
  m_pinfo=pi;
  m_nin=m_pinfo.m_ii.NExternal();
  m_nout=m_pinfo.m_fi.NExternal();
  m_flavs.resize(m_nin+m_nout);
  if (m_pinfo.m_ii.m_ps.size()>0 && m_pinfo.m_fi.m_ps.size()>0) {
    if (!(mode&1)) SortFlavours(m_pinfo);
    std::vector<Flavour> fl;
    m_pinfo.m_ii.GetExternal(fl);
    m_pinfo.m_fi.GetExternal(fl);
    if (fl.size()!=m_nin+m_nout) THROW(fatal_error,"Internal error");
    for (size_t i(0);i<fl.size();++i) m_flavs[i]=fl[i];
    m_name=GenerateName(m_pinfo.m_ii,m_pinfo.m_fi);
    m_pinfo.m_fi.BuildDecayInfos(m_nin);
    m_decins=m_pinfo.m_fi.GetDecayInfos();
    if (IsGroup()) {
      if (m_pinfo.m_nminq>0 || m_pinfo.m_nmaxq<m_nin+m_nout) 
        m_name+="__NQ_"+ToString(m_pinfo.m_nminq)+
	  "-"+ToString(m_pinfo.m_nmaxq);
    }
  }
  double massin=0.0, massout=0.0;
  for (size_t i=0;i<m_nin;++i) massin+=m_flavs[i].Mass();
  for (size_t i=m_nin;i<m_nin+m_nout;++i) massout+=m_flavs[i].Mass();
  p_int->SetISRThreshold(Max(massin,massout));
  p_int->Initialize(beamhandler,isrhandler);
  m_symfac=m_pinfo.m_fi.FSSymmetryFactor();
}

std::string Process_Base::GenerateName(const Subprocess_Info &info) 
{
  std::string name(info.m_fl.IDName());
  if (info.m_fl.Kfcode()==kf_quark && info.m_fl.IsAnti()) name+="b";
  if (info.m_ps.empty()) return name;
  name+="["+GenerateName(info.m_ps.front());
  for (size_t i(1);i<info.m_ps.size();++i) 
    name+="__"+GenerateName(info.m_ps[i]);
  if (info.m_nloqcdtype!=nlo_type::lo) 
    name+="__QCD("+ToString(info.m_nloqcdtype)+info.m_sv+")";
  if (info.m_nloewtype!=nlo_type::lo) 
    name+="__EW("+ToString(info.m_nloewtype)+info.m_sv+")";
  return name+="]";
}

std::string Process_Base::GenerateName
(const Subprocess_Info &ii,const Subprocess_Info &fi) 
{
  char nii[3], nfi[3];
  if (sprintf(nii,"%i",(int)ii.NExternal())<=0)
    THROW(fatal_error,"Conversion error");
  if (sprintf(nfi,"%i",(int)fi.NExternal())<=0)
    THROW(fatal_error,"Conversion error");
  std::string name(nii+std::string("_")+nfi);
  for (size_t i(0);i<ii.m_ps.size();++i) name+="__"+GenerateName(ii.m_ps[i]);
  for (size_t i(0);i<fi.m_ps.size();++i) name+="__"+GenerateName(fi.m_ps[i]);
  if (fi.m_nloqcdtype!=nlo_type::lo) 
    name+="__QCD("+ToString(fi.m_nloqcdtype)+fi.m_sv+")";
  if (fi.m_nloewtype!=nlo_type::lo) 
    name+="__EW("+ToString(fi.m_nloewtype)+fi.m_sv+")";
  return name;
}

class Order_NDecay {
public:
  int operator()(const Decay_Info *a,const Decay_Info *b) 
  { return IdCount(a->m_id)>IdCount(b->m_id); }
};// end of class Order_NDecay

void Process_Base::SortFlavours
(std::vector<Cluster_Leg*> &legs,FMMap *const fmm)
{
  if (legs.empty()) return;
  ATOOLS::Flavour heaviest(kf_photon);
  for (size_t i(0);i<legs.size();++i) {
    if (legs[i]->Flav().Mass()>heaviest.Mass()) heaviest=legs[i]->Flav();
    else if (legs[i]->Flav().Mass()==heaviest.Mass() &&
	     !legs[i]->Flav().IsAnti()) heaviest=legs[i]->Flav();
  }
  std::stable_sort(legs.begin(),legs.end(),Order_KF());
  std::stable_sort(legs.begin(),legs.end(),Order_IsoWeak());
  std::stable_sort(legs.begin(),legs.end(),Order_Anti());
  std::stable_sort(legs.begin(),legs.end(),Order_SVFT());
  if (fmm) 
    std::stable_sort(legs.begin(),legs.end(),Order_Multiplicity(fmm));
  if (heaviest.IsAnti()) 
    std::stable_sort(legs.begin(),legs.end(),Order_InvMass());
  else std::stable_sort(legs.begin(),legs.end(),Order_Mass());
  std::stable_sort(legs.begin(),legs.end(),Order_Coupling());
  std::stable_sort(legs.begin(),legs.end(),Order_Priority());
}

void Process_Base::SortFlavours
(Cluster_Amplitude *const ampl,const int mode)
{
  FMMap fmm;
  DecayInfo_Vector cs;
  ClusterLeg_Vector il, fl;
  std::vector<int> dec(ampl->Legs().size(),0);
  std::map<size_t,ClusterLeg_Vector> dmap;
  for (size_t j(0);j<ampl->Decays().size();++j) {
    Decay_Info *cdi(ampl->Decays()[j]);
    size_t did(cdi->m_id), ndc(IdCount(did));
    for (size_t i(ampl->NIn());i<dec.size();++i)
      if (did&ampl->Leg(i)->Id()) {
	dec[i]=1;
	dmap[cdi->m_id].push_back(ampl->Leg(i));
      }
    bool core(true);
    for (size_t i(0);i<ampl->Decays().size();++i)
      if ((ampl->Decays()[i]->m_id&did) &&
	  IdCount(ampl->Decays()[i]->m_id)>ndc) {
	core=false;
	break;
      }
    if (!core) continue;
    int kfc(cdi->m_fl.Kfcode());
    if (fmm.find(kfc)==fmm.end()) fmm[kfc]=0;
    if (cdi->m_fl.IsFermion()) {
      fmm[kfc]+=10;
      if (!cdi->m_fl.IsAnti()) fmm[kfc]+=10;
    }
    cs.push_back(cdi);
  }
  for (size_t i(0);i<ampl->Legs().size();++i)
    if (i<ampl->NIn()) {
      ampl->Leg(i)->SetFlav(ampl->Leg(i)->Flav().Bar());
      il.push_back(ampl->Leg(i));
      int kfc(ampl->Leg(i)->Flav().Kfcode());
      if (fmm.find(kfc)==fmm.end()) fmm[kfc]=0;
      if (ampl->Leg(i)->Flav().IsFermion()) {
	fmm[kfc]+=10;
	if (!ampl->Leg(i)->Flav().IsAnti()) fmm[kfc]+=10;
      }
    }
    else {
      if (dec[i]) continue;
      fl.push_back(ampl->Leg(i));
      int kfc(ampl->Leg(i)->Flav().Kfcode());
      if (fmm.find(kfc)==fmm.end()) fmm[kfc]=0;
      if (ampl->Leg(i)->Flav().IsFermion()) ++fmm[kfc];
    }
  if (mode&1) SortFlavours(il,&fmm);
  for (size_t i(0);i<cs.size();++i) {
    ampl->CreateLeg(Vec4D(),cs[i]->m_fl,ColorID(),cs[i]->m_id);
    fl.push_back(ampl->Legs().back());
  }
  SortFlavours(fl,&fmm);
  if (cs.size()) {
    cs=ampl->Decays();
    std::sort(cs.begin(),cs.end(),Order_NDecay());
    while (fl.size()<ampl->Legs().size()-ampl->NIn())
      for (ClusterLeg_Vector::iterator
	     fit(fl.begin());fit!=fl.end();++fit)
	if (dmap.find((*fit)->Id())!=dmap.end()) {
	  ClusterLeg_Vector cl(dmap[(*fit)->Id()]);
	  size_t inc(0), ncd(IdCount((*fit)->Id()));
	  for (size_t i(0);i<cs.size();++i)
	    if (IdCount(cs[i]->m_id)<ncd &&
		(cs[i]->m_id&(*fit)->Id()) && (cs[i]->m_id&inc)==0) {
	    ampl->CreateLeg(Vec4D(),cs[i]->m_fl,ColorID(),cs[i]->m_id);
	    for (ClusterLeg_Vector::iterator
		   cit(cl.begin());cit!=cl.end();)
	      if (!((*cit)->Id()&cs[i]->m_id)) ++cit;
	      else cit=cl.erase(cit);
	    cl.push_back(ampl->Legs().back());
	    inc|=cs[i]->m_id;
	  }
	  SortFlavours(cl,&fmm);
	  (*fit)->Delete();
	  fit=fl.erase(fit);
	  fl.insert(fit,cl.begin(),cl.end());
	  ampl->Legs().pop_back();
	  break;
	}
  }
  for (size_t i(0);i<ampl->NIn();++i) {
    il[i]->SetFlav(il[i]->Flav().Bar());
    ampl->Legs()[i]=il[i];
  }
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i)
    ampl->Legs()[i]=fl[i-ampl->NIn()];
}

std::string Process_Base::GenerateName(const Cluster_Amplitude *ampl)
{
  char nii[3], nfi[3];
  if (sprintf(nii,"%i",(int)ampl->NIn())<=0)
    THROW(fatal_error,"Conversion error");
  if (sprintf(nfi,"%i",(int)(ampl->Legs().size()-ampl->NIn()))<=0)
    THROW(fatal_error,"Conversion error");
  std::string name(nii+std::string("_")+nfi);
  for (size_t i(0);i<ampl->NIn();++i) 
    name+="__"+ampl->Leg(i)->Flav().Bar().IDName();
  DecayInfo_Vector decs(ampl->Decays());
  std::sort(decs.begin(),decs.end(),Order_NDecay());
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
    std::string op, cl;
    for (size_t j(0);j<decs.size();++j) {
      Int_Vector ids(ID(decs[j]->m_id));
      if (ampl->Leg(i)->Id()==(1<<ids.front()))
	op+=ToString(decs[j]->m_fl)+"[";
      else if (ampl->Leg(i)->Id()==(1<<ids.back())) cl+="]";
    }
    name+="__"+op+ampl->Leg(i)->Flav().IDName()+cl;
  }
  return name;
}

std::string Process_Base::GenerateName(const NLO_subevt *sub,const size_t &nin)
{
  char nii[3], nfi[3];
  if (sprintf(nii,"%i",(int)nin)<=0)
    THROW(fatal_error,"Conversion error");
  if (sprintf(nfi,"%i",(int)(sub->m_n-nin))<=0)
    THROW(fatal_error,"Conversion error");
  std::string name(nii+std::string("_")+nfi);
  for (size_t i(0);i<sub->m_n;++i) name+="__"+sub->p_fl[i].IDName();
  return name;
}

void Process_Base::SetGenerator(ME_Generator_Base *const gen) 
{ 
  p_gen=gen; 
}

void Process_Base::SetShower(PDF::Shower_Base *const ps)
{
  p_shower=ps; 
}

void Process_Base::FillOnshellConditions()
{
  if (!Selector()) return;
  Subprocess_Info info(m_pinfo.m_ii);
  info.Add(m_pinfo.m_fi);
  for(size_t i=0;i<m_decins.size();i++)
    if (m_decins[i]->m_osd) Selector()->AddOnshellCondition
      (PSId(m_decins[i]->m_id),sqr(m_decins[i]->m_fl.Mass()));
}

void Process_Base::FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amp,
                                  std::vector<std::vector<Complex> >& cols)
{
  THROW(fatal_error, "Virtual function called.");
}

void Process_Base::SetSelector(const Selector_Key &key)
{
  if (IsMapped()) return;
  if (p_selector==NULL) p_selector = new Combined_Selector(p_int);
  p_selector->Initialize(key);
}

bool Process_Base::Trigger(const Vec4D_Vector &p)
{
  if (IsMapped() && LookUp()) return Selector()->Result();
  return Selector()->Trigger(p);
}
 
NLO_subevtlist *Process_Base::GetSubevtList()
{
  return NULL;
}

NLO_subevtlist *Process_Base::GetRSSubevtList()
{
  return NULL;
}

ME_wgtinfo *Process_Base::GetMEwgtinfo()
{
  return NULL;
}

void Process_Base::InitCuts(Cut_Data *const cuts)
{
  cuts->Init(m_nin,m_flavs);
}

void Process_Base::BuildCuts(Cut_Data *const cuts)
{
  if (IsMapped() && LookUp()) return;
  Selector()->BuildCuts(cuts);
}

void Process_Base::SetRBMap(Cluster_Amplitude *ampl)
{
}

void Process_Base::InitPSHandler
(const double &maxerr,const std::string eobs,const std::string efunc)
{
  p_int->SetPSHandler(new Phase_Space_Handler(p_int,maxerr));
  if (eobs!="") p_int->PSHandler()->SetEnhanceObservable(eobs);
  if (efunc!="") p_int->PSHandler()->SetEnhanceFunction(efunc);
} 

double Process_Base::LastPlus()
{
  if (IsGroup()) {
    double last=0.0;
    for (size_t i(0);i<Size();++i)
      last+=(*this)[i]->LastPlus();
    return last;
  }
  double last(Last());
  return last>0.0?last:0.0;
}

double Process_Base::LastMinus()
{
  if (IsGroup()) {
    double last=0.0;
    for (size_t i(0);i<Size();++i)
      last+=(*this)[i]->LastMinus();
    return last;
  }
  double last(Last());
  return last<0.0?last:0.0;
}

void Process_Base::FillProcessMap(NLOTypeStringProcessMap_Map *apmap)
{
  p_apmap=apmap;
  if (IsGroup()) {
    for (size_t i(0);i<Size();++i) (*this)[i]->FillProcessMap(apmap);
  }
  else {
    nlo_type::code nlot(m_pinfo.m_fi.m_nloqcdtype);
    std::string fname(m_name);
    size_t pos=fname.find("EW");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    pos=fname.find("QCD");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    if (nlot&nlo_type::vsub) nlot=nlo_type::vsub;
    if (nlot&nlo_type::rsub) nlot=nlo_type::rsub;
    if (apmap->find(nlot)==apmap->end())
      (*apmap)[nlot] = new StringProcess_Map();
    StringProcess_Map *cmap((*apmap)[nlot]);
    if (cmap->find(fname)==cmap->end()) (*cmap)[fname]=this;
  }
}

size_t Process_Base::SetMCMode(const size_t mcmode)
{
  size_t cmcmode(m_mcmode);
  m_mcmode=mcmode;
  if (IsGroup())
    for (size_t i(0);i<Size();++i)
      (*this)[i]->SetMCMode(mcmode);
  return cmcmode;
}

size_t Process_Base::SetClusterMode(const size_t cmode)
{
  size_t ccmode(m_cmode);
  m_cmode=cmode;
  if (IsGroup())
    for (size_t i(0);i<Size();++i)
      (*this)[i]->SetClusterMode(cmode);
  return ccmode;
}

void Process_Base::SetSProc(Process_Base *sproc)
{
  p_sproc=sproc;
  if (IsGroup())
    for (size_t i(0);i<Size();++i)
      (*this)[i]->SetSProc(sproc);
}

void Process_Base::SetBBarMC(BBar_Multi_Channel *mc)
{
  p_mc=mc;
  if (IsGroup())
    for (size_t i(0);i<Size();++i)
      (*this)[i]->SetBBarMC(mc);
}
