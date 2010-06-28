#include "PHASIC++/Process/Process_Base.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PDF/Main/Shower_Base.H"
#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;

ME_wgtinfo::ME_wgtinfo(): m_nx(0), m_w0(0.), p_wx(0), m_y1(1.), m_y2(1.), m_renscale(0.) {}

ME_wgtinfo::~ME_wgtinfo() {if (m_nx>0) delete[] p_wx;}

ME_wgtinfo& ME_wgtinfo::operator*= (const double scal) {
  m_w0*=scal;
  for (int i=0;i<m_nx;i++) p_wx[i]*=scal;
  return *this;
}

void ME_wgtinfo::Flip() {
  std::swap<double>(m_x1,m_x2);
  std::swap<double>(m_y1,m_y2);
  if (m_nx>=10) for (int i=0;i<4;i++) std::swap<double>(p_wx[i+2],p_wx[i+6]);
  if (m_nx>=18) for (int i=0;i<4;i++) std::swap<double>(p_wx[i+10],p_wx[i+14]);
}

void ME_wgtinfo::AddMEweights(int n) {
  m_nx=n;
  p_wx=new double[m_nx];
  for (int i=0;i<m_nx;i++) p_wx[i]=0.;
}

namespace ATOOLS {
  template <> Blob_Data<PHASIC::ME_wgtinfo*>::~Blob_Data() {}
  template class Blob_Data<PHASIC::ME_wgtinfo*>;
}


Process_Base::Process_Base():
  p_parent(NULL), p_selected(this), p_mapproc(NULL),
  p_int(new Process_Integrator(this)), 
  p_selector(new Combined_Selector(p_int)),
  p_cuts(NULL), p_gen(NULL), p_shower(NULL),
  p_scale(NULL), p_kfactor(NULL),
  m_nin(0), m_nout(0), 
  m_oqcd(0), m_oew(0),
  m_lookup(false), m_trigger(true) {}

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

Process_Base *Process_Base::Parent()
{ 
  if (p_parent && p_parent!=this) return p_parent->Parent();
  return this; 
}

bool Process_Base::GeneratePoint()
{
  return true;
}

void Process_Base::SetKFactorOn(const bool on)
{
  if (p_kfactor!=NULL) p_kfactor->SetOn(on);
}

double Process_Base::Differential(const Cluster_Amplitude &ampl) 
{
  Vec4D_Vector &p(p_int->Momenta());
  for (size_t i(0);i<ampl.NIn();++i) p[i]=-ampl.Leg(i)->Mom();
  for (size_t i(ampl.NIn());i<p.size();++i) p[i]=ampl.Leg(i)->Mom();
  SetKFactorOn(false);
  double res(this->Differential(p));
  SetKFactorOn(true);
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
  { return a.m_fl.IsFermion() && b.m_fl.IsFermion()
      && (!a.m_fl.IsAnti() && b.m_fl.IsAnti()); }
  bool operator()(const Cluster_Leg *a,const Cluster_Leg *b)
  { return a->Flav().IsFermion() && b->Flav().IsFermion()
      && (!a->Flav().IsAnti() && b->Flav().IsAnti()); }
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
  for (size_t i(0);i<info.m_ps.size();++i) SortFlavours(info.m_ps[i]);
}

void Process_Base::SortFlavours(Process_Info &pi)
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
  SortFlavours(pi.m_ii,&fmm);
  SortFlavours(pi.m_fi,&fmm);
}

void Process_Base::Init(const Process_Info &pi,
			BEAM::Beam_Spectra_Handler *const beamhandler,
			PDF::ISR_Handler *const isrhandler)
{
  m_pinfo=pi;
  m_nin=m_pinfo.m_ii.NExternal();
  m_nout=m_pinfo.m_fi.NExternal();
  m_flavs.resize(m_nin+m_nout);
  if (m_pinfo.m_ii.m_ps.size()>0 && m_pinfo.m_fi.m_ps.size()>0) {
    SortFlavours(m_pinfo);
    std::vector<Flavour> fl;
    m_pinfo.m_ii.GetExternal(fl);
    m_pinfo.m_fi.GetExternal(fl);
    if (fl.size()!=m_nin+m_nout) THROW(fatal_error,"Internal error");
    for (size_t i(0);i<fl.size();++i) m_flavs[i]=fl[i];
    m_name=GenerateName(m_pinfo.m_ii,m_pinfo.m_fi);
  }
  double massin=0.0, massout=0.0;
  for (size_t i=0;i<m_nin;++i) massin+=m_flavs[i].Mass();
  for (size_t i=m_nin;i<m_nin+m_nout;++i) massout+=m_flavs[i].Mass();
  p_int->SetISRThreshold(sqr(Max(massin,massout)));
  p_int->Initialize(beamhandler,isrhandler);
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
    name+="__QCD("+ToString(info.m_nloqcdtype)+")";
  if (info.m_nloewtype!=nlo_type::lo) 
    name+="__EW("+ToString(info.m_nloewtype)+")";
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
    name+="__QCD("+ToString(fi.m_nloqcdtype)+")";
  if (fi.m_nloewtype!=nlo_type::lo) 
    name+="__EW("+ToString(fi.m_nloewtype)+")";
  return name;
}

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
}

void Process_Base::SortFlavours(Cluster_Amplitude *const ampl)
{
  FMMap fmm;
  ClusterLeg_Vector il, fl;
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
      fl.push_back(ampl->Leg(i));
      int kfc(ampl->Leg(i)->Flav().Kfcode());
      if (fmm.find(kfc)==fmm.end()) fmm[kfc]=0;
      if (ampl->Leg(i)->Flav().IsFermion()) ++fmm[kfc];
    }
  SortFlavours(il,&fmm);
  SortFlavours(fl,&fmm);
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
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) 
    name+="__"+ampl->Leg(i)->Flav().IDName();
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

void Process_Base::SetUpThreading()
{
}

void Process_Base::FillOnshellConditions()
{
  if (!Selector()) return;
  Subprocess_Info info(m_pinfo.m_ii);
  info.Add(m_pinfo.m_fi);
  Decay_Info_Vector ids(info.GetDecayInfos());
  for(size_t i=0;i<ids.size();i++)
    if (ids[i].m_osd) Selector()->AddOnshellCondition
      (PSId(ids[i].m_id),sqr(ids[i].m_fl.Mass()));  
}

void Process_Base::SetSelector(const Selector_Key &key)
{
  p_selector->Initialize(key);
}

bool Process_Base::Trigger(const Vec4D_Vector &p)
{
  if (((m_pinfo.m_fi.m_nloqcdtype)&nlo_type::real) ||
      ((m_pinfo.m_fi.m_nloqcdtype)&nlo_type::rsub)) return NoJetTrigger(p);
  if (LookUp() && IsMapped()) return true;
  m_trigger=p_selector->Trigger(p);
  return m_trigger;
}
 
bool Process_Base::NoJetTrigger(const Vec4D_Vector &p)
{
  return p_selector->NoJetTrigger(p);
}

bool Process_Base::JetTrigger
(const Vec4D_Vector &p,const Flavour_Vector &fl,int n)
{
  return p_selector->JetTrigger(p,fl,n);
}

bool Process_Base::JetTrigger(const Vec4D_Vector &p)
{
  return p_selector->JetTrigger(p);
}

NLO_subevtlist *Process_Base::GetSubevtList()
{
  return NULL;
}

void Process_Base::BuildCuts(Cut_Data *const cuts)
{
  return p_selector->BuildCuts(cuts);
}

void Process_Base::UpdateCuts(const double &sp,const double &y,
			      Cut_Data *const cuts)
{
  return p_selector->UpdateCuts(sp,y,cuts);
}

void Process_Base::AddPoint(const double &value)
{
}
