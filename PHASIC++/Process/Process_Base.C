#include "PHASIC++/Process/Process_Base.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Phase_Space_Integrator.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Channels/BBar_Multi_Channel.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Decay_Info.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_MPI.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/Color.H"
#include "ATOOLS/Math/Permutation.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <algorithm>

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

int Process_Base::s_usefmm(-1);

Process_Base::Process_Base():
  p_parent(NULL), p_selected(this), p_mapproc(NULL),
  p_sproc(NULL), p_caller(this),
  p_int(new Process_Integrator(this)), p_selector(NULL),
  p_cuts(NULL), p_gen(NULL), p_shower(NULL), p_nlomc(NULL), p_mc(NULL),
  p_scale(NULL), p_kfactor(NULL),
  m_nin(0), m_nout(0), m_maxcpl(2,99), m_mincpl(2,0), 
  m_mcmode(0), m_cmode(0),
  m_lookup(false), m_use_biweight(true),
  m_hasinternalscale(false), m_internalscale(sqr(rpa->gen.Ecms())),
  p_apmap(NULL)
{
  if (s_usefmm<0)
    s_usefmm =
      Settings::GetMainSettings()["PB_USE_FMM"].SetDefault(0).Get<int>();
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

size_t Process_Base::SynchronizeSelectedIndex(Process_Base & proc)
{
  size_t otherindex(proc.SelectedIndex());
  if (otherindex < Size()) {
    SetSelected((*this)[otherindex]);
    return otherindex;
  }
  return std::numeric_limits<size_t>::max();
}

size_t Process_Base::SelectedIndex()
{
  for (size_t i(0); i < Size(); ++i) {
    if ((*this)[i] == Selected()) {
      return i;
    }
  }
  return std::numeric_limits<size_t>::max();
}

Process_Base *Process_Base::Parent()
{ 
  if (p_parent && p_parent!=this) return p_parent->Parent();
  return this; 
}

bool Process_Base::GeneratePoint() { return true; }

void Process_Base::AddPoint(const double &value) {}

bool Process_Base::ReadIn(const std::string &pid) { return true; }

void Process_Base::WriteOut(const std::string &pid) { }

void Process_Base::EndOptimize() {}

void Process_Base::MPICollect(std::vector<double> &sv,size_t &i)
{
  if (IsGroup())
    for (size_t j(0);j<Size();++j)
      (*this)[j]->MPICollect(sv,i);
}

void Process_Base::MPIReturn(std::vector<double> &sv,size_t &i)
{
  if (IsGroup())
    for (size_t j(0);j<Size();++j)
      (*this)[j]->MPIReturn(sv,i);
}

void Process_Base::MPISync(const int mode)
{
  if (mode) return;
#ifdef USING__MPI
  size_t i(0), j(0);
  std::vector<double> sv;
  MPICollect(sv,i);
  if (mpi->Size()>1)
    mpi->Allreduce(&sv[0],sv.size(),MPI_DOUBLE,MPI_SUM);
  MPIReturn(sv,j);
#endif
}

void Process_Base::SetFixedScale(const std::vector<double> &s)
{
  if (IsMapped()) p_mapproc->SetFixedScale(s);
  if (p_scale!=NULL) p_scale->SetFixedScale(s);
}

void Process_Base::SetSelectorOn(const bool on) { Selector()->SetOn(on); }

void Process_Base::SetUseBIWeight(bool on) { m_use_biweight=on; }

Weights_Map Process_Base::Differential(const Cluster_Amplitude &ampl,
                                       Variations_Mode varmode,
                                       int mode)
{
  DEBUG_FUNC(this<<" -> "<<m_name<<", mode = "<<mode);
  msg_Debugging()<<ampl<<"\n";
  Vec4D_Vector p(ampl.Legs().size());
  for (size_t i(0);i<ampl.NIn();++i) p[i]=-ampl.Leg(i)->Mom();
  if (mode&16) THROW(not_implemented,"Invalid mode");
  for (size_t i(ampl.NIn());i<p.size();++i) p[i]=ampl.Leg(i)->Mom();
  if (mode&64) {
    if (mode&1) return {1.0};
    return {static_cast<double>(Trigger(p))};
  }
  bool selon(Selector()->On());
  if (mode&1) SetSelectorOn(false);
  if (!Trigger(p)) {
    if (Selector()->On()!=selon) SetSelectorOn(selon);
    return 0.0;
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
  if (mode&128) while (!this->GeneratePoint()); 
  else {
    std::shared_ptr<Color_Integrator> ci=
      Integrator()->ColorIntegrator();
    if (ci!=nullptr) ci->SetPoint(&ampl);
  }
  auto wgtmap = this->Differential(p, varmode);
  wgtmap/=m_issymfac;
  NLO_subevtlist *subs(this->GetSubevtList());
  if (subs) {
    (*subs)*=1.0/m_issymfac;
    (*subs).MultMEwgt(1.0/m_issymfac);
  }
  if (mode&32) {
    auto psh = Parent()->Integrator()->PSHandler();
    wgtmap*=psh->Weight(p);
  }
  if (mode&4) SetUseBIWeight(true);
  if (mode&2) SetFixedScale(std::vector<double>());
  if (Selector()->On()!=selon) SetSelectorOn(selon);
  return wgtmap;
}

bool Process_Base::IsGroup() const { return false; }
int  Process_Base::PerformTests()  { return 1; }
bool Process_Base::FillIntegrator(Phase_Space_Handler *const psh) { return false; }
void Process_Base::UpdateIntegrator(Phase_Space_Handler *const psh) {}

bool Process_Base::InitIntegrator
(Phase_Space_Handler *const psh)
{
  if (!p_sproc) return true;
  DEBUG_FUNC(m_name);
  SetBBarMC(new BBar_Multi_Channel(this,p_sproc,psh));
  psh->SetFSRIntegrator(p_mc);
  return true;
}


class Order_NDecay {
public:
  int operator()(const Decay_Info *a,const Decay_Info *b)
  { return IdCount(a->m_id)>IdCount(b->m_id); }
};// end of class Order_NDecay

void Process_Base::SortFlavours(Subprocess_Info &info,FMMap *const fmm)
{
  if (info.m_ps.empty()) return;
  ATOOLS::Flavour heaviest(kf_photon);
  for (size_t i(0);i<info.m_ps.size();++i) {
    if (info.m_ps[i].m_fl.Mass()>heaviest.Mass()) heaviest=info.m_ps[i].m_fl;
    else if (info.m_ps[i].m_fl.Mass()==heaviest.Mass() &&
	     !info.m_ps[i].m_fl.IsAnti()) heaviest=info.m_ps[i].m_fl;
  }
  std::sort(info.m_ps.begin(),info.m_ps.end(),Order_Flavour(fmm));
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
  if (mode&1) SortFlavours(pi.m_ii,s_usefmm?&fmm:NULL);
  SortFlavours(pi.m_fi,s_usefmm?&fmm:NULL);
}

bool Process_Base::InitScale()
{
  if (p_scale==NULL) return true;
  return p_scale->Initialize();
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
  m_issymfac=1.0;
  m_symfac=m_pinfo.m_fi.FSSymmetryFactor();
  if (m_nin==2 && m_flavs[0]==m_flavs[1] &&
      isrhandler->AllowSwap(m_flavs[0],m_flavs[1]))
    m_symfac*=(m_issymfac=2.0);
  m_name+=pi.m_addname;
  m_resname=m_name;
}

std::string Process_Base::BaseName
(const std::string &name,const std::string &addname)
{
  std::string fname(name);
  size_t len(addname.length());
  if (len) {
    size_t apos(fname.rfind(addname));
    if (apos!=std::string::npos) fname=fname.erase(apos,len);
  }
  size_t pos=fname.find("EW");
  if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
  pos=fname.find("QCD");
  if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
  return fname;
}

std::string Process_Base::GenerateName(const Subprocess_Info &info) 
{
  std::string name(info.m_fl.IDName());
  if (info.m_fl.Kfcode()==kf_quark && info.m_fl.IsAnti()) name+="b";
  if (info.m_ps.empty()) return name;
  name+="["+GenerateName(info.m_ps.front());
  for (size_t i(1);i<info.m_ps.size();++i) 
    name+="__"+GenerateName(info.m_ps[i]);
  if (info.m_nlotype!=nlo_type::lo) {
    if      (info.m_nlocpl[0]==1. && info.m_nlocpl[1]==0.) name+="__QCD(";
    else if (info.m_nlocpl[0]==0. && info.m_nlocpl[1]==1.) name+="__EW(";
    else                                                   name+="__QCDEW(";
    name+=ToString(info.m_nlotype)+info.m_sv+")";
  }
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
  if (fi.m_nlotype!=nlo_type::lo) {
    if      (fi.m_nlocpl[0]==1. && fi.m_nlocpl[1]==0.) name+="__QCD(";
    else if (fi.m_nlocpl[0]==0. && fi.m_nlocpl[1]==1.) name+="__EW(";
    else                                               name+="__QCDEW(";
    name+=ToString(fi.m_nlotype)+fi.m_sv+")";
  }
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
  std::sort(legs.begin(),legs.end(),Order_Flavour(fmm));
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
  if (mode&1) SortFlavours(il,s_usefmm?&fmm:NULL);
  for (size_t i(0);i<cs.size();++i) {
    ampl->CreateLeg(Vec4D(),cs[i]->m_fl,ColorID(),cs[i]->m_id);
    fl.push_back(ampl->Legs().back());
  }
  SortFlavours(fl,s_usefmm?&fmm:NULL);
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
	  SortFlavours(cl,s_usefmm?&fmm:NULL);
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

void Process_Base::SetNLOMC(PDF::NLOMC_Base *const mc)
{
  p_nlomc=mc; 
}

void Process_Base::FillOnshellConditions()
{
  if (!Selector()) return;
  Subprocess_Info info(m_pinfo.m_ii);
  info.Add(m_pinfo.m_fi);
  for(size_t i=0;i<m_decins.size();i++)
    if (m_decins[i]->m_osd) Selector()->AddOnshellCondition
      (m_decins[i]->m_id,sqr(m_decins[i]->m_fl.Mass()));
}

bool Process_Base::FillFinalState(const ATOOLS::Vec4D_Vector &p)
{
  return true;
}

std::vector<std::vector<int> > *Process_Base::Colours() const
{
  return NULL;
}

void Process_Base::FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amp,
                                  std::vector<std::vector<Complex> >& cols)
{
  THROW(fatal_error, "Virtual function called.");
}

void Process_Base::SetSelector(const Selector_Key &key)
{
  if (IsMapped()) return;
  if (p_selector==NULL) p_selector = new Combined_Selector(this);
  p_selector->Initialize(key);
}

void Process_Base::SetCaller(Process_Base *const proc)
{
  p_caller=proc;
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
  p_int->SetPSHandler(new Phase_Space_Handler(p_int, maxerr, eobs, efunc));
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
    nlo_type::code nlot(m_pinfo.m_fi.m_nlotype);
    std::string fname(m_name);
    size_t pos=m_pinfo.m_addname.length()?
      fname.find(m_pinfo.m_addname):std::string::npos;
    if (pos!=std::string::npos) fname=fname.substr(0,pos);
    pos=fname.find("EW");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    pos=fname.find("QCD");
    if (pos!=std::string::npos) fname=fname.substr(0,pos-2);
    if (nlot&nlo_type::vsub) nlot=nlo_type::vsub;
    if (nlot&nlo_type::rsub) nlot=nlo_type::rsub;
    if (apmap->find(nlot)==apmap->end())
      (*apmap)[nlot] = new StringProcess_Map();
    StringProcess_Map *cmap((*apmap)[nlot]);
    if (cmap->find(fname)!=cmap->end()) {
      if (msg_LevelIsDebugging()) {
        Process_Base* old = (*cmap)[fname];
        msg_Out()<<METHOD<<"(): replacing '"<<m_name<<"' "
                               <<Demangle(typeid(*old).name())
                               <<" -> "<<Demangle(typeid(*this).name())<<"\n";
      }
    }
    (*cmap)[fname]=this;
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

int Process_Base::NaiveMapping(Process_Base *proc) const
{
  DEBUG_FUNC(Name()<<" -> "<<proc->Name());
  const Vertex_Table *vt(s_model->VertexTable());
  std::map<Flavour,Flavour> fmap;
  std::vector<Flavour> curf(m_flavs);
  for (size_t i(0);i<curf.size();++i) fmap[curf[i]]=proc->m_flavs[i];
  for (std::map<Flavour,Flavour>::const_iterator
	 fit(fmap.begin());fit!=fmap.end();++fit)
    DEBUG_VAR(fit->first<<" -> "<<fit->second);
  for (size_t i(0);i<curf.size();++i) {
    Flavour cf(curf[i]), mf(fmap[cf]);
    if (cf==mf) continue;
    const Vertex_List &vlc(vt->find(cf)->second);
    const Vertex_List &vlm(vt->find(mf)->second);
    DEBUG_VAR(cf<<" "<<vlc.size());
    DEBUG_VAR(mf<<" "<<vlm.size());
    if (vlc.size()!=vlm.size()) return 0;
    for (size_t j(0);j<vlc.size();++j) {
      msg_Indent();
      DEBUG_VAR(*vlc[j]);
      bool match(false);
      for (size_t k(0);k<vlm.size();++k) {
	msg_Indent();
	DEBUG_VAR(*vlm[k]);
	if (vlm[k]->Compare(vlc[j])==0) {
	  msg_Indent();
	  for (int n=1;n<vlc[j]->NLegs();++n) {
	    std::map<Flavour,Flavour>::
	      const_iterator cit(fmap.find(vlc[j]->in[n]));
	    if (cit==fmap.end()) {
	      DEBUG_VAR(vlc[j]->in[n]<<" -> "<<vlm[k]->in[n]);
	      if (vlc[j]->in[n].Mass()!=vlm[k]->in[n].Mass() ||
		  vlc[j]->in[n].Width()!=vlm[k]->in[n].Width()) {
		msg_Debugging()<<"m_"<<vlc[j]->in[n]
			       <<" = "<<vlc[j]->in[n].Mass()
			       <<" / m_"<<vlm[k]->in[n]
			       <<" = "<<vlm[k]->in[n].Mass()<<"\n";
		msg_Debugging()<<"w_"<<vlc[j]->in[n]
			       <<" = "<<vlc[j]->in[n].Width()
			       <<" / w_"<<vlm[k]->in[n]
			       <<" = "<<vlm[k]->in[n].Width()<<"\n";
		return 0;
	      }
	      fmap[vlc[j]->in[n]]=vlm[k]->in[n];
	      curf.push_back(vlc[j]->in[n]);
	    }
	    else if (cit->second!=vlm[k]->in[n]) {
	      DEBUG_VAR(cit->second<<" "<<vlm[k]->in[n]);
	      return 0;
	    }
	  }
	  DEBUG_VAR(*vlc[j]);
	  DEBUG_VAR(*vlm[k]);
	  match=true;
	  break;
	}
      }
      if (!match) return 0;
    }
  }
  DEBUG_VAR("OK");
  return 1;
}

std::string Process_Base::ShellName(std::string name) const
{
  if (name.length()==0) name=m_name;
  for (size_t i(0);(i=name.find('-',i))!=std::string::npos;name.replace(i,1,"m"));
  for (size_t i(0);(i=name.find('+',i))!=std::string::npos;name.replace(i,1,"p"));
  for (size_t i(0);(i=name.find('~',i))!=std::string::npos;name.replace(i,1,"x"));
  for (size_t i(0);(i=name.find('(',i))!=std::string::npos;name.replace(i,1,"_"));
  for (size_t i(0);(i=name.find(')',i))!=std::string::npos;name.replace(i,1,"_"));
  for (size_t i(0);(i=name.find('[',i))!=std::string::npos;name.replace(i,1,"I"));
  for (size_t i(0);(i=name.find(']',i))!=std::string::npos;name.replace(i,1,"I"));
  return name;
}

void Process_Base::ConstructColorMatrix()
{
  DEBUG_FUNC(m_name);
  if (IsMapped()) return;
  std::string file(rpa->gen.Variable("SHERPA_CPP_PATH")
		   +"/Process/Sherpa/"+m_name+".col");
  My_In_File in(file);
  if (in.Open()) {
    int n;
    *in>>n;
    m_cols.m_perms.resize(n);
    m_cols.m_colfacs.resize(n,std::vector<double>(n));
    for (size_t i(0);i<m_cols.m_perms.size();++i) {
      *in>>n;
      m_cols.m_perms[i].resize(n);
      for (size_t j(0);j<n;++j)	*in>>m_cols.m_perms[i][j];
    }
    for (size_t i(0);i<m_cols.m_colfacs.size();++i)
      for (size_t j(i);j<m_cols.m_colfacs[i].size();++j) {
	*in>>m_cols.m_colfacs[i][j];
	m_cols.m_colfacs[j][i]=m_cols.m_colfacs[i][j];
      }
    std::string check;
    *in>>check;
    if (check!="eof") THROW(fatal_error,"Corrupted color file "+file);
    in.Close();
    return;
  }
  Flavour_Vector fls(m_pinfo.ExtractFlavours());
  int n(0);
  for (size_t i(0);i<fls.size();++i) if (fls[i].Strong()) ++n;
  if (n==0) return;
  for (size_t i(0);i<m_nin;++i) fls[i]=fls[i].Bar();
  m_cols=ColorMatrix(fls);
  My_Out_File out(file);
  if (!out.Open()) THROW(fatal_error,"Cannot open '"+file+"'");
  out->precision(12);
  *out<<m_cols.m_perms.size()<<"\n";
  for (size_t i(0);i<m_cols.m_perms.size();++i) {
    *out<<m_cols.m_perms[i].size();
    for (size_t j(0);j<m_cols.m_perms[i].size();++j)
      *out<<" "<<m_cols.m_perms[i][j];
    *out<<"\n";
  }
  for (size_t i(0);i<m_cols.m_colfacs.size();++i) {
    for (size_t j(i);j<m_cols.m_colfacs[i].size();++j)
      *out<<m_cols.m_colfacs[i][j]<<" ";
    *out<<"\n";
  }
  *out<<"eof\n";
}

Color_Matrix Process_Base::ColorMatrix(const Flavour_Vector &fls) const
{
  DEBUG_FUNC(fls);
  int iq(-1), np(0), nf(0);
  std::vector<int> sids, iqbs;
  for (size_t i(0);i<fls.size();++i) {
    if (fls[i].StrongCharge()>0) {
      if (fls[i].StrongCharge()==8) sids.push_back(i);
      else {
	if (iq>=0) sids.push_back(i);
	else iq=i;
	++nf;
      }
    }
    else if (fls[i].StrongCharge()<0) {
      sids.push_back(i);
      iqbs.push_back(i);
    }
  }
  bool adjoint(iq<0);
  if (iq<0) {
    iq=sids.front();
    sids.erase(sids.begin());
    iqbs.push_back(sids.back());
  }
  int unique(iqbs.size()==1?1:0);
  if (unique) {
    for (std::vector<int>::iterator
	   sidit(sids.begin());sidit!=sids.end();++sidit)
      if (*sidit==iqbs.front()) sidit=--sids.erase(sidit);
    sids.push_back(iqbs.front());
  }
  std::vector<int> idr(sids.size()+1), idc(sids.size()+1);
  idc.front()=idr.front()=iq;
  if (unique) idc.back()=idr.back()=iqbs.front();
  Permutation perms(sids.size()-unique);
  std::vector<int> act(perms.MaxNumber(),0);
  std::map<int,std::map<int,double> > cij2;
  for (size_t i(0);i<perms.MaxNumber();++i) {
    int *cur(perms.Get(i));
    for (size_t k(0);k<sids.size()-unique;++k)
      idr[k+1]=sids[cur[k]];
    int valid(false);
    for (size_t l(0);l<iqbs.size();++l)
      if (idr.back()==iqbs[l]) valid=true;
    if (!valid) continue;
    for (size_t j(i);j<perms.MaxNumber();++j) {
      int *cur(perms.Get(j));
      for (size_t k(0);k<sids.size()-unique;++k)
	idc[k+1]=sids[cur[k]];
      int valid(false), lqr(0), lqc(0);
      for (size_t l(0);l<iqbs.size();++l)
	if (idc.back()==iqbs[l]) valid=true;
      if (!valid) continue;
      Expression cij(200,100);
      cij.SetTR(1.);
      cij.pop_back();
      size_t lr(adjoint?iq:cij.FIndex()), fr(idc.back());
      for (size_t k(0);k<idr.size();++k) {
	if (fls[idr[k]].StrongCharge()==3) {
	  if (lqr>0) valid=false;
	  cij.push_back(Delta::New(idr[k],lr));
	  lqr=1;
	}
	else if (fls[idr[k]].StrongCharge()==-3) {
	  if (lqr<0) valid=false;
	  cij.push_back(Delta::New(lr,idr[k]));
	  lr=cij.FIndex();
	  lqr=-1;
	}
	else if (!adjoint) {
	  if (lqr<0) valid=false;
	  size_t nr(cij.FIndex());
	  cij.push_back(Fundamental::New(idr[k],lr,nr));
	  lr=nr;
	}
	else {
	  if (k==0 || k==idr.size()-1) continue;
	  size_t nr(k==idr.size()-2?fr:cij.AIndex());
	  cij.push_back(Adjoint::New(idr[k],lr,nr));
	  lr=nr;
	}
      }
      size_t lc(adjoint?iq:cij.FIndex()), fc(idc.back());
      for (size_t k(0);k<idc.size();++k) {
	if (fls[idc[k]].StrongCharge()==3) {
	  if (lqc>0) valid=false;
	  cij.push_back(Delta::New(lc,idc[k]));
	  lqc=1;
	}
	else if (fls[idc[k]].StrongCharge()==-3) {
	  if (lqc<0) valid=false;
	  cij.push_back(Delta::New(idc[k],lc));
	  lc=cij.FIndex();
	  lqc=-1;
	}
	else if (!adjoint) {
	  if (lqc<0) valid=false;
	  size_t nc(cij.FIndex());
	  cij.push_back(Fundamental::New(idc[k],nc,lc));
	  lc=nc;
	}
	else {
	  if (k==0 || k==idc.size()-1) continue;
	  size_t nc(k==idc.size()-2?fc:cij.AIndex());
	  cij.push_back(Adjoint::New(idc[k],lc,nc));
	  lc=nc;
	}
      }
      if (!valid) continue;
      if (act[i]==0) act[i]=++np;
      if (act[j]==0) act[j]=++np;
      if (msg->LevelIsDebugging()) cij.Print();
      cij.Evaluate();
      if (cij.Result().imag())
	msg_Error()<<METHOD<<"(): Non-real color coefficient "
		   <<cij.Result()<<std::endl;
      cij2[i][j]=cij.Result().real();
    }
  }
  double N(nf>1?sqr(Factorial(nf-1)):1.0);
  Color_Matrix cij;
  cij.m_perms.resize(np);
  cij.m_colfacs.resize(np,std::vector<double>(np));
  for (size_t i(0);i<act.size();++i)
    if (act[i]) {
      int *cur(perms.Get(i));
      for (size_t k(0);k<sids.size()-unique;++k)
	idc[k+1]=sids[cur[k]];
      cij.m_perms[act[i]-1]=idc;
      for (size_t j(i);j<act.size();++j)
	if (act[j]) {
	  cij.m_colfacs[act[i]-1][act[j]-1]=
	    cij.m_colfacs[act[j]-1][act[i]-1]=cij2[i][j]/N;
	}
    }
  for (size_t i(0);i<cij.m_perms.size();++i)
    msg_Debugging()<<i<<" "<<cij.m_perms[i]<<"\n";
  for (size_t i(0);i<cij.m_colfacs.size();++i) {
    for (size_t j(i);j<cij.m_colfacs[i].size();++j)
      msg_Debugging()<<cij.m_colfacs[i][j]<<" ";
    msg_Debugging()<<"\n";
  }
  return cij;
}
