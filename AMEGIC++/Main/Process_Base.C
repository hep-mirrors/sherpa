#include "AMEGIC++/Main/Process_Base.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/KF_Table.H"

using namespace AMEGIC;
using namespace PHASIC;
using namespace ATOOLS;

int AMEGIC::Process_Base::s_gauge=10;

AMEGIC::Process_Base::Process_Base(): 
  p_model(NULL),
  p_b(0), m_eoreset(0), p_pl(0), 
  m_print_graphs(""), p_testmoms(0),
  m_Norm(1.), m_sfactor(1.), m_lastdxs(0.), m_lastk(1.)
{
  p_subevtlist=NULL;
  p_channellibnames = new std::list<std::string>();
  static int allowmap(-1);
  if (allowmap<0) {
    Scoped_Settings amegicsettings{
      Settings::GetMainSettings()["AMEGIC"] };
    allowmap = amegicsettings["ALLOW_MAPPING"].Get<int>();
    if (allowmap!=1) msg_Info()<<METHOD<<"(): Disable process mapping.\n";
  }
  m_allowmap=allowmap;
  m_lastk=1.0;
}

AMEGIC::Process_Base::~Process_Base()
{
  delete p_channellibnames;
  if (p_pl) delete [] p_pl;
  if (p_b) delete [] p_b;
}

double AMEGIC::Process_Base::Result()
{
  return 0.0;
}

std::string AMEGIC::Process_Base::LibName()
{
  return "error";
}

void AMEGIC::Process_Base::SetPrintGraphs(std::string gpath)
{
  m_print_graphs=gpath;
}

void AMEGIC::Process_Base::Init()
{
  p_pinfo = Translate(m_pinfo);
  p_pl = new Pol_Info[NIn()+NOut()];
  for (size_t i(0);i<m_pinfo.m_ii.m_ps.size();++i)
    p_pl[i]=ExtractPolInfo(m_pinfo.m_ii.m_ps[i]);
  p_pinfo->GetTotalPolList(p_pl+NIn());
  m_mincpl.resize(m_pinfo.m_mincpl.size());
  for (size_t i(0);i<m_mincpl.size();++i) {
    m_mincpl[i]=m_pinfo.m_mincpl[i];
    if (m_mincpl[i]!=m_pinfo.m_mincpl[i])
      THROW(not_implemented,"Non-integer couplings not supported by Amegic");
  }
  m_maxcpl.resize(m_pinfo.m_maxcpl.size());
  for (size_t i(0);i<m_maxcpl.size();++i) {
    m_maxcpl[i]=m_pinfo.m_maxcpl[i];
    if (m_maxcpl[i]!=m_pinfo.m_maxcpl[i])
      THROW(not_implemented,"Non-integer couplings not supported by Amegic");
  }
  p_b    = new int[NIn()+NOut()];
  for (size_t i=0;i<NIn();i++) p_b[i] = -1;
  for (size_t i=NIn();i<NIn()+NOut();i++) p_b[i] = 1;
}


#define PTS long unsigned int
#define PT(ARG) (PTS)(ARG)

typedef PHASIC::Single_Channel *(*Lib_Getter_Function)
  (int nin,int nout,ATOOLS::Flavour* fl,
   ATOOLS::Integration_Info * const info,PHASIC::Phase_Space_Handler *psh);

PHASIC::Single_Channel *LoadChannel(int nin,int nout,ATOOLS::Flavour* fl,
			    std::string& pID,PHASIC::Phase_Space_Handler *psh)
{
  size_t pos(pID.find("/"));
  s_loader->AddPath(rpa->gen.Variable("SHERPA_LIB_PATH"));
  Lib_Getter_Function gf = (Lib_Getter_Function)
    PT(s_loader->GetLibraryFunction("Proc_"+pID.substr(0,pos),
				    "Getter_"+pID.substr(pos+1)));
  if (gf==NULL) return NULL;
  return gf(nin,nout,fl,psh->GetInfo(),psh);
}

void AMEGIC::Process_Base::RequestVariables(Phase_Space_Handler *const psh)
{
} 

bool AMEGIC::Process_Base::FillIntegrator(Phase_Space_Handler *const psh)
{
  if (p_channellibnames->empty()) return true;
  Multi_Channel *mc(psh->FSRIntegrator());
  for (std::list<std::string>::iterator it(p_channellibnames->begin());
       it!=p_channellibnames->end();++it) {
    Single_Channel *sc = LoadChannel(NIn(),NOut(),(Flavour*)&Flavours().front(),
				     *it,&*Integrator()->PSHandler());
    if (sc==0) THROW(critical_error,"PS integration channels not compiled");
    sc->SetName(*it);
    mc->Add(sc);
  }
  return false;
}

double AMEGIC::Process_Base::SymmetryFactors()
{
  Flavour* flav  = new Flavour[m_nin+m_nout];
  double sym = 1.;
  int ndecays=p_pinfo->Ndecays();
  for (int i=0;i<=ndecays;i++) {
    int j=i;
    Process_Tags* pi=p_pinfo->GetDecay(j);
    size_t nout=pi->GetStableFlavList(flav);
    sym*=SBSymmetryFactor(flav,nout);
  }
  delete[] flav;
  return sym;
}

double AMEGIC::Process_Base::SBSymmetryFactor(Flavour* flout,size_t nout)
{
  double sym = 1.;
  for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
      kfit!=s_kftable.end();++kfit) {
    Flavour hflav(kfit->first);
    if (hflav.IsHadron()) continue; 
    int cp  = 0;
    int cap = 0;
    for (size_t j=0;j<nout;j++) {
      if (flout[j]==hflav)                                      cp++;
      else {
	if ((flout[j]==hflav.Bar()) && (hflav != hflav.Bar()))  cap++;
      }
    }
    if (cp>1)  sym *= double(Factorial(cp));
    if (cap>1) sym *= double(Factorial(cap));
  } 
  return 1./sym;
}

bool AMEGIC::Process_Base::CheckMapping(const Process_Base * proc)
{
  const ATOOLS::Flavour_Vector &flavs(Flavours());
  const ATOOLS::Flavour_Vector &partner_flavs(proc->Flavours());
  // create map
  std::map<ATOOLS::Flavour,ATOOLS::Flavour> flmap;
  for (size_t i=0;i<NIn()+NOut();++i) {
    if (flmap.find(partner_flavs[i])==flmap.end()) {
      flmap[partner_flavs[i]]=flavs[i];
      if (partner_flavs[i]!=(Flavour(partner_flavs[i])).Bar()) {
	flmap[(Flavour(partner_flavs[i])).Bar()]=(Flavour(flavs[i])).Bar();
      }
    }
  }
  // check map
  for (size_t i=0;i<NIn()+NOut();++i) {
    if (flmap[partner_flavs[i]]!=flavs[i]) {
      msg_Tracking()<<" mapping test failed "<<std::endl;
      return false;
    }
  }
  return true;
}

void AMEGIC::Process_Base::InitFlavmap(const Process_Base * proc)
{
  const ATOOLS::Flavour_Vector &flavs(Flavours());
  const ATOOLS::Flavour_Vector &partner_flavs(proc->Flavours());

  for (size_t i=0;i<NIn()+NOut();++i) {
    if (m_eflmap.find(partner_flavs[i])==m_eflmap.end()) {
      m_eflmap[partner_flavs[i]]=flavs[i];
      if (partner_flavs[i]!=(Flavour(partner_flavs[i])).Bar()) {
	m_eflmap[(Flavour(partner_flavs[i])).Bar()]=(Flavour(flavs[i])).Bar();
      }
    }
  }
}

void AMEGIC::Process_Base::AddtoFlavmap(const std::string& id,const Flavour&f1)
{
  if (m_fmap.find(id)==m_fmap.end()) {
    m_fmap[id]=f1;
  }
  else if (m_fmap[id]!=f1) THROW(critical_error,"Flavour mapping not unique!");
}

void AMEGIC::Process_Base::PrintProcessSummary(int it)
{
 for(int i=0;i<it;i++) std::cout<<"  ";
  std::cout<<Name()<<std::endl;
}


ATOOLS::Flavour AMEGIC::Process_Base::ReMap(const ATOOLS::Flavour& f0,const std::string& id) const
{
  if ((Partner()==NULL)||(Partner()==this)) return f0;
  std::map<std::string,ATOOLS::Flavour>::const_iterator fit(m_fmap.find(id));
  if (fit!=m_fmap.end()) return fit->second;

  Flavour_Map::const_iterator efit(m_eflmap.find(f0));
  if (efit!=m_eflmap.end()) return efit->second;
  if (f0.IsBoson()) return f0;

  else {
    DO_STACK_TRACE;
    for (std::map<std::string,ATOOLS::Flavour>::const_iterator fit(m_fmap.begin());
	 fit!=m_fmap.end();++fit) PRINT_VAR(fit->first<<" "<<fit->second);
    for (Flavour_Map::const_iterator fit(m_eflmap.begin());
	 fit!=m_eflmap.end();++fit) PRINT_VAR(fit->first<<" "<<fit->second);
    PRINT_VAR(f0<<" "<<id);
    PRINT_VAR(this<<" "<<Name()<<" "<<Demangle(typeid(*this).name()));
    PRINT_VAR(p_mapproc<<" "<<p_mapproc->Name()<<" "<<Demangle(typeid(*p_mapproc).name()));
    PRINT_VAR(((Process_Base*)this)->Parent()<<" "<<((Process_Base*)this)->Parent()->Name());
    PRINT_VAR(p_mapproc->Parent()<<" "<<p_mapproc->Parent()->Name());
    THROW(critical_error,"Flavour map incomplete!");
  }
  return f0;
}

ATOOLS::Flavour AMEGIC::Process_Base::ReMap
(const ATOOLS::Flavour &ifl,const size_t &cid) const
{
  if (Partner()==NULL || Partner()==this) return ifl;
  bool swap(cid&((1<<m_nin)-1));
  Flavour fl(swap?ifl.Bar():ifl);
  static std::map<Flavour,std::string> s_flsmap;
  std::map<Flavour,std::string>::iterator flsit(s_flsmap.find(fl));
  if (flsit==s_flsmap.end()) flsit=s_flsmap.insert(make_pair(fl,ToString(fl))).first;
  static std::map<size_t,std::string> s_idsmap;
  std::map<size_t,std::string>::iterator idsit(s_idsmap.find(cid));
  if (idsit==s_idsmap.end()) idsit=s_idsmap.insert(make_pair(cid,ToString(cid))).first;
  std::string id(flsit->second+idsit->second);
  std::map<std::string,ATOOLS::Flavour>::const_iterator fit(m_fmap.find(id));
  if (fit!=m_fmap.end()) return swap?fit->second.Bar():fit->second;
  else {
    size_t ccid=((1<<(m_nin+m_nout))-1)-cid;
    std::map<size_t,std::string>::iterator idsit(s_idsmap.find(ccid));
    if (idsit==s_idsmap.end()) idsit=s_idsmap.insert(make_pair(ccid,ToString(ccid))).first;
    id=flsit->second+idsit->second;
    std::map<std::string,ATOOLS::Flavour>::const_iterator fit(m_fmap.find(id));
    if (fit!=m_fmap.end()) return swap?fit->second:fit->second.Bar();
    Flavour_Map::const_iterator efit(m_eflmap.find(fl));
    if (efit!=m_eflmap.end()) return swap?efit->second.Bar():efit->second;
    if (ifl.IsBoson()) return ifl;
    else {
      DO_STACK_TRACE;
      PRINT_VAR(this<<" "<<Name()<<" "<<Demangle(typeid(*this).name()));
      PRINT_VAR(p_mapproc<<" "<<p_mapproc->Name()<<" "<<Demangle(typeid(*p_mapproc).name()));
      PRINT_VAR(((Process_Base*)this)->Parent()<<" "<<((Process_Base*)this)->Parent()->Name());
      PRINT_VAR(p_mapproc->Parent()<<" "<<p_mapproc->Parent()->Name());
      THROW(critical_error,"Flavour map incomplete!");
    }
  }
  return ifl;
}

AMEGIC::Process_Base *AMEGIC::Process_Base::GetReal()
{
  return this;
}

bool AMEGIC::Process_Base::FlavCompare(PHASIC::Process_Base *const proc)
{
  if (m_nin!=proc->NIn() || m_nout!=proc->NOut()) return false;
  bool flavsok(true);
  for (size_t i(0);i<m_nin+m_nout;++i)
    if (m_flavs[i].IsAnti()!=proc->Flavours()[i].IsAnti()) flavsok=false;
  return flavsok;
}

std::string  AMEGIC::Process_Base::CreateLibName()
{
  std::string name(m_name);
  size_t apos(name.find(m_pinfo.m_addname));
  if (apos!=std::string::npos)
    name.erase(apos,m_pinfo.m_addname.length());
  size_t bpos(name.find("__QCD("));
  if (bpos!=std::string::npos) {
    size_t epos(name.find(')',bpos));
    if (epos!=std::string::npos)
      name.replace(bpos,epos-bpos+1,"");
  }
  bpos=name.find("__EW(");
  if (bpos!=std::string::npos) {
    size_t epos(name.find(')',bpos));
    if (epos!=std::string::npos)
      name.replace(bpos,epos-bpos+1,"");
  }
  name=ShellName(name+"__O");
  int sep(0);
  for (size_t i(0);i<m_pinfo.m_mincpl.size();++i) {
    name+=ToString(m_pinfo.m_mincpl[i])+"_";
    if (m_pinfo.m_mincpl[i]!=m_pinfo.m_maxcpl[i]) sep=1;
  }
  if (sep) {
    name+="_";
    for (size_t i(0);i<m_pinfo.m_maxcpl.size();++i)
      name+=ToString(m_pinfo.m_maxcpl[i])+"_";
  }
  name.erase(name.length()-1,1);
  msg_Debugging()<<"-> "<<name<<std::endl;
  return name;
}
