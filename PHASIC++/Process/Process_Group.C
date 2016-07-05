#include "PHASIC++/Process/Process_Group.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

Process_Group::~Process_Group()
{
  Clear();
}

size_t Process_Group::Size() const
{
  return m_procs.size();
}

Process_Base *Process_Group::operator[](const size_t &i)
{
  return m_procs[i];
}

Weight_Info *Process_Group::OneEvent(const int wmode,const int mode) 
{
  p_selected=NULL;
  if (p_int->TotalXS()==0.0) {
    p_selected=m_procs[int(ATOOLS::ran->Get()*m_procs.size())];
    return p_selected->OneEvent(mode);
  }
  double disc=p_int->SelectionWeight(wmode)*ATOOLS::ran->Get();
  for (size_t i=0;i<m_procs.size();++i) {
    disc-=dabs(m_procs[i]->Integrator()->SelectionWeight(wmode));
    if (disc<=0.) {
      p_selected=m_procs[i];
      return p_selected->OneEvent(mode);
    }
  }
  msg_Error()<<METHOD<<"(): Cannot select any process. xs = "
	     <<p_int->TotalXS()*rpa->Picobarn()<<" pb."<<std::endl;
  return NULL;
}

double Process_Group::Differential(const Vec4D_Vector &p)
{
  m_lastb=m_last=0.0;
  for (size_t i(0);i<m_procs.size();++i) {
    m_last+=m_procs[i]->Differential(p);
    m_lastb+=m_procs[i]->LastB();
  }
  if (IsNan(m_last))
    msg_Error()<<METHOD<<"(): "<<om::red
		<<"Cross section is 'nan'."<<om::reset<<std::endl;
  return m_last;
}

void Process_Group::SetScale(const Scale_Setter_Arguments &args)
{
  for (size_t i(0);i<m_procs.size();++i) m_procs[i]->SetScale(args);
}
  
void Process_Group::SetKFactor(const KFactor_Setter_Arguments &args)
{
  for (size_t i(0);i<m_procs.size();++i) m_procs[i]->SetKFactor(args);
}

bool Process_Group::IsGroup() const
{
  return true;
}

void Process_Group::Add(Process_Base *const proc) 
{
  if (proc==NULL) return;
  std::string name(proc->Name()), add(proc->Info().m_addname);
  if (add.length() && name.rfind(add)!=std::string::npos)
    name.erase(name.rfind(add),add.length());
  if (m_procmap.find(name)!=m_procmap.end())
    THROW(critical_error,"Doubled process '"+name+"'");
  m_procmap[name]=proc;
  if (m_maxcpl.size()<proc->MaxOrders().size()) {
    m_maxcpl.resize(proc->MaxOrders().size(),0);
    m_mincpl.resize(proc->MinOrders().size(),99);
  }
  for (size_t i(0);i<m_maxcpl.size();++i) {
    m_maxcpl[i]=Max(m_maxcpl[i],proc->MaxOrder(i));
    m_mincpl[i]=Min(m_mincpl[i],proc->MinOrder(i));
  }
  if (m_nin>0 && m_nout>0 &&
      (m_nin!=proc->NIn() || m_nout!=proc->NOut())) {
    msg_Error()<<METHOD<<"(): Cannot add process '"
	       <<proc->Name()<<"' to group '"<<m_name<<"'.\n"
	       <<"  Inconsistent number of external legs."<<std::endl; 
    return;
  }  
  m_procs.push_back(proc);
}

bool Process_Group::Remove(Process_Base *const proc) 
{
  for (std::vector<Process_Base*>::iterator xsit=m_procs.begin();
       xsit!=m_procs.end();++xsit) {
    if (*xsit==proc) {
      m_procs.erase(xsit);
      return true;
    }
  }
  return false;
}

bool Process_Group::Delete(Process_Base *const proc) 
{
  if (Remove(proc)) {
    delete proc;
    return true;
  }
  return false;
}

void Process_Group::Clear() 
{
  while (m_procs.size()>0) {
    delete m_procs.back();
    m_procs.pop_back();
  }
}

Process_Base *Process_Group::GetProcess(const std::string &name)
{
  std::map<std::string,Process_Base*>::const_iterator 
    pit(m_procmap.find(name));
  if (pit!=m_procmap.end()) return pit->second;
  if (name==m_name) return this;
  for (size_t i(0);i<m_procs.size();++i)
    if (m_procs[i]->IsGroup()) {
      Process_Base *proc(m_procs[i]->Get<Process_Group>()->GetProcess(name));
      if (proc) return proc;
    }
  return NULL;
}

bool Process_Group::CalculateTotalXSec(const std::string &resultpath,
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
    psh->InitCuts();
    for (size_t i=0;i<m_procs.size();++i)
      m_procs[i]->BuildCuts(psh->Cuts());
    p_int->ISR()->SetSprimeMin(psh->Cuts()->Smin());
  }
  psh->CreateIntegrators();
  p_int->SetResultPath(resultpath);
  p_int->ReadResults();
  p_int->SetTotal(0);
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
  if (p_int->TotalXS()!=0.0) {
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

void Process_Group::SetLookUp(const bool lookup)
{
  m_lookup=lookup;
  for (size_t i(0);i<m_procs.size();++i) m_procs[i]->SetLookUp(lookup);
}

bool Process_Group::CheckFlavours
(const Subprocess_Info &ii,const Subprocess_Info &fi,int mode) const
{
  int charge(0), strong(0);
  size_t quarks(0), nin(ii.m_ps.size()), nout(fi.m_ps.size());
  for (size_t i(0);i<nin;++i) {
    const Flavour &fl(ii.m_ps[i].m_fl);
    charge+=-fl.IntCharge();
    if (abs(fl.StrongCharge())!=8)
      strong+=-fl.StrongCharge();
    quarks+=fl.IsQuark();
    if (mode==0 && quarks>m_pinfo.m_nmaxq) {
      msg_Debugging()<<METHOD<<"(): '"<<GenerateName(ii,fi)<<"': n_q > "
		     <<m_pinfo.m_nmaxq<<". Skip process.\n";
      return false;
    }
  }
  for (size_t i(0);i<nout;++i) {
    const Flavour &fl(fi.m_ps[i].m_fl);
    charge+=fl.IntCharge();
    if (abs(fl.StrongCharge())!=8)
      strong+=fl.StrongCharge();
    quarks+=fl.IsQuark();
    if (mode==0 && quarks>m_pinfo.m_nmaxq) {
      msg_Debugging()<<METHOD<<"(): '"<<GenerateName(ii,fi)<<"': n_q > "
		     <<m_pinfo.m_nmaxq<<". Skip process.\n";
      return false;
    }
  }
  if (mode==0 && quarks<m_pinfo.m_nminq) {
    msg_Debugging()<<METHOD<<"(): '"<<GenerateName(ii,fi)<<"': n_q < "
		   <<m_pinfo.m_nminq<<". Skip process.\n";
    return false;
  }
  if (charge!=0 || strong!=0) return false;
  bool res(true);
  for (size_t i(0);i<fi.m_ps.size();++i) {
    if (fi.m_ps[i].m_ps.empty()) continue;
    Subprocess_Info cii;
    cii.m_ps.push_back(fi.m_ps[i]);
    if (!CheckFlavours(cii,fi.m_ps[i],1)) {
      res=false;
      break;
    }
  }
  return res;
}

void Process_Group::SetFlavour(Subprocess_Info &cii,Subprocess_Info &cfi,
			       const ATOOLS::Flavour &fl,const size_t i) const
{
  if (i<m_nin) cii.SetExternal(fl,i);
  else cfi.SetExternal(fl,i-m_nin);
}

bool Process_Group::ConstructProcess(Process_Info &pi)
{
  if (!CheckFlavours(pi.m_ii,pi.m_fi)) return false;
  Process_Info cpi(pi);
  SortFlavours(cpi);
  std::string name(GenerateName(cpi.m_ii,cpi.m_fi));
  if (m_procmap.find(name)!=m_procmap.end()) return false;
  Process_Base *proc(GetProcess(cpi));
  if (!proc) return false;
  proc->SetGenerator(Generator());
  proc->Init(pi,p_int->Beam(),p_int->ISR());
  if (!Initialize(proc)) {
    msg_Debugging()<<METHOD<<"(): Init failed for '"
		   <<proc->Name()<<"'\n";
    delete proc;
    m_procmap[name]=NULL;
    return false;
  }
  msg_Debugging()<<METHOD<<"(): Init ok '"
		 <<proc->Name()<<"'\n";
  Add(proc);
  return true;
}

bool Process_Group::ConstructProcesses(Process_Info &pi,const size_t &ci)
{
  if (ci==m_nin+m_nout) {
    if (!ConstructProcess(pi)) return false;
#ifdef USING__MPI
    if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
    std::string mapfile(rpa->gen.Variable("SHERPA_CPP_PATH")
			+"/Process/Sherpa/"+m_name+".map");
    std::string str, tmp;
    My_In_File in(mapfile);
    if (in.Open())
      for (getline(*in,tmp);in->good();
	   getline(*in,tmp)) str+=tmp+"\n";
    in.Close();
    My_Out_File out(mapfile);
    if (!out.Open()) THROW(fatal_error,"Cannot open '"+mapfile+"'");
    *out<<str;
    Flavour_Vector fl(m_procs.back()->Info().m_ii.GetExternal());
    for (size_t i(0);i<fl.size();++i) *out<<(long int)fl[i]<<" ";
    fl=m_procs.back()->Info().m_fi.GetExternal();
    for (size_t i(0);i<fl.size();++i) *out<<(long int)fl[i]<<" ";
    *out<<"0\n";
    out.Close();
#ifdef USING__MPI
    }
#endif
    return true;
  }
  bool one(false);
  for (size_t j(0);j<m_flavs[ci].Size();++j) {
    SetFlavour(pi.m_ii,pi.m_fi,m_flavs[ci][j],ci);
    if (ConstructProcesses(pi,ci+1)) one=true;
  }
  return one;
}

bool Process_Group::ConstructProcesses()
{
  Process_Info cpi(m_pinfo);
  std::string mapfile(rpa->gen.Variable("SHERPA_CPP_PATH")
		      +"/Process/Sherpa/"+m_name+".map");
  msg_Debugging()<<"checking for '"<<mapfile<<"' ... "<<std::flush;
  if (FileExists(mapfile)) {
    msg_Debugging()<<"found"<<std::endl;
    My_In_File map(mapfile);
    if (!map.Open()) THROW(fatal_error,"Corrupted map file '"+mapfile+"'");
    long int cfl, cnt;
    *map>>cfl;
    while (!map->eof()) {
      for (cnt=0;cnt<m_nin+m_nout && !map->eof();++cnt) {
        SetFlavour(cpi.m_ii,cpi.m_fi,Flavour(abs(cfl),cfl<0),cnt);
        *map>>cfl;
      }
      if (cnt!=m_nin+m_nout || cfl || !ConstructProcess(cpi))
	THROW(fatal_error,"Corrupted map file '"+mapfile+"'");
      *map>>cfl;
    }
    return m_procs.size();
  }
  msg_Debugging()<<"not found"<<std::endl;
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()==0)
#endif
  My_In_File::ExecDB
    (rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/","begin");
  bool res(ConstructProcesses(cpi,0));
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()==0)
#endif
  My_In_File::ExecDB
    (rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/","commit");
  return res;
}

void Process_Group::SetGenerator(ME_Generator_Base *const gen) 
{ 
  Process_Base::SetGenerator(gen);
  for (size_t i(0);i<m_procs.size();++i) 
    m_procs[i]->SetGenerator(gen);
}

void Process_Group::SetShower(PDF::Shower_Base *const ps) 
{ 
  Process_Base::SetShower(ps);
  for (size_t i(0);i<m_procs.size();++i) 
    m_procs[i]->SetShower(ps);
}

void Process_Group::SetVariationWeights(SHERPA::Variation_Weights *const vw)
{
  Process_Base::SetVariationWeights(vw);
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->SetVariationWeights(vw);
}

void Process_Group::SetSelector(const Selector_Key &key)
{
  Process_Base::SetSelector(key);
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->SetSelector(key);
}

void Process_Group::SetFixedScale(const std::vector<double> &s)
{
  Process_Base::SetFixedScale(s);
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->SetFixedScale(s);
}

void Process_Group::SetSelectorOn(const bool on)
{
  Process_Base::SetSelectorOn(on);
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->SetSelectorOn(on);
}

void Process_Group::SetUseBIWeight(bool on)
{
  Process_Base::SetUseBIWeight(on);
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->SetUseBIWeight(on);
}

bool Process_Group::Trigger(const Vec4D_Vector &p)
{
  bool trigger=false;
  for (size_t i(0);i<Size();++i)
    if (m_procs[i]->Trigger(p)) trigger=true;
  return trigger;
}
 
void Process_Group::FillOnshellConditions()
{
  Process_Base::FillOnshellConditions();
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->FillOnshellConditions();
}
