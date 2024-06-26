#include "AMEGIC++/Main/Process_Group.H"
#include "AMEGIC++/Main/Single_Process.H"
#include "AMEGIC++/Main/Single_Process_MHV.H"
#include "AMEGIC++/Main/Single_Process_External.H"
#include "AMEGIC++/Main/Single_Process_Combined.H"
#include "AMEGIC++/DipoleSubtraction/Single_Virtual_Correction.H"
#include "AMEGIC++/DipoleSubtraction/Single_Real_Correction.H"
#include "AMEGIC++/Main/Process_Tags.H"

#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Library_Loader.H"

using namespace AMEGIC;
using namespace MODEL;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

AMEGIC::Process_Group::Process_Group()
{
  p_testmoms=NULL;
  p_channellibnames = new std::list<std::string>();
}

AMEGIC::Process_Group::~Process_Group()
{
  if (p_testmoms) delete[] p_testmoms;
  delete p_channellibnames;
}

PHASIC::Process_Base *AMEGIC::Process_Group::GetProcess(const PHASIC::Process_Info &pi) const
{
  int typechk(0);
  if (pi.m_fi.m_nlotype&nlo_type::real) typechk++;
  if (pi.m_fi.m_nlotype&nlo_type::vsub ||
      pi.m_fi.m_nlotype&nlo_type::loop ||
      pi.m_fi.m_nlotype&nlo_type::born) typechk++;
  if (typechk>1 && pi.m_nlomode!=nlo_mode::yfs)
    THROW(fatal_error,"NLO_Parts 'RS' and 'BVI' must be assigned separately!");

  nlo_type::code nlotype=pi.m_fi.m_nlotype;
  // QCD/EW subtraction
  // dont do for YFS!
  // if(pi.m_nlomode!=nlo_mode::yfs){
  if (nlotype!=nlo_type::lo && pi.m_nlomode!=nlo_mode::yfs) {
    if (nlotype&nlo_type::real || nlotype&nlo_type::rsub) {
      Single_Real_Correction *src = new Single_Real_Correction();
      if (!(nlotype&nlo_type::real)) src->SetNoTree(true);
      return src;
    }
    else if (nlotype&nlo_type::born ||
             nlotype&nlo_type::vsub || nlotype&nlo_type::loop) {
      return new Single_Virtual_Correction();
    }
    else {
      THROW(fatal_error, "Unknown NLO type" + ToString(nlotype) + ".");
    }
  }
  // LO processes
  else if (nlotype==nlo_type::lo || nlotype==nlo_type::real || pi.m_nlomode==nlo_mode::yfs) {
    if (pi.m_amegicmhv>0) {
      if (pi.m_amegicmhv==10 ||
          pi.m_amegicmhv==12) return new Single_Process_External();
      if (pi.m_amegicmhv==11) return new Single_Process_Combined();
      if (CF.MHVCalculable(pi)) return new Single_Process_MHV();
      if (pi.m_amegicmhv==2) return NULL;
    }
    return new Single_Process();
  }
  else {
    return NULL;
  }
}

bool AMEGIC::Process_Group::Initialize(PHASIC::Process_Base *const proc)
{
  if (m_whitelist.size()>0) {
    if (m_whitelist.find(proc->Name())==m_whitelist.end()) return false;
  }
  if (!p_testmoms) {
    if (!p_pinfo) p_pinfo=Translate(m_pinfo);
    p_testmoms = new Vec4D[m_nin+m_nout];
    Phase_Space_Handler::TestPoint(p_testmoms,&Info(),Generator());
    PrepareTestMoms(p_testmoms,m_nin,m_nout);
  }
  AMEGIC::Process_Base* apb=proc->Get<AMEGIC::Process_Base>();
  apb->SetPrintGraphs(m_pinfo.m_gpath);
  apb->SetTestMoms(p_testmoms);
  int res=apb->InitAmplitude(p_model,p_top,m_umprocs,m_errprocs); 
  if (s_partcommit)
    My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/",0);
  if (res) proc->SetParent((PHASIC::Process_Base*)this);
  return res;
}

int AMEGIC::Process_Group::InitAmplitude(Amegic_Model * model,Topology * top)
{
  p_model=model;
  p_top=top;
  m_mfname = "P"+ToString(m_nin)+"_"+ToString(m_nout)+"/"+m_name+".map";

  string name = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_mfname;
  string buffer;
  My_In_File file(name);
  file.Open();
  for (;*file;) {
    getline(*file,buffer);
    if (buffer.length()>0) {
      m_whitelist.insert(buffer);
    }
  }
  file.Close();
  if (m_whitelist.size()>0) m_mfname="";

  return true;
}

void AMEGIC::Process_Group::WriteMappingFile()
{
  if (m_mfname==string("")) return;
  std::string name = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/"+m_mfname;

  std::string str, tmp;
  My_In_File in(name);
  if (in.Open())
    for (getline(*in,tmp);in->good(); getline(*in,tmp)) str+=tmp+"\n";
  in.Close();
  My_Out_File out(name);
  out.Open();
  *out<<str;
  for (size_t i=0;i<m_procs.size();i++) *out<<m_procs[i]->Name()<<"\n";
  out.Close();
}

bool AMEGIC::Process_Group::SetUpIntegrator()
{
  if (p_parent==NULL || (*p_parent)[0]->IsGroup()/* this is fudgy, need mode ... */) {
    for (size_t i(0);i<m_procs.size();i++) {
      int res=m_procs[i]->Get<AMEGIC::Process_Base>()->SetUpIntegrator();
      if (s_partcommit)
	My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/",0);
      if (!res) return false;
    }
  }
  if (m_nin==2) {
    if ( (m_flavs[0].Mass() != p_int->ISR()->Flav(0).Mass()) ||
	 (m_flavs[1].Mass() != p_int->ISR()->Flav(1).Mass()) ) 
      p_int->ISR()->SetPartonMasses(m_flavs);
  }
  for (size_t i=0;i<m_procs.size();i++) 
    m_procs[i]->Get<AMEGIC::Process_Base>()->AddChannels(p_channellibnames);
  return true;
}

void AMEGIC::Process_Group::SetPrintGraphs(std::string gpath) 
{
 for (size_t i=0;i<m_procs.size();i++) 
   m_procs[i]->Get<AMEGIC::Process_Base>()->SetPrintGraphs(gpath);
}


#define PTS long unsigned int
#define PT(ARG) (PTS)(ARG)

typedef PHASIC::Single_Channel *(*Lib_Getter_Function)
  (int nin,int nout,ATOOLS::Flavour* fl,
   ATOOLS::Integration_Info * const info,PHASIC::Phase_Space_Handler *psh);

PHASIC::Single_Channel *LoadChannels(int nin,int nout,ATOOLS::Flavour* fl,
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

bool AMEGIC::Process_Group::FillIntegrator
(PHASIC::Phase_Space_Handler *const psh)
{
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->Get<AMEGIC::Process_Base>()->RequestVariables(psh);
  Multi_Channel *mc(psh->FSRIntegrator());
  if (mc==NULL) return true;
  My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");
  if (!SetUpIntegrator()) THROW(fatal_error,"No integrator");
  My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");
  if (p_channellibnames->empty()) return true;
  for (std::list<std::string>::iterator it(p_channellibnames->begin());
       it!=p_channellibnames->end();++it) {
    Single_Channel *sc = LoadChannels(NIn(),NOut(),(Flavour*)&Flavours().front(),
				      *it,&*Integrator()->PSHandler());
    if (sc==0) THROW(critical_error,"PS integration channels not compiled");
    sc->SetName(*it);
    mc->Add(sc);
  }
  return false;
}

void AMEGIC::Process_Group::EndOptimize()
{
  int reset(0);
  for (size_t i(0);i<m_procs.size();++i)
    if (m_procs[i]->Get<AMEGIC::Process_Base>()->EOReset()) reset=1;
  if (reset) p_int->Reset();
}

AMEGIC::Process_Base *AMEGIC::Process_Group::Partner() const  
{ 
  return 0; 
}

Amplitude_Handler *AMEGIC::Process_Group::GetAmplitudeHandler() 
{
  return 0;
} 

Helicity *AMEGIC::Process_Group::GetHelicity() 
{
  return 0;
}

bool AMEGIC::Process_Group::NewLibs() 
{
  for (size_t i(0);i<m_procs.size();++i) 
    if (m_procs[i]->Get<AMEGIC::Amegic_Base>()->NewLibs()) return true;
  return false;
}

std::string AMEGIC::Process_Group::PSLibName() 
{
  return "";
}        

void AMEGIC::Process_Group::Minimize()
{
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->Get<AMEGIC::Amegic_Base>()->Minimize();  
}

void AMEGIC::Process_Group::PrintProcessSummary(int it)
{
 if (it==0) cout<<"============================================"<<endl;
  if (it==1) cout<<"  ------------------------------------------"<<endl;
  if (it==2) cout<<"   - - - - - - - - - - - - - - - - - - - -"<<endl;
  for(int i=0;i<it;i++) std::cout<<"  ";
  std::cout<<Name()<<std::endl;

  for (size_t i=0;i<m_procs.size();++i) m_procs[i]->Get<AMEGIC::Process_Base>()->PrintProcessSummary(it+1);
} 

void AMEGIC::Process_Group::FillAlphaHistogram(ATOOLS::Histogram* histo,double weight)
{
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->Get<AMEGIC::Process_Base>()->FillAlphaHistogram(histo,weight);
}
