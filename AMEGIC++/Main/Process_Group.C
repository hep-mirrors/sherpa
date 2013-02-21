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
#include "ATOOLS/Org/Data_Reader.H"
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
  if (pi.m_fi.m_nloqcdtype&nlo_type::real) typechk++;
  if (pi.m_fi.m_nloqcdtype&nlo_type::vsub||
      pi.m_fi.m_nloqcdtype&nlo_type::loop||
      pi.m_fi.m_nloqcdtype&nlo_type::born) typechk++;    
  if (typechk>1) THROW(fatal_error,"NLO_QCD_Parts 'RS', and 'BVI' must be assigned separately!");

  nlo_type::code nloqcd=pi.m_fi.m_nloqcdtype;
  if (nloqcd&nlo_type::real && nloqcd&nlo_type::rsub) {
    return new Single_Real_Correction();
  }
  else if (nloqcd&nlo_type::vsub || nloqcd&nlo_type::loop) {
    return new Single_Virtual_Correction();
  }
  else if (nloqcd==nlo_type::lo || nloqcd==nlo_type::born || nloqcd==nlo_type::real) {
    if (pi.m_amegicmhv>0) {
      if (pi.m_amegicmhv==10) return new Single_Process_External();
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
    Vec4D sum;
    Poincare lab(Vec4D(sqrt(10.0),0.0,0.0,1.0));
    msg_Debugging()<<"After boost:\n";
    for (size_t i(0);i<m_nin+m_nout;++i) {
      lab.Boost(p_testmoms[i]);
      sum+=i<m_nin?-p_testmoms[i]:p_testmoms[i];
      msg_Debugging()<<"  p["<<i<<"] = "<<p_testmoms[i]<<"\n";
    }
    msg_Debugging()<<"} -> sum = "<<sum<<"\n";
    Poincare rot(Vec4D::ZVEC,Vec4D(sqrt(14.0),1.0,2.0,3.0));
    msg_Debugging()<<"After rotation:\n";
    for (size_t i(0);i<m_nin+m_nout;++i) {
      rot.Rotate(p_testmoms[i]);
      sum+=i<m_nin?-p_testmoms[i]:p_testmoms[i];
      msg_Debugging()<<"  p["<<i<<"] = "<<p_testmoms[i]<<"\n";
    }
    msg_Debugging()<<"} -> sum = "<<sum<<"\n";
  }
  AMEGIC::Process_Base* apb=proc->Get<AMEGIC::Process_Base>();
  apb->SetPrintGraphs(m_pinfo.m_gpath);
  apb->SetTestMoms(p_testmoms);
  if (!apb->InitAmplitude(p_model,p_top,m_umprocs,m_errprocs)) return false;
  proc->SetParent((PHASIC::Process_Base*)this);
  return true;
}

int AMEGIC::Process_Group::InitAmplitude(Model_Base * model,Topology * top)
{
  p_model=model;
  p_top=top;
  m_oew=m_oqcd=0;
  m_mfname = "P"+ToString(m_nin)+"_"+ToString(m_nout)+"/"+m_name+".map";

  string name = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_mfname;
  string buffer;
  ifstream file;
  file.open(name.c_str()); 
  for (;file;) {
    getline(file,buffer); 
    if (buffer.length()>0) {
      m_whitelist.insert(buffer);
    }
  }
  file.close();
  if (m_whitelist.size()>0) m_mfname="";

  return true;
}

void AMEGIC::Process_Group::WriteMappingFile()
{
  if (m_mfname==string("")) return;
  std::string name = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/"+m_mfname;
  fstream file;
  file.open(name.c_str(),ios::out|ios::app); 

  for (size_t i=0;i<m_procs.size();i++) file<<m_procs[i]->Name()<<endl;

  file.close();
}

bool AMEGIC::Process_Group::SetUpIntegrator()
{
  if (p_parent==NULL || (*p_parent)[0]->IsGroup()/* this is fudgy, need mode ... */) {
    for (size_t i(0);i<m_procs.size();i++)
      if (!(m_procs[i]->Get<AMEGIC::Process_Base>()->SetUpIntegrator())) return false;
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


int AMEGIC::Process_Group::PerformTests()
{
  int res(1);
  for (size_t i=0;i<m_procs.size();i++) 
    if (!m_procs[i]->Get<AMEGIC::Amegic_Base>()->PerformTests()) res=0;
  return res;
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
  if (!SetUpIntegrator()) THROW(fatal_error,"No integrator");
  if (p_channellibnames->empty()) return true;
  Multi_Channel *mc(psh->FSRIntegrator());
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
