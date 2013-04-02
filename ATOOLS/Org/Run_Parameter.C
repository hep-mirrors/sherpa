#include <iostream> 
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Data_Writer.H"
#include "ATOOLS/Org/SVN_Info.H"
#include "ATOOLS/Org/binreloc.h"
#include <stdlib.h>
#include <unistd.h>
#include <pwd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <limits>

namespace ATOOLS {
  Run_Parameter *rpa(NULL);
}

double getpmem()
{
#if defined(ARCH_LINUX) || defined(ARCH_UNIX)
  unsigned long int ps(getpagesize());
  unsigned long int sc(sysconf(_SC_PHYS_PAGES));
  return double(sc)*double(ps);
#endif
#ifdef ARCH_DARWIN
  int mib[2]={CTL_HW,HW_PHYSMEM};
  unsigned long int miblen(2);
  unsigned long int pmem(0), len(sizeof(pmem));
  if (sysctl(mib,miblen,&pmem,&len,NULL,0)!=0) {
    std::cerr<<"sysctl failed"<<std::endl;
    return 0.0;
  }
  return double(pmem);
#endif
  std::cerr<<"cannot determine physical memory"<<std::endl;
  return 1.e15;
}

int getncpu()
{
#if defined(ARCH_LINUX) || defined(ARCH_UNIX)
  return sysconf(_SC_NPROCESSORS_ONLN);
#endif
#ifdef ARCH_DARWIN
  int mib[2]={CTL_HW,HW_AVAILCPU};
  unsigned long int miblen(2);
  unsigned long int ncpu(1), len(sizeof(ncpu));
  if (sysctl(mib,miblen,&ncpu,&len,NULL,0)!=0) {
    mib[1]=HW_NCPU;
    if (sysctl(mib,miblen,&ncpu,&len,NULL,0)!=0) {
      std::cerr<<"sysctl failed"<<std::endl;
      ncpu = 1;
    }
  }
  return ncpu;
#endif
  std::cerr<<"cannot determine number of cpus"<<std::endl;
  return 1;
}

using namespace ATOOLS;
using namespace std;

std::map<const std::string,const SVN_Info*> *
ATOOLS::SVN_Info::s_objects=NULL;

SVN_Info::SVN_Info(const std::string &name,
		   const std::string &branch,
		   const std::string &revision,
		   const std::string &checksum):
  m_name(name), m_branch(branch), m_revision(revision),
  m_checksum(checksum)
{
  static bool init(false);
  if (!init || s_objects==NULL) {
    s_objects = new std::map<const std::string,const SVN_Info*>();
    init=true;
  }
  s_objects->insert(make_pair(name,this));
}

SVN_Info::~SVN_Info()
{
  for (std::map<const std::string,const SVN_Info*>::iterator
	 it(s_objects->begin());it!=s_objects->end();++it)
    if (it->second==this) {
      s_objects->erase(it);
      break;
    }
  if (s_objects->empty()) delete s_objects;
}

Run_Parameter::Run_Parameter() 
{
  AnalyseEnvironment();
  gen.m_analysis   = 0;
  gen.m_nevents   = 0;
  gen.m_cutscheme = 0;
  gen.m_ecms      = gen.m_accu        = 0.;
  gen.m_beam1     = gen.m_beam2      = Flavour(kf_none);
  gen.m_pdfset[0] = gen.m_pdfset[1] = NULL;
  gen.m_ngenevents = 0;
  gen.m_batchmode = 1;
  gen.SetTimeOut(3600);
  gen.m_softsc = 0;
  gen.m_hardsc = 0;
  gen.m_pbeam[0] = Vec4D(0.,0.,0.,0.);
  gen.m_pbeam[1] = Vec4D(0.,0.,0.,0.);
  gen.m_clevel=100;
} 

std::ostream &ATOOLS::operator<<(std::ostream &str,const Run_Parameter &rp)
{ 
  return str<<"("<<&rp<<"): {\n}"; 
}

void Run_Parameter::AnalyseEnvironment() 
{
#ifdef __GNUC__
#if __GNUC__ == 2 && __GNUC_MINOR__ == 96
#error Sherpa was not designed for gcc 2.96
#endif
#endif
  char *var=NULL;
  gen.m_variables["SHERPASYS"]=std::string(((var=getenv("SHERPASYS"))==NULL?"":var));
  gen.m_variables["SHERPA_CPP_PATH"]=std::string(((var=getenv("SHERPA_CPP_PATH"))==NULL?"":var));
  gen.m_variables["SHERPA_LIB_PATH"]=std::string(((var=getenv("SHERPA_LIB_PATH"))==NULL?"":var));
  gen.m_variables["SHERPA_DAT_PATH"]=std::string(((var=getenv("SHERPA_DAT_PATH"))==NULL?"":var));
  gen.m_variables["LD_LIBRARY_PATH"]=std::string(((var=getenv("LD_LIBRARY_PATH"))==NULL?"":var));
  gen.m_variables["SHERPA_RUN_PATH"]=GetCWD();
  gen.m_variables["HOME"]=std::string(((var=getenv("HOME"))==
				       NULL?gen.m_variables["SHERPA_RUN_PATH"]:var));

  // The paths are determined with the following fallback route:
  // 1. Environment variable
  // 2. binreloc (if enabled during configure)
  // 3. Hard coded value in installation directory
  // set share path
  string sharepath=SHERPA_SHARE_PATH;
  string includepath=SHERPA_INCLUDE_PATH;
  string librarypath=SHERPA_LIBRARY_PATH;
  BrInitError error;
  if (br_init_lib(&error)) {
    string BR_prefix=br_find_prefix(SHERPA_PREFIX);
    sharepath=BR_prefix+"/share/SHERPA-MC";
    includepath=BR_prefix+"/include/SHERPA-MC";
    librarypath=BR_prefix+"/lib/SHERPA-MC";
  }
  
  gen.m_variables["SHERPA_SHARE_PATH"]=
    (var=getenv("SHERPA_SHARE_PATH"))==NULL?sharepath:var;

  // set include path
  gen.m_variables["SHERPA_INC_PATH"]=
    (var=getenv("SHERPA_INCLUDE_PATH"))==NULL?includepath:var;

  // set library path 
  gen.m_variables["SHERPA_LIBRARY_PATH"]=
    (var=getenv("SHERPA_LIBRARY_PATH"))==NULL?librarypath:var;

}

void Run_Parameter::Init(std::string path,std::string file,int argc,char* argv[])
{
  m_path = path;
  path=gen.m_variables["PATH_PIECE"];
  gen.m_timer.Start();
  struct passwd* user_info = getpwuid(getuid());
  if (!user_info) gen.m_username="<unknown user>";
  else gen.m_username=user_info->pw_gecos;
  size_t pos(gen.m_username.find(','));
  if (pos<std::string::npos)
    gen.m_username.erase(pos,gen.m_username.length()-pos);
  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(m_path);
  dr.SetInputFile(file);
  std::string color=dr.GetValue<std::string>("PRETTY_PRINT","On");
  if (color=="Off") msg->SetModifiable(false);
  std::string outputlevel = dr.GetValue<std::string>("OUTPUT","2");
  std::string logfile = dr.GetValue<std::string>("LOG_FILE","");
  msg->Init(outputlevel,logfile);
  msg->SetMPIMode(dr.GetValue<int>("MPI_OUTPUT",0));
  if (msg->LevelIsInfo()) 
    msg_Out()<<"Welcome to "<<exh->ProgramName()<<", "<<gen.m_username
	     <<". Initialization of framework underway."<<std::endl;
  msg_Info()<<"The local time is "<<rpa->gen.Timer().TimeString(0)<<"."<<std::endl;
  // make path nice
  if (path.length()>0) {
    if (path[0]!='/') path=gen.m_variables["SHERPA_RUN_PATH"]+"/"+path;
    while (path.length()>0 && 
	   (path[path.length()-1]=='/' || path[path.length()-1]=='.')) 
      path=path.substr(0,path.length()-1);
  }

  // set cpp path
  std::string cpppath=dr.GetValue<std::string>("SHERPA_CPP_PATH",std::string(""));
  if (cpppath.length()>0 && cpppath[0]=='/') gen.m_variables["SHERPA_CPP_PATH"]=cpppath;
  else if (path!=gen.m_variables["SHERPA_RUN_PATH"]) gen.m_variables["SHERPA_CPP_PATH"]=path;
  else if (gen.m_variables["SHERPA_CPP_PATH"].length()==0) 
    gen.m_variables["SHERPA_CPP_PATH"]=gen.m_variables["SHERPA_RUN_PATH"];
  // set lib path
  std::string libpath=dr.GetValue<std::string>("SHERPA_LIB_PATH",std::string(""));
  if (libpath.length()>0 && libpath[0]=='/') gen.m_variables["SHERPA_LIB_PATH"]=libpath;
  else if (gen.m_variables["SHERPA_LIB_PATH"].length()==0) 
    gen.m_variables["SHERPA_LIB_PATH"]=gen.m_variables["SHERPA_CPP_PATH"]
      +std::string("/Process/lib");
  if (gen.m_variables["SHERPA_DAT_PATH"].length()==0) {
    if (path.length()>0 && path[0]=='/') gen.m_variables["SHERPA_DAT_PATH"]=path;
    else gen.m_variables["SHERPA_DAT_PATH"]=gen.m_variables["SHERPA_RUN_PATH"]+"/"+path;
  }
  msg_Tracking()<<METHOD<<"(): Paths are {\n"
		<<"   SHERPA_INC_PATH = "<<gen.m_variables["SHERPA_INC_PATH"]<<"\n"
		<<"   SHERPA_SHARE_PATH = "<<gen.m_variables["SHERPA_SHARE_PATH"]<<"\n"
		<<"   SHERPA_CPP_PATH = "<<gen.m_variables["SHERPA_CPP_PATH"]<<"\n"
		<<"   SHERPA_LIB_PATH = "<<gen.m_variables["SHERPA_LIB_PATH"]<<"\n"
		<<"   SHERPA_DAT_PATH = "<<gen.m_variables["SHERPA_DAT_PATH"]<<"\n"
		<<"}"<<std::endl;
#ifndef __sgi
  setenv("LD_LIBRARY_PATH",(gen.m_variables["LD_LIBRARY_PATH"]+std::string(":")+
			    gen.m_variables["SHERPA_LIB_PATH"]).c_str(),1);
#endif
  MakeDir(gen.m_variables["HOME"]+"/.sherpa/",true);
  gen.m_analysis           = dr.GetValue<int>("ANALYSIS",0);
  dr.SetAllowUnits(true);
  gen.m_nevents            = dr.GetValue<long int>("EVENTS",100);
  dr.SetAllowUnits(false);
  s_loader->AddPath(rpa->gen.Variable("SHERPA_RUN_PATH"));

  // read only if defined (no error message if not defined)
  std::vector<long int> seeds;
  std::vector<long int> seed1;
  std::vector<long int> seed2;
  gen.m_seed2 = -1;
  if (dr.VectorFromFile(seeds,"RANDOM_SEED")) {
    gen.m_seed = seeds[0];
    // if 2nd seed is given, store it
    if (seeds.size() == 2) { gen.m_seed2 = seeds[1]; } 
  } 
  else if (dr.VectorFromFile(seed1,"RANDOM_SEED1") && 
	   dr.VectorFromFile(seed2,"RANDOM_SEED2")) {
    gen.m_seed  = seed1[0];
    gen.m_seed2 = seed2[0]; 
  }
  else gen.m_seed=1234;

#ifdef USING__MPI
  int rank=MPI::COMM_WORLD.Get_rank();
  gen.m_seed*=rank+1;
  if (gen.m_seed2>0) gen.m_seed2*=rank+1;
#endif

  if (gen.m_seed2<0) gen.SetVariable("RNG_SEED",ToString(gen.m_seed));
  else gen.SetVariable("RNG_SEED",ToString(gen.m_seed)+" "+ToString(gen.m_seed2));

  gen.m_timeout = dr.GetValue<double>("TIMEOUT",std::numeric_limits<double>::max());
  if (gen.m_timeout<0.) gen.m_timeout=0.;
  rpa->gen.m_timer.Start();
  gen.m_batchmode = dr.GetValue<int>("BATCH_MODE",logfile==""?1:3);
  gen.m_clevel= dr.GetValue<int>("CITATION_DEPTH",1);
  int ncpus(getncpu());
  msg_Tracking()<<METHOD<<"(): Getting number of CPU cores: "
		<<ncpus<<"."<<std::endl;
  gen.SetVariable("NUMBER_OF_CPUS",ToString(ncpus));
#ifdef RLIMIT_AS
  rlimit lims;
  getrlimit(RLIMIT_AS,&lims);
  double slim(getpmem());
#ifdef USING__LHAPDF
  if (slim+400000000.0 < double((1<<32)-1)) {
    slim+=400000000.0;
  }
#endif
  msg_Tracking()<<METHOD<<"(): Getting memory limit: "
		<<slim/double(1<<30)<<" GB."<<std::endl;
  std::vector<std::string> aspars;
  if (!dr.VectorFromFile(aspars,"RLIMIT_AS")) {
    lims.rlim_cur=(rlim_t)(slim-double(100*(1<<20)));
  }
  else {
    if (aspars.size()==1) {
      lims.rlim_cur=(rlim_t)
	(ToType<double>(aspars[0])*slim);
    }
    else if (aspars.size()==2) {
      if (aspars[1]=="MB") lims.rlim_cur=(rlim_t)
	(ToType<double>(aspars[0])*(1<<20));
      else if (aspars[1]=="GB") lims.rlim_cur=(rlim_t)
	(ToType<double>(aspars[0])*(1<<30));
      else if (aspars[1]=="%") lims.rlim_cur=(rlim_t)
	(ToType<double>(aspars[0])*slim/100.0);
      else
	THROW(fatal_error,"Invalid syntax in '"+m_file+"'");
    }
    else {
      THROW(fatal_error,"Invalid syntax in '"+m_file+"'");
    }
  }
  int limpercpu=dr.GetValue<int>("RLIMIT_BY_CPU",0);
  if (limpercpu) lims.rlim_cur/=(double)ncpus;
  if (setrlimit(RLIMIT_AS,&lims)!=0)
    msg_Error()<<METHOD<<"(): Cannot set memory limit."<<std::endl;
  getrlimit(RLIMIT_AS,&lims);
  msg_Info()<<METHOD<<"(): Setting memory limit to "
	    <<lims.rlim_cur/double(1<<30)<<" GB."<<std::endl;
#endif
  int stacktrace = dr.GetValue<int>("STACK_TRACE",1);
  exh->SetStackTrace(stacktrace);
  gen.m_accu = dr.GetValue<double>
    ("Num._Accuracy",dr.GetValue<double>("NUM_ACCURACY",1.e-10));
  //gen.m_runtime            = dr.GetValue<std::string>("Runtime"); // Time
  if (gen.m_seed2!=-1) { ran->SetSeed(gen.m_seed, gen.m_seed2); }
                  else { ran->SetSeed(gen.m_seed); }
  msg_Debugging()<<METHOD<<"(): Set global tags {\n";
  const String_Map &gtags(Read_Write_Base::GlobalTags());
  for (String_Map::const_iterator tit(gtags.begin());tit!=gtags.end();++tit)
    msg_Debugging()<<"  '"<<tit->first<<"' -> '"<<tit->second<<"'\n";
  msg_Debugging()<<"}\n";
}

Run_Parameter::~Run_Parameter() 
{ 
  if (msg->Level()>=1) gen.m_timer.PrintTime();
}

bool Run_Parameter::Gen::CheckTime(const double limit)
{ 
  if (limit==0.) {
    if (m_timeout>0.) return m_timer.UserTime()<m_timeout;
  }
  else {
    return m_timer.UserTime()<limit;
  }
  return false;
}

void Run_Parameter::Gen::AddCitation(const size_t &level,
                                     const std::string &cite)
{
  if (level<=m_clevel) {
    for (size_t i=0; i<m_cites.size(); ++i) if (m_cites[i]==cite) return;
    m_cites.push_back(cite);
  }
}

void Run_Parameter::Gen::WriteCitationInfo()
{
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()) return;
#endif
  if (Citations().empty()) return;
  std::string refname("Sherpa_References.tex");
  std::ofstream f((rpa->gen.Variable("SHERPA_RUN_PATH")+"/"+refname).c_str());
  f<<"%% Citation summary file generated by Sherpa "
   <<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION<<std::endl;
  f<<"%% PID "+ToString(getpid())+" on "+rpa->gen.Timer().TimeString(0)<<std::endl;
  f<<"\n\\documentclass{article}\n\n\\begin{document}\n"<<std::endl;
  for (size_t i=0; i<Citations().size(); ++i) {
    f<<Citations()[i]<<std::endl;
  }
  f<<"\n\\end{document}\n\n"<<std::endl;
  f<<"%% You have used the following configuration:\n";
  PrintSVNVersion(f,1,"%% ");
  std::cout<<std::string(72,'-')<<"\n"
	   <<om::bold<<"Please cite the publications listed in '"
	   <<om::red<<refname<<om::reset<<om::bold<<"'."<<om::reset
	   <<"\n  Extract the bibtex list by running 'get_bibtex "
	   <<refname<<"'\n  or email the file to 'slaclib2@slac.stanford.edu'"
	   <<", subject 'generate'.\n"<<std::string(72,'-')<<std::endl;
}

void  Run_Parameter::Gen::SetEcms(double _ecms)     { 
  m_ecms    = _ecms;
}
void  Run_Parameter::Gen::SetCplScale(double _cplscale)     { 
  m_cplscale = _cplscale;
}
void  Run_Parameter::Gen::SetPBeam(short unsigned int i,Vec4D pbeam) { 
  m_pbeam[i]=pbeam;
}
void  Run_Parameter::Gen::SetBeam1(const Flavour b) { 
  m_beam1  = b;   
}
void  Run_Parameter::Gen::SetBeam2(const Flavour b) { 
  m_beam2  = b;   
}

std::string Run_Parameter::Gen::Variable(const std::string &key,const std::string &def) 
{ 
  return m_variables.find(key)!=m_variables.end()?m_variables[key]:def; 
}

void Run_Parameter::Gen::PrintSVNVersion(std::ostream &str,const int mode,
					 const std::string &prefix)
{
  const std::map<const std::string,const SVN_Info*> &info(*SVN_Info::Infos());
  if (info.empty()) THROW(fatal_error,"No SVN information");
  std::string branch(info.begin()->second->Branch());
  std::string revision(info.begin()->second->Revision());
  str<<prefix<<"SVN branch "<<branch<<", revision "<<revision;
  if (mode&1) str<<" {\n";
  else str<<"."<<std::endl;
  for (std::map<const std::string,const SVN_Info*>::const_iterator
	 iit(info.begin());iit!=info.end();++iit) {
    if (mode&1) str<<prefix<<" "<<iit->second->Checksum()
		   <<"  "<<iit->second->Name()<<"\n";
    if (iit->second->Branch()!=branch) str<<prefix
      <<"===> "<<iit->second->Name()<<" has branch "<<iit->second->Branch()
      <<", first seen was "<<branch<<" <===\n";
    if (iit->second->Revision()!=revision) str<<prefix
      <<"===> "<<iit->second->Name()<<" has revision "<<iit->second->Revision()
      <<", first seen was "<<revision<<" <===\n";
  }
  if (mode&1) str<<prefix<<"}\n";
  str<<std::endl;
}
