#include <iostream>
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Git_Info.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/binreloc.h"
#include "ATOOLS/Org/My_File.H"
#include <stdlib.h>
#include <unistd.h>
#include <pwd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#ifdef ARCH_DARWIN
#include <sys/sysctl.h>
#endif
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

std::map<const std::string,const Git_Info*> *
ATOOLS::Git_Info::s_objects=NULL;

bool ATOOLS::Git_Info::s_check=false;

Git_Info::Git_Info(const std::string &name,
		   const std::string &branch,
		   const std::string &revision,
		   const std::string &checksum):
  m_name(name), m_branch(branch), m_revision(revision),
  m_checksum(checksum)
{
  static bool init(false);
  if (!init || s_objects==NULL) {
    s_objects = new std::map<const std::string,const Git_Info*>();
    init=true;
  }
  s_objects->insert(make_pair(name,this));
  if (s_check) {
    std::string branch(s_objects->begin()->second->Branch());
    std::string revision(s_objects->begin()->second->Revision());
    if (m_branch!=branch || m_revision!=revision)
      msg_Info()<<"===> "<<m_name<<" has local modifications "
		<<m_checksum<<" <===\n";
  }
}

Git_Info::~Git_Info()
{
  for (std::map<const std::string,const Git_Info*>::iterator
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
  gen.m_nevents   = 0;
  gen.m_ecms      = gen.m_accu = gen.m_sqrtaccu = 0.;
  gen.m_beam1     = gen.m_beam2      = Flavour(kf_none);
  gen.m_pdfset[0] = gen.m_pdfset[1] = NULL;
  gen.m_ngenevents = 0;
  gen.m_ntrials   = 0;
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
  char *var=NULL;
  gen.m_variables["SHERPASYS"]=std::string(((var=getenv("SHERPASYS"))==NULL?"":var));
  gen.m_variables["SHERPA_CPP_PATH"]=std::string(((var=getenv("SHERPA_CPP_PATH"))==NULL?"":var));
  gen.m_variables["SHERPA_LIB_PATH"]=std::string(((var=getenv("SHERPA_LIB_PATH"))==NULL?"":var));
  gen.m_variables[LD_PATH_NAME]=std::string(((var=getenv(LD_PATH_NAME))==NULL?"":var));
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

void Run_Parameter::RegisterDefaults()
{
  Settings& s = Settings::GetMainSettings();
  s["PRETTY_PRINT"].SetDefault("On");
  s["LOG_FILE"].SetDefault("");
  s["SHERPA_CPP_PATH"].SetDefault("");
  s["SHERPA_LIB_PATH"].SetDefault("");
  s["EVENTS"].SetDefault(100);

  std::vector<long int> seeds = {-1, -1, -1, -1};
  s["RANDOM_SEED"].SetDefault(seeds);
  for (size_t i = 0; i < seeds.size(); ++i)
    s["RANDOM_SEED" + ToString(i + 1)].SetDefault(seeds[i]);

  s["MEMLEAK_WARNING_THRESHOLD"].SetDefault(1<<24);
  s["TIMEOUT"].SetDefault(-1.0);

  std::string logfile{ s["LOG_FILE"].Get<std::string>() };
  s["BATCH_MODE"].SetDefault(logfile==""?1:3);

  s["CITATION_DEPTH"].SetDefault(1);
  s["MPI_SEED_MODE"].SetDefault(0);
  s["MPI_EVENT_MODE"].SetDefault(0);

  s["RLIMIT_BY_CPU"].SetDefault(false);
  s["STACK_TRACE"].SetDefault(1);
  s["NUM_ACCURACY"].SetDefault(1.e-10);
}

void Run_Parameter::Init()
{
  RegisterDefaults();
  Settings& s = Settings::GetMainSettings();

  // set path
  std::string path=s.GetPath();
  if (path[0]!='/') path=gen.m_variables["SHERPA_RUN_PATH"]+"/"+path;
  while (path.length()>0
         && (path.back()=='/' || path.back()=='.')) {
    path.pop_back();
  }

  gen.m_timer.Start();

  // set user name
  struct passwd* user_info = getpwuid(getuid());
  if (!user_info) gen.m_username="<unknown user>";
  else gen.m_username=user_info->pw_gecos;
  size_t pos(gen.m_username.find(','));
  if (pos<std::string::npos)
    gen.m_username.erase(pos,gen.m_username.length()-pos);
  char hn[32];
  if (gethostname(hn,32)) gen.m_hostname="<unknown host>";
  else gen.m_hostname=std::string(hn);

  // initialise output
  std::string color = s["PRETTY_PRINT"].Get<std::string>();
  if (color=="Off") msg->SetModifiable(false);
  msg->Init();

  // print welcome message
  if (msg->LevelIsInfo())  {
    msg_Out() << "Welcome to Sherpa, " << gen.m_username << " on "
              << gen.m_hostname
              << ".\nInitialization of framework underway ...\n";
  }
  msg_Info()
    <<"Local time: "<<rpa->gen.Timer().TimeString(0)<<std::endl;

  // set cpp path
  std::string cpppath=s["SHERPA_CPP_PATH"].Get<std::string>();
  if (cpppath.length()==0 || cpppath[0]!='/') {
    if (path!=gen.m_variables["SHERPA_RUN_PATH"]) gen.m_variables["SHERPA_CPP_PATH"]=path;
    else if (gen.m_variables["SHERPA_CPP_PATH"].length()==0)
      gen.m_variables["SHERPA_CPP_PATH"]=gen.m_variables["SHERPA_RUN_PATH"];
  }
  if (cpppath.length()) gen.m_variables["SHERPA_CPP_PATH"]+=(cpppath[0]=='/'?"":"/")+cpppath;

  // set lib path
  std::string libpath=s["SHERPA_LIB_PATH"].Get<std::string>();
  if (libpath.length()>0 && libpath[0]=='/') gen.m_variables["SHERPA_LIB_PATH"]=libpath;
  else if (gen.m_variables["SHERPA_LIB_PATH"].length()==0)
    gen.m_variables["SHERPA_LIB_PATH"]=gen.m_variables["SHERPA_CPP_PATH"]
      +std::string("/Process/Amegic/lib");

  msg_Tracking()
    <<METHOD<<"(): Paths are {\n"
    <<"   SHERPA_INC_PATH = "  <<gen.m_variables["SHERPA_INC_PATH"]  <<"\n"
    <<"   SHERPA_SHARE_PATH = "<<gen.m_variables["SHERPA_SHARE_PATH"]<<"\n"
    <<"   SHERPA_CPP_PATH = "  <<gen.m_variables["SHERPA_CPP_PATH"]  <<"\n"
    <<"   SHERPA_LIB_PATH = "  <<gen.m_variables["SHERPA_LIB_PATH"]  <<"\n"
    <<"}"<<std::endl;


  // configure event generation
  gen.m_variables["EVENT_GENERATION_MODE"]="-1";
  gen.m_nevents = s["EVENTS"].Get<long int>();

  s_loader->AddPath(rpa->gen.Variable("SHERPA_RUN_PATH"));

  // read only if defined (no error message if not defined)
  long int seed;
  std::vector<long int> seeds = s["RANDOM_SEED"].GetVector<long int>();
  for (int i(0);i<4;++i) gen.m_seeds[i] = -1;
  for (int i(0);i<Min((int)seeds.size(),4);++i) gen.m_seeds[i] = seeds[i];
  for (int i(0);i<4;++i) {
    seed = s["RANDOM_SEED" + ToString(i + 1)].Get<long int>();
    if (seed != -1) {
      gen.m_seeds[i] = seed;
    }
  }
  int nseed=0;
  for (int i(0);i<4;++i) if (gen.m_seeds[i]>0) ++nseed;
  if (nseed==0) {
    gen.m_seeds[0]=1234;
  }
  else if (nseed>1) {
    if (gen.m_seeds[0]<0) gen.m_seeds[0]=12345;
    if (gen.m_seeds[1]<0) gen.m_seeds[1]=65435;
    if (gen.m_seeds[2]<0) gen.m_seeds[2]=34221;
    if (gen.m_seeds[3]<0) gen.m_seeds[3]=12345;
  }

#ifdef USING__MPI
  int rank=mpi->Rank();
  int size=mpi->Size();
  if (s["MPI_EVENT_MODE"].Get<int>()==1) {
    gen.m_nevents = (gen.m_nevents%size == 0) ? (gen.m_nevents/size) : (gen.m_nevents/size+1);
  }
  if (s["MPI_SEED_MODE"].Get<int>()==0) {
    msg_Info()<<"Seed mode: '*'\n";
    for (int i(0);i<4;++i)
      if (gen.m_seeds[i]>0) gen.m_seeds[i]*=rank+1;
  }
  else {
    msg_Info()<<"Seed mode: '+'\n";
    for (int i(0);i<4;++i)
      if (gen.m_seeds[i]>0) gen.m_seeds[i]+=rank;
  }
#endif

  std::string seedstr;
  if (gen.m_seeds[1]>0)
    for (int i(1);i<4;++i) seedstr+="_"+ToString(gen.m_seeds[i]);
  gen.SetVariable("RNG_SEED",ToString(gen.m_seeds[0])+seedstr);
  gen.SetVariable("MEMLEAK_WARNING_THRESHOLD",
		  ToString(s["MEMLEAK_WARNING_THRESHOLD"].Get<int>()));
  gen.m_timeout = s["TIMEOUT"].Get<double>();
  if (gen.m_timeout<0.) gen.m_timeout=0.;
  rpa->gen.m_timer.Start();
  gen.m_batchmode = s["BATCH_MODE"].Get<int>();
  gen.m_clevel= s["CITATION_DEPTH"].Get<int>();
  int ncpus(getncpu());
  msg_Tracking()<<METHOD<<"(): Getting number of CPU cores: "
		<<ncpus<<"."<<std::endl;
#ifdef RLIMIT_AS
  rlimit lims;
  getrlimit(RLIMIT_AS,&lims);
  double slim(getpmem());
  msg_Tracking()<<METHOD<<"(): Getting memory limit: "
		<<slim/double(1<<30)<<" GB."<<std::endl;
  const auto rawrlim = s["RLIMIT_AS"]
    .SetDefault((rlim_t)(slim-double(100*(1<<20))))
    .Get<double>();
  lims.rlim_cur = rlim_t(rawrlim < 1.0 ? slim * rawrlim : rawrlim);
  if (s["RLIMIT_BY_CPU"].Get<bool>()) lims.rlim_cur/=(double)ncpus;
  if (setrlimit(RLIMIT_AS, &lims) != 0) {
    msg_Error() << om::brown << om::bold << "WARNING" << om::reset
                << ": Memory limit can not be set." << std::endl;
  }
  getrlimit(RLIMIT_AS,&lims);
  msg_Info() << "Memory limit: " << lims.rlim_cur / double(1 << 30) << " GB"
             << std::endl;
#endif
  gen.m_accu = s["NUM_ACCURACY"].Get<double>();
  gen.m_sqrtaccu = sqrt(gen.m_accu);
  if (gen.m_seeds[1]>0) {
    ran->SetSeed(gen.m_seeds[0],gen.m_seeds[1],gen.m_seeds[2],gen.m_seeds[3]);
  }
  else { ran->SetSeed(gen.m_seeds[0]); }
  msg_Debugging()<<METHOD<<"(): Global tags {\n";
  const String_Map &gtags{ s.GetTags() };
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
    if (m_timeout==0.)
      return true;
    else
      return m_timer.UserTime()<m_timeout;
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
  if (Citations().empty()) return;
  if (!Settings::GetMainSettings()["WRITE_REFERENCES_FILE"].Get<bool>())
    return;
  std::string refname("References.tex");
  My_Out_File f((rpa->gen.Variable("SHERPA_RUN_PATH")+"/"+refname).c_str());
  f.Open();
  *f<<"%% This is a citation summary file generated by Sherpa "
    <<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION<<"\n";
  *f<<"%% on "+rpa->gen.Timer().TimeString(0)<<".\n";
  *f<<"%% It contains LaTeX-style citations for Sherpa and any external\n";
  *f<<"%% scientific software or physics results used to generate the given result with\n";
  *f<<"%% Sherpa. Upload this file to https://inspirehep.net/bibliography-generator to\n";
  *f<<"%% generate a bibliography for your publication.\n";
  *f<<"\n\\documentclass{article}\n\n\\begin{document}\n"<<std::endl;
  for (size_t i=0; i<Citations().size(); ++i) {
    *f<<Citations()[i]<<std::endl;
  }
  *f<<"\n\\end{document}\n\n"<<std::endl;
  *f<<"%% You have used the following configuration:\n";
  PrintGitVersion(*f,1,"%% ");
#ifdef USING__MPI
  if (mpi->Rank()==0) {
#endif
    const int framewidth {60};
    msg_Out() << Frame_Header{framewidth};
    msg_Out() << Frame_Line{MyStrStream{}
                                << om::bold << om::green
                                << "Please cite the publications listed in '"
                                << refname << "'." << om::reset,
                            framewidth};
    msg_Out() << Frame_Footer{framewidth};
#ifdef USING__MPI
  }
#endif
}

void  Run_Parameter::Gen::SetEcms(double _ecms) {
  m_ecms = _ecms;
}

void  Run_Parameter::Gen::SetPBeam(short unsigned int i,Vec4D pbeam) {
  m_pbeam[i]=pbeam;
}
void  Run_Parameter::Gen::SetPBunch(short unsigned int i,Vec4D pbunch) {
  m_pbunch[i]=pbunch;
}
void  Run_Parameter::Gen::SetBeam1(const Flavour b) {
  m_beam1  = b;
}
void  Run_Parameter::Gen::SetBeam2(const Flavour b) {
  m_beam2  = b;
}

std::string Run_Parameter::Gen::Variable(const std::string &key)
{
  const auto it = m_variables.find(key);
  if (it == m_variables.end()) {
    THROW(fatal_error,
          "Runtime parameter \"" + key
          + "\" not registered and no default given to fall back on");
  }
  return it->second;
}

void Run_Parameter::Gen::PrintGitVersion(std::ostream &str,
                                         const bool& shouldprintversioninfo,
					 const std::string &prefix)
{
  const std::map<const std::string,const Git_Info*> &info(*Git_Info::Infos());
  if (info.empty()) THROW(fatal_error,"No Git information");
  std::string branch(info.begin()->second->Branch());
  std::string revision(info.begin()->second->Revision());
  if (branch.find("rel-")!=0)
    msg_Info() << om::bold << om::brown << "WARNING" << om::reset
               << ": You are using an unsupported development branch.\n";
  str<<prefix<<"Git branch: "<<branch<<"\n"<<prefix<<"Revision: "<<revision;
  if (shouldprintversioninfo) str<<" {\n";
  else str<<"."<<std::endl;
  for (std::map<const std::string,const Git_Info*>::const_iterator
	 iit(info.begin());iit!=info.end();++iit) {
    if (shouldprintversioninfo) str<<prefix<<"  "<<iit->second->Checksum()
                                   <<"  "<<iit->second->Name()<<"\n";
    if (iit->second->Revision()!=revision) str<<prefix
      <<"  -> "<<iit->second->Name()<<" has local modifications: "
      <<iit->second->Checksum()<<"\n";
  }
  if (shouldprintversioninfo) str<<prefix<<"}\n";
  Git_Info::SetCheck(true);
}
