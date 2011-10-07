#include "ATOOLS/Org/Exception.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include <sys/types.h>
#include <unistd.h>

#define USING_Stack_Trace
#ifndef __USE_GNU
#ifdef __GNUC__
#define __USE_GNU
#ifdef ARCH_DARWIN
#undef USING_Stack_Trace
#endif
#else 
#undef USING_Stack_Trace
#endif
#endif

#ifdef USING_Stack_Trace
#include <execinfo.h>
#include <dlfcn.h>
#define MAX_BACKTRACE_DEPTH 128
#endif

#ifdef USING__MPI
#include "mpi.h"
#endif

using namespace ATOOLS;

ATOOLS::Exception_Handler *ATOOLS::exh(NULL);

Exception_Handler::Exception_Handler():
  m_active(true), m_prepared(false), m_stacktrace(true), 
  m_print(true), m_noremove(false),
  m_signal(0), m_exitcode(0), m_exception(0),
  m_nbus(0), m_nsegv(0), m_mpi_timeout(3600),
  m_progname("Sherpa")
{
  std::set_terminate(ATOOLS::Terminate);
  std::set_unexpected(ATOOLS::Terminate);
  signal(SIGSEGV,ATOOLS::SignalHandler);
  signal(SIGINT,ATOOLS::SignalHandler);
  signal(SIGPIPE,ATOOLS::SignalHandler);
  signal(SIGBUS,ATOOLS::SignalHandler);
  signal(SIGFPE,ATOOLS::SignalHandler);
  signal(SIGABRT,ATOOLS::SignalHandler);
  signal(SIGTERM,ATOOLS::SignalHandler);
  signal(SIGXCPU,ATOOLS::SignalHandler);
  signal(SIGUSR1,ATOOLS::SignalHandler);
#ifdef USING__MPI
  m_mpi.resize(MPI::COMM_WORLD.Get_size(),1);
#endif
}

Exception_Handler::~Exception_Handler()
{
}

bool Exception_Handler::ReadInStatus(const std::string &path)
{
  bool success(true);
  msg_Info()<<METHOD<<"(): Reading status from '"<<path<<"' {"<<std::endl;
  for (size_t i=0;i<m_terminatorobjects.size();++i) 
    if (!m_terminatorobjects[i]->ReadInStatus(path)) success=false;
  msg_Info()<<"}"<<std::endl;
  return success;
}

bool Exception_Handler::ApproveTerminate()
{
  static size_t inttrials=0;
  if (++inttrials>2) kill(getpid(),9);
  if (m_print) msg_Tracking()<<"Exception_Handler::ApproveTerminate(): "
			     <<"Asking for termination ..."<<std::endl;
  if (m_testerfunctions.size()==0 && m_testerobjects.size()==0) {
    if (m_print) msg_Tracking()<<"... approved."<<std::endl;
    return true;
  }
  bool approved=true;
  m_noremove=true;
  for (size_t i=0;i<m_testerfunctions.size();++i) 
    if (!m_testerfunctions[i]()) approved=false;
  for (size_t i=0;i<m_testerobjects.size();++i) 
    if (!m_testerobjects[i]->ApproveTerminate()) approved=false;
  m_noremove=false;
  if (approved && m_print) {
    msg_Tracking()<<"... approved."<<std::endl;
  }
  else {
    if (m_print) msg_Tracking()<<"... refused."<<std::endl;
  }
  return approved;
}

void Exception_Handler::PrepareTerminate()
{
  static size_t trials=0;
  if (++trials>3) kill(getpid(),9);
  if (m_print) msg_Tracking()<<"Exception_Handler::PrepareTerminate(): "
			     <<"Preparing termination ..."<<std::endl;
  while (m_terminatorobjects.size()>0) {
    m_noremove=true;
    m_terminatorobjects.back()->PrepareTerminate();
    m_noremove=false;
    std::vector<Terminator_Object*>::iterator end=m_terminatorobjects.end();
    RemoveTerminatorObject(*--end);
  }
  while (m_terminatorfunctions.size()>0) {
    m_noremove=true;
    m_terminatorfunctions.back()();
    m_noremove=false;
    std::vector<Terminator_Function>::iterator end=m_terminatorfunctions.end();
    RemoveTerminatorFunction(*--end);
  }
  if (m_print) msg_Tracking()<<"... prepared."<<std::endl;
}

void Exception_Handler::Exit(int exitcode)
{
  rpa->gen.WriteCitationInfo();
  if (m_print) msg_Error()<<om::bold<<"Exception_Handler::Exit: "
			  <<om::reset<<om::blue<<"Exiting "
			  <<m_progname<<" with code "
			  <<om::reset<<om::bold<<"("
			  <<om::red<<exitcode<<om::reset<<om::bold<<")"
			  <<om::reset<<tm::curon<<std::endl;
  MPISuspend(1);
  exit(exitcode);
}

void Exception_Handler::MPISuspend(int mode)
{
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  int rank=MPI::COMM_WORLD.Get_rank();
  if (rank>0) {
    if (mode) {
      msg_Error()<<METHOD<<"(): MPI rank "<<rank<<", pid "<<getpid()<<" on "
		 <<rpa->gen.Variable("HOSTNAME")<<" is killed."<<std::endl;
    }
    if (m_mpi.size()>rank) {
      int flag=-1;
      MPI_Request req;
      MPI_Irsend(&flag,1,MPI::INT,0,10*size+rank,MPI::COMM_WORLD,&req);
      MPI::Finalize();
    }
  }
  else {
    if (mode) {
      msg_Error()<<METHOD<<"(): MPI master is killed. Abort."<<std::endl;
    }
  }
#endif
}

void Exception_Handler::MPISync()
{
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  int rank=MPI::COMM_WORLD.Get_rank();
  if (rank==0) {
    for (int tag=1;tag<size;++tag) {
      if (m_mpi[tag]==0) continue;
      int flag, test;
      unsigned int time=0;
      for (;time<m_mpi_timeout;++time) {
	MPI_Request req;
	MPI_Irecv(&flag,1,MPI::INT,tag,10*size+tag,MPI::COMM_WORLD,&req);
	MPI_Test(&req,&test,MPI_STATUS_IGNORE);
	if (test) {
	  m_mpi[tag]=flag<0?0:1;
	  if (m_mpi[tag]==0)
	    msg_Error()<<METHOD<<"(): MPI rank "<<tag
		       <<" of "<<size<<" exited."<<std::endl;
	  break;
	}
	MPI_Cancel(&req);
	m_mpi[tag]=0;
	sleep(1);
      }
      if (time==m_mpi_timeout) {
      	msg_Error()<<METHOD<<"(): MPI rank "<<tag
      		   <<" of "<<size<<" does not respond."<<std::endl;
      }
    }
  }
  else {
    int flag=0;
    MPI_Request req;
    MPI_Isend(&flag,1,MPI::INT,0,10*size+rank,MPI::COMM_WORLD,&req);
  }
#endif
}

int Exception_Handler::MPIStat(int rank)
{
  if (rank<0 || rank>m_mpi.size())
    THROW(fatal_error,"Index out of range");
  return m_mpi[rank];
}

void Exception_Handler::Reset()
{
  m_exception=NULL;
  m_signal=0;
}

void ATOOLS::Terminate() 
{
  exh->Terminate();
  exh->Reset();
}

void Exception_Handler::Terminate() 
{
  SetExitCode();
  if ((m_signal!=SIGTERM && m_signal!=SIGINT &&
       m_signal!=SIGXCPU && m_signal!=SIGPIPE) &&
      (m_exception==NULL || 
       (m_exception->Type()!=ex::normal_exit &&
	m_exception->Type()!=ex::missing_input &&
	m_exception->Type()!=ex::missing_module))) {
    if (m_print) {
      if (m_stacktrace) GenerateStackTrace(msg->Error());
    }
    rpa->gen.SetVariable
      ("SHERPA_STATUS_PATH",rpa->gen.Variable("SHERPA_RUN_PATH")+
       "/Status__"+rpa->gen.Timer().TimeString(3));
    msg_Error()<<METHOD<<"(): Pre-crash status saved to '"
	       <<rpa->gen.Variable("SHERPA_STATUS_PATH")<<"'."<<std::endl;
    MakeDir(rpa->gen.Variable("SHERPA_STATUS_PATH"));
  }
  if (!ApproveTerminate()) {
    m_exception=NULL;
    return;
  }
  PrepareTerminate();
  m_prepared=true;
  if (!m_active) abort();
  Exit(m_exitcode);
}

void Exception_Handler::AddTesterFunction(bool (*function)(void))
{ 
  m_testerfunctions.push_back(function); 
}

void Exception_Handler::AddTerminatorFunction(void (*function)(void))
{ 
  m_terminatorfunctions.push_back(function); 
}

void Exception_Handler::AddTesterObject(Tester_Object *const object)
{ 
  m_testerobjects.push_back(object); 
}

void Exception_Handler::AddTerminatorObject(Terminator_Object *const object)
{ 
  m_terminatorobjects.push_back(object); 
}

void Exception_Handler::RemoveTesterObject(Tester_Object *const testerobject)
{
  if (m_noremove) return;
  for (std::vector<Tester_Object*>::iterator toit=m_testerobjects.begin();
       toit!=m_testerobjects.end();) {
    if (*toit==testerobject) toit=m_testerobjects.erase(toit); 
    else ++toit;
  }
}

void Exception_Handler::
RemoveTerminatorObject(Terminator_Object *const terminatorobject)
{
  if (m_noremove) return;
  for (std::vector<Terminator_Object*>::iterator 
	 toit=m_terminatorobjects.begin();
       toit!=m_terminatorobjects.end();) {
    if (*toit==terminatorobject) toit=m_terminatorobjects.erase(toit); 
    else ++toit;
  }
}

void Exception_Handler::RemoveTesterFunction(bool (*function)(void))
{
  if (m_noremove) return;
  for (std::vector<Tester_Function>::iterator tfit=m_testerfunctions.begin();
       tfit!=m_testerfunctions.end();) {
    if (*tfit==function) tfit=m_testerfunctions.erase(tfit); 
    else ++tfit;
  }
}

void Exception_Handler::RemoveTerminatorFunction(void (*function)(void))
{
  if (m_noremove) return;
  for (std::vector<Terminator_Function>::iterator 
	 tfit=m_terminatorfunctions.begin();
       tfit!=m_terminatorfunctions.end();) {
    if (*tfit==function) tfit=m_terminatorfunctions.erase(tfit); 
    else ++tfit;
  }
}

void Exception_Handler::SetExitCode()
{
  m_print=true;
  if (m_exception==NULL) return;
  if (m_exception->m_class=="Amegic") m_exitcode=201;
  else m_exitcode=1;
  if (m_exception->m_type==ex::normal_exit ||
      m_exception->m_type==ex::missing_input) m_print=false;
}

void ATOOLS::SignalHandler(int signal) 
{
  exh->SignalHandler(signal);
  exh->Reset();
}

void Exception_Handler::SignalHandler(int signal) 
{
  m_signal=signal;
  m_print=true;
  msg_Error()<<std::endl<<om::bold<<"Exception_Handler::SignalHandler: "
	     <<om::reset<<om::blue<<"Signal "<<om::reset<<om::bold
	     <<"("<<om::red<<signal<<om::reset<<om::bold<<")"
	     <<om::reset<<om::blue<<" caught. "<<om::reset<<std::endl;
  switch (signal) {
  case SIGSEGV:
    ++m_nsegv;
    GenerateStackTrace(std::cout,false);
    if (m_nsegv>3) {
      msg_Error()<<om::reset<<"   Abort immediately."<<om::reset<<std::endl;
      kill(getpid(),9);
    }
  case SIGABRT:
    if (!m_active && m_prepared) abort();
  case SIGTERM:
  case SIGXCPU:
    msg_Error()<<om::reset<<"   Cannot continue."<<om::reset<<std::endl;
    m_exitcode=2;
    Terminate();
    break;
  case SIGINT:
    m_exitcode=1;
    Terminate();
    break;
  case SIGBUS:
    ++m_nbus;
    if (m_nbus>3) {
      msg_Error()<<om::reset<<"   Abort immediately."<<om::reset<<std::endl;
      kill(getpid(),9);
    }
    GenerateStackTrace(std::cout,false);
    msg_Error()<<om::reset<<"   Cannot continue."<<om::reset<<std::endl;
    m_exitcode=3;
    Terminate();
    break;
  case SIGFPE:
    msg_Error()<<"   Floating point exception."<<om::reset<<std::endl;
    m_exitcode=1;
    Terminate();
    break;
  case SIGPIPE:
    msg_Error()<<"   Pipe closed. Will stop writing."<<om::reset<<std::endl;
    m_exitcode=0;
    Terminate();
    break;
  case SIGUSR1:
#ifdef USING__MPI
    if (MPI::COMM_WORLD.Get_rank()) {
      m_exitcode=1;
      Terminate();
    }
    else {
      msg_Error()<<"   SIGUSR1 ignored by MPI master."<<om::reset<<std::endl;
    }
#endif
    break;
  default:
    msg_Error()<<"   Cannot handle signal."<<om::reset<<std::endl;
    m_exitcode=1;
    Terminate();
  }
}

void Exception_Handler::GenerateStackTrace(std::ostream &ostr,
					   const bool endline,
					   const std::string &comment)
{
#ifdef USING_Stack_Trace
  ostr<<comment<<om::bold<<"Exception_Handler::GenerateStackTrace(..): "
      <<om::reset<<om::blue<<"Generating stack trace "
      <<om::reset<<om::bold<<"\n{"<<om::reset<<std::endl;
  // adapted from root version 3.10 TUnixSystem.cxx
  void *trace[MAX_BACKTRACE_DEPTH];
  int depth=backtrace(trace,MAX_BACKTRACE_DEPTH);
  for (int n=0; n<depth;++n) {
    unsigned long addr=(unsigned long)trace[n];
    Dl_info info;
    if (dladdr(trace[n],&info) && info.dli_fname && info.dli_fname[0]) {
      unsigned long symaddr=(unsigned long)info.dli_saddr;
      if (symaddr==(unsigned long)NULL) continue;
      const char *symname=info.dli_sname;
      if (!info.dli_sname || !info.dli_sname[0]) symname="<unknown function>";
      if (!msg->LevelIsDebugging()) {
	if (std::string(symname).find
	    ("Exception_Handler")!=std::string::npos ||
	    std::string(symname).find
	    ("SignalHandler")!=std::string::npos) continue;
      }
      std::string linfo;
      unsigned long libaddr=(unsigned long)info.dli_fbase;
      unsigned long offset=(addr>=libaddr)?addr-libaddr:libaddr-libaddr;
      char cmd[4096];
      sprintf(cmd,"addr2line -se %s 0x%016lx 2>/dev/null",
	      info.dli_fname,offset);
      if (FILE *pf=popen(cmd,"r")) {
	char buf[2048];
	if (fgets(buf,2048,pf)) {
	  linfo=buf;
	  linfo=linfo.substr(0,linfo.length()-1);
	}
	if (linfo=="??:0") {
	  pclose(pf);
	  sprintf(cmd,"addr2line -se %s 0x%016lx 2>/dev/null",
		  info.dli_fname,addr);
	  pf=popen(cmd,"r");
	  if (fgets(buf,2048,pf)) {
	    linfo=buf;
	    linfo=linfo.substr(0,linfo.length()-1);
	  }
	  if (linfo=="??:0") linfo="";
	}
	pclose(pf);
      }
      ostr<<comment<<"  "<<std::setiosflags(std::ios::left)
	  <<std::setw(15)<<trace[n]<<std::dec
	  <<" in '"<<om::red<<Demangle(symname)<<om::reset<<"' ";
      if (linfo!="") ostr<<"("<<om::lblue<<linfo<<om::reset<<")";
      ostr<<"\n";
      if (msg->LevelIsDebugging()) ostr<<"                from '"<<
	om::brown<<info.dli_fname<<om::reset<<"'\n";
      ostr<<std::flush;
      if (std::string(info.dli_sname)=="main") break;
    } 
    else {
      ostr<<comment<<"   "<<addr<<" in   <unknown function>"<<std::endl;
    }
  }
  ostr<<comment<<om::bold<<"}"<<om::reset;
  if (endline) ostr<<std::endl;
#endif
}
