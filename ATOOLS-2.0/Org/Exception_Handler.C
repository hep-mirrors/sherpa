#include "Exception.H"

#include "Run_Parameter.H"
#include "Scaling.H"
#include <sys/types.h>
#include <unistd.h>
#include <execinfo.h>
#include <dlfcn.h>

#define MAX_BACKTRACE_DEPTH 128

using namespace ATOOLS;

bool Exception_Handler::s_active=true;
bool Exception_Handler::s_prepared=false;
bool Exception_Handler::s_stacktrace=true;

unsigned int Exception_Handler::s_exitcode=0;
Exception *Exception_Handler::s_exception=0;

unsigned int Exception_Handler::s_nbus=0;
unsigned int Exception_Handler::s_nsegv=0;

std::vector<Exception_Handler::Tester_Function> 
Exception_Handler::s_testerfunctions=std::vector<Exception_Handler::Tester_Function>();
std::vector<Exception_Handler::Terminator_Function> 
Exception_Handler::s_terminatorfunctions=std::vector<Exception_Handler::Terminator_Function>();
std::vector<Tester_Object*> 
Exception_Handler::s_testerobjects=std::vector<Tester_Object*>();
std::vector<Terminator_Object*> 
Exception_Handler::s_terminatorobjects=std::vector<Terminator_Object*>();

bool Exception_Handler::ApproveTerminate()
{
  msg.Error()<<"Exception_Handler::ApproveTerminate(): Asking for termination ..."<<std::endl;
  if (s_testerfunctions.size()==0 && s_testerobjects.size()==0) {
      msg.Error()<<"... approved."<<std::endl;
      return true;
  }
  if (s_testerfunctions.size()>0) {
    for (size_t i=0;i<s_testerobjects.size();++i) if (s_testerfunctions[i]()) {
      msg.Error()<<"... approved."<<std::endl;
      return true;
    }
  }
  if (s_testerobjects.size()>0) {
    for (size_t i=0;i<s_testerobjects.size();++i) if (s_testerobjects[i]->ApproveTerminate()) {
      msg.Error()<<"... approved."<<std::endl;
      return true;
    }
  }
  msg.Error()<<"... refused."<<std::endl;
  return false;
}

void Exception_Handler::PrepareTerminate()
{
  msg.Error()<<"Exception_Handler::PrepareTerminate(): Preparing termination ..."<<std::endl;
  for (size_t i=0;i<s_terminatorobjects.size();++i) s_terminatorobjects[i]->PrepareTerminate(); 
  for (size_t i=0;i<s_terminatorfunctions.size();++i) s_terminatorfunctions[i](); 
  msg.Error()<<"... prepared."<<std::endl;
}

void Exception_Handler::Exit(int exitcode)
{
  msg.Error()<<om::bold<<"Exception_Handler::Exit: "<<om::reset<<om::blue
	     <<"Exiting Sherpa with code "<<om::reset<<om::bold<<"("<<om::red<<exitcode
	     <<om::reset<<om::bold<<")"<<om::reset<<std::endl;
  msg.LogFile()<<"Exception_Handler::Exit: Exiting Sherpa with code ("
	       <<exitcode<<")"<<std::endl;
  exit(exitcode);
}

void Exception_Handler::Terminate() 
{
  GenerateStackTrace(msg.LogFile(),true,"! ");
  if (s_stacktrace) GenerateStackTrace(msg.Error());
  if (!ApproveTerminate()) {
    s_exception=NULL;
    return;
  }
  PrepareTerminate();
  s_prepared=true;
  if (!s_active) abort();
  SetExitCode();
  Exit(s_exitcode);
}

void Exception_Handler::RemoveTesterObject(Tester_Object *const testerobject)
{
  for (std::vector<Tester_Object*>::iterator toit=s_testerobjects.begin();
       toit!=s_testerobjects.end();) {
    if (*toit==testerobject) toit=s_testerobjects.erase(toit); 
    else ++toit;
  }
}

void Exception_Handler::RemoveTerminatorObject(Terminator_Object *const terminatorobject)
{
  for (std::vector<Terminator_Object*>::iterator toit=s_terminatorobjects.begin();
       toit!=s_terminatorobjects.end();) {
    if (*toit==terminatorobject) toit=s_terminatorobjects.erase(toit); 
    else ++toit;
  }
}

void Exception_Handler::SetExitCode()
{
  if (s_exception==NULL) return;
  if (s_exception->m_class=="ISR_Handler")                 s_exitcode=151;
  else if (s_exception->m_class=="MI_Base")                s_exitcode=211;
  else if (s_exception->m_class=="Simple_Chain")           s_exitcode=212;
  else if (s_exception->m_class=="Matrix_Element_Handler") s_exitcode=201;
  else s_exitcode=1;
}

void Exception_Handler::SignalHandler(int signal) 
{
  std::string input="y";
  msg.Error()<<om::bold<<"Exception_Handler::SignalHandler: "<<om::reset<<om::blue
	     <<"Signal "<<om::reset<<om::bold<<"("<<om::red<<signal
	     <<om::reset<<om::bold<<")"<<om::reset<<om::blue<<" caught. "<<om::reset<<std::endl;
  switch (signal) {
  case SIGSEGV:
    ++s_nsegv;
    if (!rpa.gen.BatchMode()) {
      msg.Error()<<"   Do you want to debug the program (y/n)? "<<om::reset;
      std::cin>>input;
      if (input=="y" || input=="Y") {
	system((std::string("gdb Sherpa ")+ATOOLS::ToString(getpid())).c_str());
      }
    }
    if (s_nsegv>3) {
      msg.Error()<<om::reset<<"   Abort immediately."<<om::reset<<std::endl;
      kill(getpid(),9);
    }
  case SIGABRT:
    if (!s_active && s_prepared) abort();
  case SIGTERM:
  case SIGXCPU:
    msg.Error()<<om::reset<<"   Cannot continue."<<om::reset<<std::endl;
    s_exitcode=2;
    Terminate();
    break;
  case SIGINT:
    if (!rpa.gen.BatchMode()) {
      msg.Error()<<"   Do you want to stop the program (y/n/s/d)? "<<om::reset;
      std::cin>>input;
    }
    if (input=="d" || input=="D") {
      system((std::string("gdb Sherpa ")+ATOOLS::ToString(getpid())).c_str());
    }
    else if (input=="s" || input=="S") {
      GenerateStackTrace(std::cout,false);
      std::cout<<" Hit <return> to continue. ";      
      std::cin.get();
      std::cin.get();
    }
    if (input!="y" && input!="Y") return;
    s_exitcode=1;
    Terminate();
    break;
  case SIGBUS:
    ++s_nbus;
    if (s_nbus>3) {
      msg.Error()<<om::reset<<"   Abort immediately."<<om::reset<<std::endl;
      kill(getpid(),9);
    }
    msg.Error()<<om::reset<<"   Cannot continue."<<om::reset<<std::endl;
    s_exitcode=3;
    Terminate();
    break;
  case SIGFPE:
    msg.Error()<<"   Sherpa does not throw floating point exceptions."<<om::reset<<std::endl;
    break;
  default:
    msg.Error()<<"   Cannot handle signal."<<om::reset<<std::endl;
    s_exitcode=1;
    Terminate();
  }
}

void Exception_Handler::GenerateStackTrace(std::ostream &ostr,const bool endline,
					   const std::string &comment)
{
  ostr<<comment<<om::bold<<"Exception_Handler::GenerateStackTrace(..): "
      <<om::reset<<om::blue<<"Generating stack trace "<<om::reset
      <<"(adapted from ROOT version 3.10) "<<om::bold<<"{"<<om::reset<<std::endl;
  // adapted from root version 3.10 TUnixSystem.cxx
  void *trace[MAX_BACKTRACE_DEPTH];
  int depth=backtrace(trace,MAX_BACKTRACE_DEPTH);
  for (int n=0; n<depth;++n) {
    unsigned long addr=(unsigned long)trace[n];
    Dl_info info;
    if (dladdr(trace[n],&info) && info.dli_fname && info.dli_fname[0]) {
      unsigned long symaddr=(unsigned long)info.dli_saddr;
      const char *symname=info.dli_sname;
      if (!info.dli_sname || !info.dli_sname[0]) symname="<unknown>";
      ostr<<comment<<"   0x"<<std::hex<<symaddr<<std::dec<<" in "
	  <<symname<<" from "<<info.dli_fname<<std::endl;
    } 
    else {
      ostr<<comment<<"   "<<addr<<" in <unknown function>"<<std::endl;
    }
  }
  ostr<<comment<<om::bold<<"}"<<om::reset;
  if (endline) ostr<<std::endl;
}
