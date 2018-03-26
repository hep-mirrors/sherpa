#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Terminator_Objects.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_MPI.H"

#include <sys/types.h>
#include <signal.h>
#include <unistd.h>

using namespace ATOOLS;

ATOOLS::Terminator_Object_Handler *ATOOLS::exh(NULL);

void ATOOLS::Terminate() 
{
  exh->Terminate(1);
}

Terminator_Object_Handler::Terminator_Object_Handler():
  m_noremove(false),
  m_nbus(0), m_nsegv(0)
{
  std::set_terminate(ATOOLS::Terminate);
  std::set_unexpected(ATOOLS::Terminate);
  signal(SIGSEGV,ATOOLS::HandleSignal);
  signal(SIGINT,ATOOLS::HandleSignal);
  signal(SIGPIPE,ATOOLS::HandleSignal);
  signal(SIGBUS,ATOOLS::HandleSignal);
  signal(SIGFPE,ATOOLS::HandleSignal);
  signal(SIGABRT,ATOOLS::HandleSignal);
  signal(SIGTERM,ATOOLS::HandleSignal);
  signal(SIGXCPU,ATOOLS::HandleSignal);
  signal(SIGUSR1,ATOOLS::HandleSignal);
}

bool Terminator_Object_Handler::ReadInStatus(const std::string &path)
{
  bool success(true);
  msg_Info()<<METHOD<<"(): Reading status from '"<<path<<"' {"<<std::endl;
  for (size_t i=0;i<m_terminatorobjects.size();++i) 
    if (!m_terminatorobjects[i]->ReadInStatus(path)) success=false;
  msg_Info()<<"}"<<std::endl;
  return success;
}

void Terminator_Object_Handler::PrepareTerminate()
{
  static size_t trials=0;
  if (++trials>3) Abort(1);
  msg_Tracking()<<"Terminator_Object_Handler::PrepareTerminate(): "
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
  msg_Tracking()<<"... prepared."<<std::endl;
}

void Terminator_Object_Handler::Terminate(unsigned int excode) 
{
  PrepareTerminate();
  rpa->gen.WriteCitationInfo();
  msg_Error()<<om::bold<<"Terminator_Object_Handler::Exit: "
	     <<om::reset<<om::blue<<"Exiting Sherpa with code "
	     <<om::reset<<om::bold<<"("
	     <<om::red<<excode<<om::reset<<om::bold<<")"
	     <<om::reset<<tm::curon<<std::endl;
  exit(excode);
}

void Terminator_Object_Handler::AddTerminatorFunction(void (*function)(void))
{ 
  m_terminatorfunctions.push_back(function); 
}

void Terminator_Object_Handler::AddTerminatorObject(Terminator_Object *const object)
{ 
  m_terminatorobjects.push_back(object); 
}

void Terminator_Object_Handler::
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

void Terminator_Object_Handler::RemoveTerminatorFunction(void (*function)(void))
{
  if (m_noremove) return;
  for (std::vector<Terminator_Function>::iterator 
	 tfit=m_terminatorfunctions.begin();
       tfit!=m_terminatorfunctions.end();) {
    if (*tfit==function) tfit=m_terminatorfunctions.erase(tfit); 
    else ++tfit;
  }
}

Terminator_Object_Handler::~Terminator_Object_Handler()
{
}

void ATOOLS::HandleSignal(int signal) 
{
  exh->HandleSignal(signal);
}

void Terminator_Object_Handler::HandleSignal(int signal) 
{
  msg_Error()<<std::endl<<om::bold<<"Terminator_Object_Handler::HandleSignal: "
	     <<om::reset<<om::blue<<"Signal "<<om::reset<<om::bold
	     <<"("<<om::red<<signal<<om::reset<<om::bold<<")"
	     <<om::reset<<om::blue<<" caught. "<<om::reset<<std::endl;


  if (signal!=SIGINT)
    {
      GenerateStackTrace(msg->Error());
      rpa->gen.SetVariable
	("SHERPA_STATUS_PATH",rpa->gen.Variable("SHERPA_RUN_PATH")+
	 "/Status__"+rpa->gen.Timer().TimeString(3));
      msg_Error()<<METHOD<<"(): Pre-crash status saved to:\n'"
		 <<rpa->gen.Variable("SHERPA_STATUS_PATH")<<"'\n"<<std::endl;
      MakeDir(rpa->gen.Variable("SHERPA_STATUS_PATH"));
    }


  switch (signal) {

  case SIGSEGV:
    ++m_nsegv;
    GenerateStackTrace(std::cout,false);
    if (m_nsegv>3) {
      msg_Error()<<om::reset<<"   Abort immediately."<<om::reset<<std::endl;
      Abort(1);
    }
    break;
    
  case SIGABRT:
    Abort();
    break;

  case SIGTERM:

  case SIGXCPU:
    msg_Error()<<om::reset<<"   Cannot continue."<<om::reset<<std::endl;
    Terminate(2);
    break;

  case SIGINT:
    Terminate(1);
    break;

  case SIGBUS:
    ++m_nbus;
    if (m_nbus>3) {
      msg_Error()<<om::reset<<"   Abort immediately."<<om::reset<<std::endl;
      Abort(1);
    }
    GenerateStackTrace(std::cout,false);
    msg_Error()<<om::reset<<"   Cannot continue."<<om::reset<<std::endl;
    Terminate(3);
    break;

  case SIGFPE:
    msg_Error()<<"   Floating point exception."<<om::reset<<std::endl;
    Terminate(1);
    break;

  case SIGPIPE:
    msg_Error()<<"   Pipe closed. Will stop writing."<<om::reset<<std::endl;
    Terminate(0);
    break;

  case SIGUSR1:
    Terminate(1);
    break;

  default:
    msg_Error()<<"   Cannot handle signal."<<om::reset<<std::endl;
    Terminate(1);
  }
  
}
