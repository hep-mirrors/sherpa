#include "Exception.H"

#include "Run_Parameter.H"
#include <iostream>

using namespace ATOOLS;

unsigned int Exception::s_exitcode=1;
std::vector<Exception::Tester> Exception::s_testers=std::vector<Exception::Tester>();
std::vector<Exception::Terminator> Exception::s_terminators=std::vector<Exception::Terminator>();
std::vector<Terminator_Base*> Exception::s_objects=std::vector<Terminator_Base*>();

void Terminator_Base::Terminate()
{
  msg.Error()<<"Terminator_Base::Terminate(): "
	     <<"Vitual function called !"<<std::endl;
}

Exception::Exception(const ex::type type,const std::string info):
  m_type(type),
  m_info(info)
{
  SetExitCode();
}

Exception::Exception(const ex::type type,const std::string info,
		     const std::string cclass,const std::string cmethod):
  m_type(type),
  m_info(info),
  m_class(cclass),
  m_method(cmethod) 
{
  SetExitCode();
}

void Exception::Terminate() 
{
  msg.Error()<<"ATOOLS::Terminate(): Preparing termination ..."<<std::endl;
  PrepareTerminate();
  Exit(s_exitcode);
}

void Exception::Exit(int exitcode)
{
  msg.Error()<<"Irregularly terminating Sherpa: exit code ("
	     <<exitcode<<")"<<std::endl;
  exit(exitcode);
}

void Exception::SEGVHandler(int signal) 
{
  static unsigned int segv_signals;
  ++segv_signals;
  if (segv_signals>3) {
    msg.Error()<<"ATOOLS::SEGVHandler("<<signal<<"): Caught SIGSEGV 3 times. "<<std::endl
	       <<"   Abort immediately."<<std::endl;
    abort();
  }
  msg.Error()<<"ATOOLS::SEGVHandler("<<signal<<"): Signal SIGSEGV caught. "<<std::endl
	     <<"   Cannot run further. Preparing termination ..."<<std::endl;
  PrepareTerminate();
  Exit(2);
}

void Exception::INTHandler(int signal) 
{
  std::string input="y";
  if (!rpa.gen.BatchMode()) {
    msg.Error()<<"ATOOLS::INTHandler("<<signal<<"): Signal SIGINT caught."<<std::endl
	       <<"   Do you want to stop the program ? ";
    std::cin>>input;
  }
  if (input!="y" && input!="Y") return;
  msg.Error()<<"ATOOLS::INTHandler("<<signal<<"): "
	     <<"Preparing termination ..."<<std::endl;
  PrepareTerminate();
  Exit(1);
}

void Exception::FPEHandler(int signal) 
{
  msg.Error()<<"ATOOLS::FPEHandler("<<signal<<"): "<<std::endl
	     <<"   Sherpa does not throw floating point exceptions."<<std::endl;
}

bool Exception::ApproveTerminate() const
{
  if (s_testers.size()==0) return true;
  for (size_t i=0;i<s_testers.size();++i) if (s_testers[i](*this)) return true;
  return false;
}

void Exception::PrepareTerminate()
{
  for (size_t i=0;i<s_objects.size();++i) s_objects[i]->Terminate(); 
  for (size_t i=0;i<s_terminators.size();++i) s_terminators[i]();
}

void Exception::UpdateLogFile() const 
{
  msg.LogFile()<<"Sherpa";
  if (m_class.length()>0) {
    msg.LogFile()<<" : "<<m_class<<"::"<<m_method;
  }
  msg.LogFile()<<" throws "<<m_type<<": "<<m_info;
}

void Exception::RemoveObject(Terminator_Base *const object)
{
  for (std::vector<Terminator_Base*>::iterator oit=s_objects.begin();
       oit!=s_objects.end();++oit) if (*oit==object) s_objects.erase(oit); 
}

void Exception::SetExitCode()
{
  s_exitcode=1;
  if (m_class=="ISR_Handler") s_exitcode=151;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const ex::type &type)
{
  switch (type) {
  case ex::normal_exit   : return str<<"normal exit";
  case ex::critical_error: return str<<"critical error";
  case ex::fatal_error   : return str<<"fatal error";
  case ex::unknown       : return str<<"unknown exception";
  }
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const Exception &exception)
{
  str<<om::bold<<"Sherpa";
  if (exception.m_class.length()>0) {
    str<<": "<<om::reset<<om::blue<<exception.m_class
       <<"::"<<exception.m_method;
  }
  return str<<om::reset<<" throws "<<om::bold<<om::red
	    <<exception.m_type<<om::reset<<om::bold<<": "<<om::reset<<om::red
	    <<exception.m_info<<om::reset;
}

