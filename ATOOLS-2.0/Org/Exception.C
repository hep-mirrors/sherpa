#include "Exception.H"

#include <iostream>

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &str,const ex::type &type)
{
  switch (type) {
  case ex::normal_exit         : return str<<"normal exit";
  case ex::unknown_option      : return str<<"unknown option";
  case ex::inconsistent_option : return str<<"inconsistent option";
  case ex::not_implemented     : return str<<"not implemented";
  case ex::critical_error      : return str<<"critical error";
  case ex::fatal_error         : return str<<"fatal error";
  case ex::unknown             : return str<<"unknown exception";
  }
  return str;
}

Tester_Object::~Tester_Object()
{
}

bool Tester_Object::ApproveTerminate()
{
  msg.Error()<<"Tester_Oject::ApproveTerminate(): "
	     <<"Virtual function called !"<<std::endl;
  return true;
}

Terminator_Object::~Terminator_Object()
{
}

void Terminator_Object::PrepareTerminate()
{
  msg.Error()<<"Terminator_Object::PrepareTerminate(): "
	     <<"Virtual function called !"<<std::endl;
}

Exception::Exception(const ex::type type,const std::string info):
  m_type(type),
  m_info(info)
{
  Exception_Handler::s_exception=this;
}

Exception::Exception(const ex::type type,const std::string info,
		     const std::string cclass,const std::string cmethod):
  m_type(type),
  m_info(info),
  m_class(cclass),
  m_method(cmethod) 
{
  Exception_Handler::s_exception=this;
}

Exception::~Exception() 
{
  if (Exception_Handler::s_exception==this) {
    Exception_Handler::SetExitCode();
    Exception_Handler::s_exception=NULL;
  }
}

void Exception::UpdateLogFile() const 
{
  msg.LogFile()<<"Sherpa";
  if (m_class.length()>0) {
    msg.LogFile()<<" : "<<m_class<<"::"<<m_method;
  }
  msg.LogFile()<<" throws "<<m_type<<": "<<m_info;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const Exception &exception)
{
  str<<om::bold<<"Sherpa";
  if (exception.m_class.length()>0) {
    str<<": "<<om::reset<<om::blue<<exception.m_class
       <<"::"<<exception.m_method;
  }
  return str<<om::reset<<" throws "<<om::bold<<om::red
	    <<exception.m_type<<om::reset<<om::bold<<": "<<std::endl<<"   "<<om::reset<<om::red
	    <<exception.m_info<<om::reset;
}

