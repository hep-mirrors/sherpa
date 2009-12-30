#include "ATOOLS/Org/Exception.H"

#include "ATOOLS/Math/MathTools.H"
#include <iostream>
#include <typeinfo>

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
  case ex::missing_input       : return str<<"missing input";
  case ex::missing_module      : return str<<"missing module";
  case ex::unknown             : return str<<"unknown exception";
  }
  return str;
}

Tester_Object::~Tester_Object()
{
}

bool Tester_Object::ApproveTerminate()
{
  msg_Error()<<METHOD<<"() ["<<typeid(this).name()
	     <<"]: Virtual function called !"<<std::endl;
  exh->GenerateStackTrace(std::cout);
  return true;
}

Terminator_Object::~Terminator_Object()
{
}

bool Terminator_Object::ReadInStatus(const std::string &path)
{
  return true;
}

void Terminator_Object::PrepareTerminate()
{
  msg_Error()<<METHOD<<"() ["<<typeid(this).name()
	     <<"]: Virtual function called !"<<std::endl;
}

Exception::Exception(const ex::type type,const std::string info):
  m_type(type),
  m_info(info)
{
  exh->m_exception=this;
}

Exception::Exception(const ex::type type,const std::string info,
		     std::string cmethod):
  m_type(type),
  m_info(info)
{
  cmethod=cmethod.substr(0,ATOOLS::Min(cmethod.length(),cmethod.find("(")));
  size_t pos;
  while ((pos=cmethod.find(" "))!=std::string::npos) 
    cmethod=cmethod.substr(pos+1);
  pos=cmethod.find("::");
  while (pos!=std::string::npos) {
    m_class=cmethod.substr(0,pos);
    cmethod=cmethod.substr(pos+2);
    pos=cmethod.find("::");
    m_method=cmethod.substr(0,ATOOLS::Min(cmethod.length(),pos));
  }
  exh->m_exception=this;
}

Exception::Exception(const ex::type type,const std::string info,
		     const std::string cclass,const std::string cmethod):
  m_type(type),
  m_info(info),
  m_class(cclass),
  m_method(cmethod) 
{
  exh->m_exception=this;
}

Exception::~Exception() 
{
  if (exh->m_exception==this) {
    exh->SetExitCode();
    exh->m_exception=NULL;
  }
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

