#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"
#include <iostream>
#include <typeinfo>

using namespace ATOOLS;

Exception::Exception(const std::string& type,
		     const std::string& info,
		     const std::string& cmethod,
         const std::string& file,
         const int& line):
  m_info(info), m_type(type), m_file(file) ,m_line(line)
{
  std::string cmethod_
    =cmethod.substr(0,ATOOLS::Min(cmethod.length(),cmethod.find("(")));
  size_t pos;
  while ((pos=cmethod_.find(" "))!=std::string::npos) 
    cmethod_=cmethod_.substr(pos+1);
  pos=cmethod_.find("::");
  while (pos!=std::string::npos) {
    m_class=cmethod_.substr(0,pos);
    cmethod_=cmethod_.substr(pos+2);
    pos=cmethod_.find("::");
    m_method=cmethod_.substr(0,ATOOLS::Min(cmethod_.length(),pos));
  }
}

std::ostream &ATOOLS::operator<<(std::ostream &str,
				 const Exception &exception)
{
    str<<om::bold<<om::red
        <<exception.TypeName()<<om::reset
        <<om::bold<<" thrown" <<om::reset;

    if (exception.m_class.length()>0)
        str <<om::bold<< " in "<<om::reset
        <<om::blue
        <<exception.m_class<<"::"
        <<exception.m_method
        <<om::reset;

    if (exception.m_file.length() || exception.m_line > 0){
        str << om::green << "(in ";
        if (exception.m_file.length()) 
            str << exception.m_file;
        if (exception.m_line > 0) 
            str << " line " << exception.m_line;
        str << ")" << om::reset;
    }

    str<<":\n"<<om::red
        <<exception.m_info
        <<om::reset;

    return str;
}
