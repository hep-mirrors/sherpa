#include "Message.H"
#include "Run_Parameter.H"

namespace ATOOLS {
  Message msg;
}

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &str,const om::code modifier) 
{
  if (!msg.Modifiable()) return str;
  switch (modifier) {
#ifdef USING__COLOUR
  case om::reset:   return str<<"\e[0m";
  case om::bold:    return str<<"\e[1m";
  case om::blink:   return str<<"\e[5m";
  case om::backgnd: return str<<"\e[6m";
  case om::red:     return str<<"\e[31m";
  case om::green:   return str<<"\e[32m";
  case om::brown:   return str<<"\e[33m";
  case om::blue:    return str<<"\e[34m";
  case om::violet:  return str<<"\e[35m";
  case om::lblue:   return str<<"\e[36m";
  case om::grey:    return str<<"\e[37m";
  case om::none:    return str;
#else
  default: return str;
#endif
  }
}
 
Message::Message() 
{      
  m_output = &std::cout;
  m_error = &std::cerr;
  m_no = new std::ofstream("/dev/null",std::ios::app);
  m_file = 0;
  m_level = 0;
  m_modifiable = false;
}

void Message::Init(const int level) 
{ 
  m_level = level; 
  if (m_level>2) (*m_output)<<"Initialize output module Message. Level "<<m_level<<std::endl;
}

void Message::SetFile(const char *file="Amegic.log") 
{
  if (m_file==0) {
    m_output = new std::ofstream(file);
    m_error = m_output;
  }
  m_file = 1;
}

void Message::SetStandard() 
{
  if (m_file==1) delete m_output;
  m_output = &std::cout;
  m_error = &std::cerr;
  m_file = 0;
}

void Message::SetPrecision(const int precision) 
{
  if (m_output) m_output->precision(precision);
}

