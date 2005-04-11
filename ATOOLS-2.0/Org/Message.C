#include "Message.H"

#include "Run_Parameter.H"
#include "MyStrStream.H"
#include <sys/stat.h>

namespace ATOOLS {
  Message msg;
}

template <class Value_Type>
std::string ToString(const Value_Type value) {
  MyStrStream converter;
  std::string converted;
  converter<<value;
  converter>>converted;
  return converted;
}

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &str,const bm::code modifier) 
{
  switch (modifier) {
  case bm::back: return msg.Modifiable()?str<<"\b":str<<" \\b ";
  case bm::cr:   return msg.Modifiable()?str<<"\r":str<<"\n";
  case bm::bell: return msg.Modifiable()?str<<"\a":str<<" \\a ";
  case bm::none: return str;
  }
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const om::code modifier) 
{
  if (!msg.Modifiable()) return str;
  switch (modifier) {
#ifdef USING__COLOUR
  case om::reset:    return str<<"\e[0m";
  case om::bold:     return str<<"\e[1m";
  case om::underln:  return str<<"\e[4m";
  case om::blink:    return str<<"\e[5m";
  case om::blackbg:  return str<<"\e[7m";
  case om::red:      return str<<"\e[31m";
  case om::green:    return str<<"\e[32m";
  case om::brown:    return str<<"\e[33m";
  case om::blue:     return str<<"\e[34m";
  case om::violet:   return str<<"\e[35m";
  case om::lblue:    return str<<"\e[36m";
  case om::grey:     return str<<"\e[37m";
  case om::redbg:    return str<<"\e[41m";
  case om::greenbg:  return str<<"\e[42m";
  case om::brownbg:  return str<<"\e[43m";
  case om::bluebg:   return str<<"\e[44m";
  case om::violetbg: return str<<"\e[45m";
  case om::lbluebg:  return str<<"\e[46m";
  case om::greybg:   return str<<"\e[47m";
  case om::none:     return str;
#else
  default: return str;
#endif
  }
  return str;
}
 
std::ostream &ATOOLS::operator<<(std::ostream &str,const mm modifier)
{
  if (!msg.Modifiable()) return str;
  switch (modifier.m_code) {
#ifdef USING__COLOUR
  case mm::up:    return str<<"\e["<<modifier.m_num<<"A";
  case mm::down:  return str<<"\e["<<modifier.m_num<<"B";
  case mm::right: return str<<"\e["<<modifier.m_num<<"C";
  case mm::left:  return str<<"\e["<<modifier.m_num<<"D";
  case mm::none:  return str;
#else
  default: return str;
#endif
  }
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const tm::code modifier) 
{
  if (!msg.Modifiable()) return str;
  switch (modifier) {
#ifdef USING__COLOUR
  case tm::curon:  return str<<"\e[?25h";
  case tm::curoff: return str<<"\e[?25l";
  case mm::none:   return str;
#else
  default: return str;
#endif
  }
  return str;
}

Message::Message() 
{      
  p_output = &std::cout;
  p_error = &std::cerr;
  p_no = new std::ofstream("/dev/null",std::ios::app);
  p_logfile = p_no;
  m_file = 0;
  m_level = 0;
  m_modifiable = false;
  m_indent = 0;
}

Message::~Message() 
{      
  if (p_logfile!=p_no) delete p_logfile;
}

void Message::Init(const int level,const std::string &logfile) 
{ 
  m_level = level; 
  if (m_level&16) {
    InitLogFile(logfile);
    Out()<<"Initialize output module Message. Level "<<m_level<<std::endl;
  }
}

void Message::InitLogFile(const std::string &logfile) 
{
  if (p_logfile!=p_no) return;
  std::string name=logfile;
  if (name==std::string("")) {
    int i=0;
    do { 
      name=std::string("./sherpa_log_")+ToString(++i)+std::string(".log"); 
      std::ifstream testfile(name.c_str());
      if (!testfile.is_open()) break;
    } while (true);
  }
  std::string command="echo \"! ****************************************\n";
  command+="! *           Sherpa Log File            *\n";
  command+="! ****************************************\n\" > ";
  system((command+name).c_str());
  command="echo \"! starting Sherpa 1.0.6 at '$HOSTNAME' ";
  command+="(architecture $HOSTTYPE) \n! on `date`\n\" >> ";
  system((command+name).c_str());
  command="echo \"! +--------------------------------------+\n";
  command+="! |            shell settings            |\n";
  command+="! +--------------------------------------+\n\" >> ";
  system((command+name).c_str());
  command="echo \"! PATH=$PATH\n\n! LD_LIBRARY_PATH=$LD_LIBRARY_PATH\n\n! CLHEPDIR=$CLHEPDIR\n";
  command+="! PWD=$PWD\n\n! Sherpa was compiled for gcc `gcc -dumpversion`\n\" >> ";
  system((command+name).c_str());
  p_logfile = new std::ofstream(name.c_str(),std::ios::app);
  (*p_logfile)<<"! +--------------------------------------+"<<std::endl;
  (*p_logfile)<<"! |            run parameters            |"<<std::endl;
  (*p_logfile)<<"! +--------------------------------------+\n"<<std::endl;
}

void Message::SetStandard() 
{
  if (m_file==1) delete p_output;
  p_output = &std::cout;
  p_error = &std::cerr;
  m_file = 0;
}

void Message::SetPrecision(const int precision) 
{
  if (p_output) p_output->precision(precision);
}

std::ostream &Message::Out() const 
{ 
  return *p_output<<Indent(); 
}

std::ostream &Message::Error() const     
{ 
  if (m_level >= 0) return *p_output<<Indent(); 
  return *p_no; 
}

std::ostream &Message::Events() const    
{ 
  if (m_level & 1) return *p_output<<Indent(); 
  return *p_no;  
}

std::ostream &Message::Info() const      
{ 
  if (m_level & 2) return *p_output<<Indent(); 
  return *p_no;  
}

std::ostream &Message::Tracking() const  
{ 
  if (m_level & 4) return *p_output<<Indent(); 
  return *p_no;  
}

std::ostream &Message::Debugging() const 
{ 
  if (m_level & 8) return *p_output<<Indent(); 
  return *p_no;  
}

std::ostream &Message::LogFile() const   
{ 
  if (m_level & 16) return *p_logfile<<Indent(); 
  return *p_no; 
}
