#include "Message.H"

#include "Run_Parameter.H"
#include "MyStrStream.H"
#include <sys/stat.h>

namespace ATOOLS {
  Message msg;
}

template <class Value_Type>
std::string ToString(const Value_Type value) {
  std::stringstream converter;
  std::string converted;
  converter<<value;
  converter>>converted;
  return converted;
}

using namespace ATOOLS;

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
 
Message::Message() 
{      
  p_logfile = NULL;
  p_output = &std::cout;
  p_error = &std::cerr;
  p_no = new std::ofstream("/dev/null",std::ios::app);
  m_file = 0;
  m_level = 0;
  m_modifiable = false;
}

Message::~Message() 
{      
  if (p_logfile!=NULL) delete p_logfile;
}

void Message::Init(const int level,const std::string &logfile) 
{ 
  m_level = level; 
  if (m_level&&4) {
    InitLogFile(logfile);
    (*p_output)<<"Initialize output module Message. Level "<<m_level<<std::endl;
  }
}

void Message::InitLogFile(const std::string &logfile) 
{
  if (p_logfile!=NULL) return;
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
  command="echo \"! starting Sherpa 1.0.3 at '$HOSTNAME' ";
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

