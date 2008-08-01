#include "Message.H"

#include "Run_Parameter.H"
#include "MyStrStream.H"

#include <sys/stat.h>
#include <iterator>

namespace ATOOLS {
  Message *msg(new Message());
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
  case bm::back: return msg->Modifiable()?str<<"\b":str<<" \\b ";
  case bm::cr:   return msg->Modifiable()?str<<"\r":str<<"\n";
  case bm::bell: return msg->Modifiable()?str<<"\a":str<<" \\a ";
  case bm::none: return str;
  }
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const om::code modifier) 
{
  if (!msg->Modifiable()) return str;
  switch (modifier) {
#ifdef USING__COLOUR
  case om::reset:    return str<<"\033[0m";
  case om::bold:     return str<<"\033[1m";
  case om::underln:  return str<<"\033[4m";
  case om::blink:    return str<<"\033[5m";
  case om::blackbg:  return str<<"\033[7m";
  case om::red:      return str<<"\033[31m";
  case om::green:    return str<<"\033[32m";
  case om::brown:    return str<<"\033[33m";
  case om::blue:     return str<<"\033[34m";
  case om::violet:   return str<<"\033[35m";
  case om::lblue:    return str<<"\033[36m";
  case om::grey:     return str<<"\033[37m";
  case om::redbg:    return str<<"\033[41m";
  case om::greenbg:  return str<<"\033[42m";
  case om::brownbg:  return str<<"\033[43m";
  case om::bluebg:   return str<<"\033[44m";
  case om::violetbg: return str<<"\033[45m";
  case om::lbluebg:  return str<<"\033[46m";
  case om::greybg:   return str<<"\033[47m";
  case om::none:     return str;
#else
  default: return str;
#endif
  }
  return str;
}
 
std::ostream &ATOOLS::operator<<(std::ostream &str,const mm modifier)
{
  if (!msg->Modifiable()) return str;
  switch (modifier.m_code) {
#ifdef USING__COLOUR
  case mm::up:    return str<<"\033["<<modifier.m_num<<"A";
  case mm::down:  return str<<"\033["<<modifier.m_num<<"B";
  case mm::right: return str<<"\033["<<modifier.m_num<<"C";
  case mm::left:  return str<<"\033["<<modifier.m_num<<"D";
  case mm::none:  return str;
#else
  default: return str;
#endif
  }
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const tm::code modifier) 
{
  if (!msg->Modifiable()) return str;
  switch (modifier) {
#ifdef USING__COLOUR
  case tm::curon:  return str<<"\033[?25h";
  case tm::curoff: return str<<"\033[?25l";
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
  delete p_no;
}

void Message::Init(const std::string& level,const std::string &logfile) 
{ 
  Data_Reader dr("|","[","!");
  dr.AddLineSeparator("]");
  dr.AddComment("#");
  dr.SetString(level);
  std::vector<std::vector<std::string> > svv;
  dr.MatrixFromString(svv);

  m_level=2;
  for (size_t i=0;i<svv.size();++i) {
    if (svv[i].size()==1) m_level=ATOOLS::ToType<int>(svv[i][0]);
    else if (svv[i].size()==2) {
      int lev=ToType<int>(svv[i][1]);
      if (lev&1) m_contextevents.insert(svv[i][0]);
      if (lev&2) m_contextinfo.insert(svv[i][0]);
      if (lev&4) m_contexttracking.insert(svv[i][0]);
      if (lev&8) m_contextdebugging.insert(svv[i][0]);
    }
    else {
      Out()<<METHOD<<": Do not understand output specification: "<<std::endl;
      copy(svv[i].begin(), svv[i].end(),
           std::ostream_iterator<std::string>(Out(), " "));
    }
  }
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
  std::string gccv;
  FILE *pout(popen("gcc -dumpversion","r"));
  int cch;
  while ((cch=fgetc(pout))!=EOF) gccv+=(char)cch;
  pclose(pout);
  if (gccv.find('\n')!=std::string::npos) 
    gccv=gccv.substr(0,gccv.find('\n'));
  char *var=NULL;
  p_logfile = new std::ofstream(name.c_str(),std::ios::app);
  (*p_logfile)<<"! ****************************************\n"
	      <<"! *           Sherpa Log File            *\n"
	      <<"! ****************************************\n";
  (*p_logfile)<<"! starting Sherpa 1.0.6 at '"
	      <<std::string(((var=getenv("HOSTNAME"))==NULL?"":var))
	      <<"' (architecture "
	      <<std::string(((var=getenv("HOSTTYPE"))==NULL?"":var))
	      <<") \n! on "<<rpa.gen.Timer().TimeString()<<"\n";
  (*p_logfile)<<"! +--------------------------------------+\n"
	      <<"! |            shell settings            |\n"
	      <<"! +--------------------------------------+\n"
	      <<"! PATH="<<rpa.gen.Variable("PATH")<<"\n\n"
	      <<"! LD_LIBRARY_PATH="<<rpa.gen.Variable("LD_LIBRARY_PATH")
	      <<"\n\n! PWD="<<rpa.gen.Variable("SHERPA_RUN_PATH")
	      <<"\n\n! Sherpa was compiled for gcc "<<gccv<<"\n\n";
  (*p_logfile)<<"! +--------------------------------------+\n"
	      <<"! |            run parameters            |\n"
	      <<"! +--------------------------------------+\n!\n";
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
  if (m_level >= 0) return *p_output; 
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

std::string Message::ExtractMethodName(std::string cmethod) const   
{ 
  std::string cclass("<no class>"), method("<no method>");
  cmethod=cmethod.substr(0,ATOOLS::Min(cmethod.length(),cmethod.find("(")));
  size_t pos;
  while ((pos=cmethod.find(" "))!=std::string::npos) 
    cmethod=cmethod.substr(pos+1);
  pos=cmethod.find("::");
  while (pos!=std::string::npos) {
    cclass=cmethod.substr(0,pos);
    cmethod=cmethod.substr(pos+2);
    pos=cmethod.find("::");
    method=cmethod.substr(0,ATOOLS::Min(cmethod.length(),pos));
  }
  return cclass+"::"+cmethod;
}

bool Message::LevelIsEvents(const std::string& context) const
{
  for (std::set<std::string>::reverse_iterator rit=m_contextevents.rbegin();
       rit!=m_contextevents.rend(); ++rit) {
    if (context.find(*rit)!=std::string::npos) return true;
  }
  return false;
}

bool Message::LevelIsInfo(const std::string& context) const
{
  for (std::set<std::string>::reverse_iterator rit=m_contextinfo.rbegin();
       rit!=m_contextinfo.rend(); ++rit) {
    if (context.find(*rit)!=std::string::npos) return true;
  }
  return false;
}

bool Message::LevelIsTracking(const std::string& context) const
{
  for (std::set<std::string>::reverse_iterator rit=m_contexttracking.rbegin();
       rit!=m_contexttracking.rend(); ++rit) {
    if (context.find(*rit)!=std::string::npos) return true;
  }
  return false;
}

bool Message::LevelIsDebugging(const std::string& context) const
{
  for (std::set<std::string>::reverse_iterator rit=m_contextdebugging.rbegin();
       rit!=m_contextdebugging.rend(); ++rit) {
    if (context.find(*rit)!=std::string::npos) return true;
  }
  return false;
}
