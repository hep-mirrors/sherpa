#include "ATOOLS/Org/Message.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include <sys/stat.h>
#include <iterator>

namespace ATOOLS {
  Message *msg(NULL);
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


indentbuf::indentbuf(std::streambuf* basebuf) :
  m_basebuf(basebuf), m_indent(0), at_start(true)
{
}

indentbuf::~indentbuf()
{
}

void indentbuf::Indent(size_t i)
{
  m_indent+=i;
}

void indentbuf::DeIndent(size_t i)
{
  if (m_indent>=i) m_indent-=i;
}

std::streambuf::int_type indentbuf::overflow(int_type ch)
{
  if (traits_type::eq_int_type(ch, traits_type::to_int_type('\r'))) {
  }
  if (ch == traits_type::eof())
    return traits_type::not_eof(ch);

  if (traits_type::not_eof(ch)) {
    if (at_start)
      for (size_t i = 0; i < m_indent; ++i)
        m_basebuf->sputc(traits_type::to_char_type(' '));
    m_basebuf->sputc(traits_type::to_char_type(ch));
    if (traits_type::eq_int_type(ch, traits_type::to_int_type('\n')))
      at_start = true;
    else
      at_start = false;
  }
  return ch; 
}

Message::Message() :
  m_devnull("/dev/null", std::ios::app),
  m_buf(std::cout.rdbuf()),
  p_log(NULL),
  m_output(std::cout.rdbuf()),
  m_error(std::cerr.rdbuf())
{
  m_logfile = "";
  m_level = 0;
  m_modifiable = true;
  m_mpimode = 0;
}

Message::~Message() 
{
  SetOutStream(m_buf);
  if (p_log) delete p_log;
}

void Message::Init()
{ 
  Settings& s = Settings::GetMainSettings();

  std::string logfile = s["LOG_FILE"].Get<std::string>();
  if (logfile!="") {
    m_logfile = logfile;
    p_log = new std::ofstream(logfile.c_str(),std::ios::app);
    SetOutStream(*p_log);
  }

  m_buf.SetBaseBuf(m_output.rdbuf());
  SetOutStream(m_buf);

  // set general output level
  m_level = s["OUTPUT"].SetDefault(2).Get<int>();

  // set function-specific output level
  auto fctoutput = s["FUNCTION_OUTPUT"];
  for (const auto& fctname : fctoutput.GetKeys()) {
    const auto level = fctoutput[fctname].SetDefault(m_level).Get<int>();
    if (level & 1)  m_contextevents.insert(fctname);
    if (level & 2)  m_contextinfo.insert(fctname);
    if (level & 4)  m_contexttracking.insert(fctname);
    if (level & 8)  m_contextdebugging.insert(fctname);
    if (level & 32) m_contextiodebugging.insert(fctname);
  }

  m_mpimode = s["MPI_OUTPUT"].SetDefault(0).Get<int>();
}

void Message::SetStandard() 
{
  SetOutStream(std::cout);
  SetErrStream(std::cerr);
}

std::ostream &Message::Out()
{ 
#ifdef USING__MPI
  if (!m_mpimode && 
      mpi->Rank()) return m_devnull;
#endif
  return m_output; 
}

std::ostream &Message::Error()
{ 
#ifdef USING__MPI
  if (!m_mpimode && 
      mpi->Rank()) return m_devnull;
#endif
  if (m_level >= 0) return m_output; 
  return m_devnull; 
}

std::ostream &Message::Events()
{ 
#ifdef USING__MPI
  if (!m_mpimode && 
      mpi->Rank()) return m_devnull;
#endif
  if (m_level & 1) return m_output; 
  return m_devnull;  
}

std::ostream &Message::Info()
{ 
#ifdef USING__MPI
  if (!m_mpimode && 
      mpi->Rank()) return m_devnull;
#endif
  if (m_level & 2) return m_output; 
  return m_devnull;  
}

std::ostream &Message::Tracking()
{ 
#ifdef USING__MPI
  if (!m_mpimode && 
      mpi->Rank()) return m_devnull;
#endif
  if (m_level & 4) return m_output; 
  return m_devnull;  
}

std::ostream &Message::Debugging()
{ 
#ifdef USING__MPI
  if (!m_mpimode && 
      mpi->Rank()) return m_devnull;
#endif
  if (m_level & 8) return m_output; 
  return m_devnull;  
}

std::ostream &Message::IODebugging()
{
#ifdef USING__MPI
  if (!m_mpimode && 
      mpi->Rank()) return m_devnull;
#endif
  if (m_level & 32) return m_output;
  return m_devnull;
}

std::string Message::ExtractMethodName(std::string cmethod) const   
{ 
  for (size_t pos(cmethod.find(", "));
       pos!=std::string::npos;pos=cmethod.find(", ")) cmethod.erase(pos+1,1);
  for (size_t pos(cmethod.find("> >"));
       pos!=std::string::npos;pos=cmethod.find("> >")) cmethod.erase(pos+1,1);
  std::string cclass("<no class>"), method("<no method>");
  cmethod=cmethod.substr(0,ATOOLS::Min(cmethod.length(),cmethod.find("(")));
  size_t pos;
  while ((pos=cmethod.find(" "))!=std::string::npos) 
    cmethod=cmethod.substr(pos+1);
  pos=cmethod.find("::");
  for (size_t bpos(cmethod.find("<"));pos!=std::string::npos && pos<bpos;bpos-=pos+2) {
    cclass=cmethod.substr(0,pos);
    cmethod=cmethod.substr(pos+2);
    pos=cmethod.rfind("::");
    method=cmethod.substr(0,ATOOLS::Min(cmethod.length(),pos));
  }
  if (cclass=="<no class>") return cmethod;
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

bool Message::LevelIsIODebugging(const std::string& context) const
{
  for (std::set<std::string>::reverse_iterator rit=m_contextiodebugging.rbegin();
       rit!=m_contextiodebugging.rend(); ++rit) {
    if (context.find(*rit)!=std::string::npos) return true;
  }
  return false;
}
