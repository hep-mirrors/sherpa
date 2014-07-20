#include "MODEL/Interaction_Models/Color_Function.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace MODEL;
using namespace ATOOLS;

bool Color_Function::operator==(const Color_Function &c) const
{
  if (m_type!=c.m_type) return false;
  if (m_partarg[0]!=c.m_partarg[0]) return false;
  if (m_partarg[1]!=c.m_partarg[1]) return false;
  if (m_partarg[2]!=c.m_partarg[2]) return false;
  if (m_strarg[0]!=c.m_strarg[0]) return false;
  if (m_strarg[1]!=c.m_strarg[1]) return false;
  if (m_strarg[2]!=c.m_strarg[2]) return false;
  if (((bool)p_next)^((bool)c.p_next)) return false;
  if (p_next) return *p_next==*c.p_next;
  return true;
}

Color_Function & Color_Function::operator=(const Color_Function & c)
{
  if (this!=&c) {
    m_type       = c.m_type;
    m_partarg[0] = c.m_partarg[0];
    m_partarg[1] = c.m_partarg[1];
    m_partarg[2] = c.m_partarg[2];
    m_strarg[0]  = c.m_strarg[0];
    m_strarg[1]  = c.m_strarg[1];
    m_strarg[2]  = c.m_strarg[2]; 
    if (p_next) delete p_next; 
    if (c.p_next!=0)
      p_next=new Color_Function(*(c.p_next));
    else
      p_next = 0;
  }
  return *this;
}

std::string Color_Function::PID() const
{
  std::string pid(ATOOLS::ToString(m_type));
  if (p_next) pid+="*"+p_next->PID();
  return pid;
}

std::string Color_Function::FullString() const
{
  std::string pid(String());
  if (p_next) pid+="*"+p_next->FullString();
  return pid;  
}

std::string Color_Function::String() const {
      
  if (m_type==cf::None) return std::string("1");
      
  std::string help;
  switch (m_type) {
  case 0 :  help = std::string("T[ , , ]");break;
  case 1 :  help = std::string("F[ , , ]");break;
  case 2 :  help = std::string("D[ , ]");break;
  case 4 :  help = std::string("G[ , ]");break;
  default : return std::string("1");
  }
  for (short int i=0;i<3;i++) {
    if ((m_type==cf::D || m_type==cf::G) && i==2) break;
    help[2+i*2] = m_strarg[i];
  }
  return help;
}

Color_Function::~Color_Function() {
  if (p_next) delete p_next;
}

std::ostream &MODEL::operator<<(std::ostream &str,const cf::code &c)
{
  switch (c) {
  case cf::T: return str<<"T"; 
  case cf::F: return str<<"F"; 
  case cf::D: return str<<"D"; 
  case cf::None: return str<<"None"; 
  case cf::G: return str<<"G"; 
  case cf::Unknown: return str<<"Unknown";
  }
  THROW(fatal_error,"Invalid code '"+ToString((int)c)+"'");
  return str;
}
