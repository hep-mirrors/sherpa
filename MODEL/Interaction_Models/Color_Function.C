#include "MODEL/Interaction_Models/Color_Function.H"


using namespace MODEL;

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

std::string Color_Function::String() {
      
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
