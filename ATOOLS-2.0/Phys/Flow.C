#include "Flow.H"

using namespace ATOOLS;

unsigned int Flow::s_qcd_counter=600;

std::ostream& ATOOLS::operator<<(std::ostream &ostr,const Flow &flow)
{
  ostr<<"[";
  for (std::map<unsigned int,unsigned int>::const_iterator fit=flow.m_code.begin();
       fit!=flow.m_code.end();++fit) ostr<<"("<<fit->first<<"="<<fit->second<<")";
  return ostr<<"]";
}

Flow::Flow(Particle *owner): 
  p_owner(owner) 
{ 
  for (short unsigned int i=1;i<3;++i) m_code[i]=0;
}

Flow::~Flow() 
{ 
  m_code.clear(); 
}

unsigned int Flow::Code(const unsigned int index) const
{
  std::map<unsigned int,unsigned int>::const_iterator count=m_code.find(index);
  if (count!=m_code.end()) return count->second;
  return 0;
}

void Flow::SetCode(const unsigned int index,const int code) 
{
  if (code==-1) m_code[index]=++s_qcd_counter; 
  else m_code[index]=(unsigned int)code;
}
















