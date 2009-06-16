#include "ATOOLS/Phys/Flow.H"

using namespace ATOOLS;

unsigned int Flow::s_qcd_counter=600;

namespace ATOOLS {
std::ostream& operator<<(std::ostream &ostr,const Flow &flow)
{
  ostr<<"[";
  for (std::map<unsigned int,unsigned int>::const_iterator 
	 fit=flow.m_code.begin();fit!=flow.m_code.end();++fit) 
    ostr<<"("<<fit->first<<"="<<fit->second<<")";
  return ostr<<"]";
}
}

Flow::Flow(Particle *owner): 
  p_owner(owner) 
{ 
  for (short unsigned int i=1;i<3;++i) m_code[i]=0;
}

Flow::Flow(const Flow &flow): 
  p_owner(flow.p_owner) 
{ 
  for (unsigned int i=1;i<3;++i) m_code[i]=flow.m_code.find(i)->second;
}

Flow::~Flow() {}

void Flow::SetCode(const unsigned int index,const int code) 
{
  if (code==-1) m_code[index]=++s_qcd_counter; 
  else m_code[index]=(unsigned int)code;
}

void Flow::SetCode(const Flow &flow)
{
  m_code=flow.m_code;
}

unsigned int Flow::Code(const unsigned int index) const
{
  std::map<unsigned int,unsigned int>::const_iterator cit=m_code.find(index);
  if (cit!=m_code.end()) return cit->second;
  return 0;
}

int Flow::Index(const unsigned int code) const
{
  for (std::map<unsigned int,unsigned int>::const_iterator 
	 cit=m_code.begin();cit!=m_code.end();++cit) {
    if (cit->second==code) return cit->first;
  }
  return -1;
}















