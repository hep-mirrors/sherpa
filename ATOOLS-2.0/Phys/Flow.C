#include "Flow.H"
#include "Particle.H"


namespace ATOOLS {
  long Flow::s_qcd_counter = 600;
}

using namespace ATOOLS;

Flow::Flow(Particle * _owner) : p_owner(_owner) { 
  m_code.insert(std::make_pair<int,int>(1,0)); 
  m_code.insert(std::make_pair<int,int>(2,0)); 
}

Flow::~Flow() { 
  m_code.clear(); 
}


const Particle * Flow::Owner() const { return p_owner; }

int Flow::Code(int _index) {
  int count = m_code.count(_index);
  if (count>0) return m_code[_index];
  return 0;
}

void Flow::SetCode(int _index,int _code) {
  if (_code==-1) _code = ++s_qcd_counter; 
  //  m_code.insert(std::make_pair<int,int>(_index,_code)); 
  //    caused problems with gcc 2.95.3 20010315 !

  m_code[_index]=_code;
}
















