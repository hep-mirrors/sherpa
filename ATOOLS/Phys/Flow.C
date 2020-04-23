#include "ATOOLS/Phys/Flow.H"

using namespace ATOOLS;

unsigned int Flow::s_qcd_counter=600;

namespace ATOOLS {
std::ostream& operator<<(std::ostream &ostr,const Flow &flow)
{
  ostr << "[";
  for (int i = 0; i < 2; i++)
    ostr << "(" << i+1 << "=" << flow.m_codes[i] << ")";
  return ostr << "]";
}
}

Flow::Flow()
{
  for (int i = 0; i < 2; i++)
    m_codes[i] = 0;
}

Flow::Flow(const Flow &flow)
{
  for (int i = 0; i < 2; i++)
    m_codes[i] = flow.m_codes[i];
}

Flow::~Flow() {}

void Flow::SetCode(const unsigned int index,const int code) 
{
  if (code == -1)
    m_codes[index-1] = ++s_qcd_counter;
  else
    m_codes[index-1] = code;
}

void Flow::SetCode(const Flow &flow)
{
  for (int i = 0; i < 2; i++)
    m_codes[i] = flow.m_codes[i];
}

unsigned int Flow::Code(const unsigned int index) const
{
  return m_codes[index-1];
}

int Flow::Index(const unsigned int code) const
{
  for (int i = 0; i < 2; i++)
    if (m_codes[i] == code)
      return i+1;
  return -1;
}

void Flow::SwapColourIndices() {
  unsigned int help(m_codes[0]);
  m_codes[0] = m_codes[1];
  m_codes[1] = help;
}
