#include "ATOOLS/Phys/Color_Tester.H"

using namespace ATOOLS;

Color_Tester::Color_Tester(const unsigned int index,
			   const unsigned int color):
  m_color(std::pair<unsigned int,unsigned int>(index,color)) {}

void Color_Tester::Turn()
{
  m_color.first=3-m_color.first;
}

bool Color_Tester::Test(const Particle *parton) const
{
  if (parton->GetFlow()->Code(m_color.first)==m_color.second) return true;
  return false;
}
