#include "Jet_Finder_Base.H"

#include "Message.H"

namespace ATOOLS {

  std::vector<int> ID(size_t id)
  {
    std::vector<int> ids;
    for (size_t n(0);id>0;++n) {
      if (id&(1<<n)) {
	ids.push_back(n);
	id-=1<<n;
      }
    }
    return ids;
  }

  size_t IDCount(size_t id)
  {
    size_t idc(0);
    for (size_t n(0);id>0;++n) {
      if (id&(1<<n)) {
	++idc;
	id-=1<<n;
      }
    }
    return idc;
  }

}

using namespace ATOOLS;

Jet_Finder_Base::Jet_Finder_Base():
  m_ycut(-1.0), m_gycut(-1.0), m_cycut(-1.0), m_gcycut(-1.0), m_delta_r(1.0)
{
}

Jet_Finder_Base::~Jet_Finder_Base() 
{
}

void Jet_Finder_Base::SetDeltaR(double dr) 
{ 
  if (dr<=1.e-6) {
    msg_Error()<<METHOD<<"(): \\delta_R to small, ignore and set to "
	       <<m_delta_r<<"."<<std::endl;
    return;
  }
  m_delta_r=dr; 
}

void Jet_Finder_Base::FillCombinations()
{
}
