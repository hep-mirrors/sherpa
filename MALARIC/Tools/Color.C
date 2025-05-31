#include "MALARIC/Tools/Color.H"

namespace MALARIC {

  std::ostream &operator<<(std::ostream &s,const Color &c)
  {
    s<<'('<<c.m_i<<','<<c.m_j<<"){"<<c.m_w<<'|'<<c.m_n<<'}';
    return s;
  }

}
