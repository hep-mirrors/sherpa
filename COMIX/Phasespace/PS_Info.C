#include "COMIX/Phasespace/PS_Info.H"

using namespace COMIX;

std::ostream &COMIX::operator<<(std::ostream &str,const PS_Info &s)
{
  return str<<'{'<<s.H()<<';'<<s(0)<<","<<s(1)<<'|'<<s[0]<<'}';
}

