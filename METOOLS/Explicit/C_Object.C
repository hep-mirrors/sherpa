#include "METOOLS/Explicit/C_Object.H"

using namespace METOOLS;

CObject::~CObject() {}

std::ostream &METOOLS::operator<<(std::ostream &str,const CObject &s)
{
  return str<<'['<<s(0)<<','<<s(1)<<'|'<<s.H()<<","<<s.S()<<']';
}

