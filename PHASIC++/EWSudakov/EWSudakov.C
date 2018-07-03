#include "PHASIC++/EWSudakov/EWSudakov.H"
#include "ATOOLS/Phys/Flavour.H"

#include <iostream>

using namespace ATOOLS;

namespace PHASIC {

  std::ostream& operator<<(std::ostream& os, const EWSudakov_Log_Type& t)
  {
    switch (t) {
      case EWSudakov_Log_Type::Ls:
        return os << "L(s)";
      case EWSudakov_Log_Type::lZ:
        return os << "l_Z";
      case EWSudakov_Log_Type::lSSC:
        return os << "l_s";
      case EWSudakov_Log_Type::lC:
        return os << "l_C";
    }
  }

  std::ostream& operator<<(std::ostream& os, const Leg_Set& legset)
  {
    os << "leg:kf_code list: { ";
    for (const auto& leg : legset)
      os << leg.first << ":" << Flavour{leg.second} << " ";
    return os << '}';
  }

}

