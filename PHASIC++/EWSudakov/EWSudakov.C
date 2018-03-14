#include "PHASIC++/EWSudakov/EWSudakov.H"

#include <iostream>

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
    }
  }

}

