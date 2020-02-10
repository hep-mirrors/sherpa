#include "PHASIC++/EWSudakov/EWSudakov.H"
#include "ATOOLS/Phys/Flavour.H"

#include <iostream>

using namespace ATOOLS;

namespace PHASIC {
  EWSudakov_Log_Type convert(const std::string& logt)
  {
      if(logt == "Ls")
        return EWSudakov_Log_Type::Ls ;
      if(logt ==  "lZ")
        return EWSudakov_Log_Type::lZ;
      if(logt ==  "lSSC")
        return EWSudakov_Log_Type::lSSC;
      if(logt ==  "lC")
        return EWSudakov_Log_Type::lC;
      if(logt == "lYuk")
        return EWSudakov_Log_Type::lYuk;
      if(logt == "lPR")
        return EWSudakov_Log_Type::lPR;
  }
  
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
      case EWSudakov_Log_Type::lYuk:
        return os << "l_Yuk";
      case EWSudakov_Log_Type::lPR:
        return os << "l_PR";
    }
  }

  std::ostream& operator<<(std::ostream& os, const Leg_Kfcode_Map& legmap)
  {
    os << "leg:kf_code list: { ";
    for (const auto& leg : legmap)
      os << leg.first << ":" << Flavour{static_cast<long>(leg.second)} << " ";
    return os << '}';
  }

  std::ostream& operator<<(std::ostream& os, const Leg_Kfcode_Map_Signed& legmap)
  {
    os << "leg:signed kf_code list: { ";
    for (const auto& leg : legmap)
      os << leg.first << ":" << Flavour{leg.second} << " ";
    return os << '}';
  }
  
}
