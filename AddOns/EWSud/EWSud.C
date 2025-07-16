#include "AddOns/EWSud/EWSud.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Flavour.H"

#include <iostream>

using namespace ATOOLS;

namespace EWSud {

  EWSudakov_Log_Type EWSudakovLogTypeFromString(const std::string& logt)
  {
    if (logt == "LSC")
      return EWSudakov_Log_Type::Ls;
    else if (logt == "Z")
      return EWSudakov_Log_Type::lZ;
    else if (logt == "SSC")
      return EWSudakov_Log_Type::lSSC;
    else if (logt == "C")
      return EWSudakov_Log_Type::lC;
    else if (logt == "Yuk")
      return EWSudakov_Log_Type::lYuk;
    else if (logt == "PR")
      return EWSudakov_Log_Type::lPR;
    else if (logt == "I")
      return EWSudakov_Log_Type::lI;
    else
      THROW(fatal_error,
            "Can not convert " + logt + " to EW Sudakov log type.");
  }

  std::ostream& operator<<(std::ostream& os, const EWSudakov_Log_Type& t)
  {
    switch (t) {
      case EWSudakov_Log_Type::Ls:
        return os << "LSC";
      case EWSudakov_Log_Type::lZ:
        return os << "Z";
      case EWSudakov_Log_Type::lSSC:
        return os << "SSC";
      case EWSudakov_Log_Type::lC:
        return os << "C";
      case EWSudakov_Log_Type::lYuk:
        return os << "Yuk";
      case EWSudakov_Log_Type::lPR:
        return os << "PR";
      case EWSudakov_Log_Type::lI:
        return os << "I";
      default:
        return os;
    }
  }

  double EWSudakov_Log_Corrections_Map::KFactor() const
  {
    double kfac = 1.0;
    for (const auto &kv : *this) {
      kfac += kv.second;
    }
    return kfac;
  }


  std::ostream& operator<<(std::ostream& o,
                           const EWSudakov_Log_Corrections_Map& m)
  {
    o << "1 - K_EWSud = " << m.KFactor() << " (";
    bool is_first {true};
    for (const auto &kv : m) {
      o << kv.first << ": " << (is_first ? "" : ", ") << kv.second;
      is_first = false;
    }
    return o << ')';
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

  Leg_Kfcode_Map ConvertToPhysicalPhase(Leg_Kfcode_Map legs) {
    for (auto& kv : legs) {
      if (kv.second == kf_phiplus)
        kv.second = kf_Wplus;
      else if (kv.second == kf_chi)
        kv.second = kf_Z;
    }
    return legs;
  }

}
