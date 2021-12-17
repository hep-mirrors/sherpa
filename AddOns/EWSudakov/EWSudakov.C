#include "AddOns/EWSudakov/EWSudakov.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Flavour.H"

#include <iostream>

using namespace ATOOLS;

namespace EWSud {

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
