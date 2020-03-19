#include "PHASIC++/EWSudakov/EWSudakov.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Flavour.H"

#include <iostream>

using namespace ATOOLS;

namespace PHASIC {

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
