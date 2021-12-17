#ifndef EWSudakov_EWSudakov_H
#define EWSudakov_EWSudakov_H

#include "ATOOLS/Phys/EWSudakov_Info.H"
#include "ATOOLS/Phys/Flavour.H"

#include <array>
#include <set>
#include <unordered_map>

/// definitions for types and enums used in the module
/** Note that this module refers in many places to a publication by
 *  Denner/Pozzorini, which can be found under arXiv:hep-ph/0010201.
 *  Also the definitions and calculations follow their conventions
 *  as far as possible to make comparisons and checks more straightforward.
 */

namespace EWSud {

  using Two_Leg_Indizes = std::array<size_t, 2>;
  using Leg_Index_Set = std::set<size_t>;
  // NOTE: the following map can not be unordered, since it is being used as a
  // key and therefore needs to support comparison operations
  using Leg_Kfcode_Map = std::map<size_t, kf_code>;
  // NOTE: Flavour has a (long int)-conversion operator
  using Leg_Kfcode_Map_Signed = std::unordered_map<size_t, long int>;
  using Couplings = std::unordered_map<long int, Complex>;

  std::ostream& operator<<(std::ostream& os, const Leg_Kfcode_Map&);
  std::ostream& operator<<(std::ostream& os, const Leg_Kfcode_Map_Signed&);

  Leg_Kfcode_Map ConvertToPhysicalPhase(Leg_Kfcode_Map legs);

}

#endif
