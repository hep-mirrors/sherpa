#include "ATOOLS/Phys/Decay_Info.H"

#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <iomanip>

using namespace ATOOLS;

namespace ATOOLS {

  std::ostream &operator<<(std::ostream &ostr,const Decay_Info &di)
  {
    ostr<<ToString(ID(di.m_id))<<"["<<di.m_fl<<"|"
	<<di.m_nmax<<","<<di.m_osd<<"]";
    ostr<<" ("<<&di<<")";
    if (di.m_subsequentdecays.size()>0) {
      ostr<<" -> ";
      for (size_t i(di.m_subsequentdecays.size());i>0;--i) {
        ostr<<" "<<di.m_subsequentdecays[i-1];
      }
    }
    return ostr;
  }

}
