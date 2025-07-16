#include "ATOOLS/Phys/Decay_Info.H"

#include "ATOOLS/Phys/Blob.H"
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

  template DecayInfo_Vector &Blob_Data_Base::Get<DecayInfo_Vector>();

  std::ostream & operator<<(std::ostream &s,const DecayInfo_Vector &ds)
  {
    if (ds.empty()) return s<<"{NULL}";
    s<<"{"<<ID(ds.front()->m_id);
    for (size_t i(1);i<ds.size();++i) s<<','<<ID(ds[i]->m_id);
    return s<<"}";
  }

  template <> Blob_Data<DecayInfo_Vector>::~Blob_Data() {}

  template class Blob_Data<DecayInfo_Vector>;

}
