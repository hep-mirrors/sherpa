#ifndef ATOOLS_Phys_Decay_Info_H
#define ATOOLS_Phys_Decay_Info_H

#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {

  struct Decay_Info {
    size_t m_id;
    Flavour m_fl;
    int m_nmax, m_osd;
    inline Decay_Info(const size_t &id,const Flavour &fl,
		      const int nmax,const int osd):
      m_id(id), m_fl(fl), m_nmax(nmax), m_osd(osd) {}
  };// end of struct Decay_Info

  typedef std::vector<Decay_Info> DecayInfo_Vector;

  std::ostream &operator<<(std::ostream &ostr,const Decay_Info &di);

}// end of namespace ATOOLS

#endif
