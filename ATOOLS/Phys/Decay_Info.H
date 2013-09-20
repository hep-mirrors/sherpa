#ifndef ATOOLS_Phys_Decay_Info_H
#define ATOOLS_Phys_Decay_Info_H

#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {

  struct Decay_Info;

  typedef std::vector<Decay_Info* > DecayInfo_Vector;

  struct Decay_Info {
    size_t m_id;
    Flavour m_fl;
    int m_nmax, m_osd;
    DecayInfo_Vector m_subsequentdecays;
    inline Decay_Info(const size_t &id,const Flavour &fl,
		      const int nmax,const int osd):
      m_id(id), m_fl(fl), m_nmax(nmax), m_osd(osd) {}
    inline void SetSubsequentDecayInfos(DecayInfo_Vector sdi)
    { m_subsequentdecays=sdi; }
    inline const DecayInfo_Vector& SubsequentDecayInfos()
    { return m_subsequentdecays; }
  };// end of struct Decay_Info


  std::ostream &operator<<(std::ostream &ostr,const Decay_Info &di);

  std::ostream &operator<<(std::ostream &s,const DecayInfo_Vector &ds);

}// end of namespace ATOOLS

#endif
