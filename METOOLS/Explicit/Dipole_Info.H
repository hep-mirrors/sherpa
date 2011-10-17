#ifndef METOOLS_Explicit_Dipole_Info_H
#define METOOLS_Explicit_Dipole_Info_H

namespace METOOLS {

  class Dipole_Info {
  private:

    int m_mode, m_stat, m_drmode;

    double m_amin, m_amax;

  public:

    inline void SetMode(const int mode)   { m_mode=mode;   }
    inline void SetDRMode(const int mode) { m_drmode=mode; }

    inline void SetStat(const int stat) { m_stat=stat; }

    inline void SetAMin(const double &amin) { m_amin=amin; }
    inline void SetAMax(const double &amax) { m_amax=amax; }

    inline int Mode() const   { return m_mode;   }
    inline int DRMode() const { return m_drmode; }

    inline int Stat() const { return m_stat; }

    inline double AMin() const { return m_amin; }
    inline double AMax() const { return m_amax; }

  };// end of class Dipole_Info

}// end of namespace METOOLS

#endif
