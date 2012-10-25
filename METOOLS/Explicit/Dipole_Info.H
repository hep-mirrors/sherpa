#ifndef METOOLS_Explicit_Dipole_Info_H
#define METOOLS_Explicit_Dipole_Info_H

namespace METOOLS {

  class Dipole_Info {
  private:

    int m_mode, m_stat, m_drmode, m_nf;

    double m_amin, m_kappa, m_amax, m_kt2max;

  public:

    inline void SetMode(const int mode)   { m_mode=mode;   }
    inline void SetDRMode(const int mode) { m_drmode=mode; }

    inline void SetStat(const int stat) { m_stat=stat; }

    inline void SetNf(const int nf) { m_nf=nf; }

    inline void SetAMin(const double &amin) { m_amin=amin; }
    inline void SetKappa(const double &kap) { m_kappa=kap; }

    inline void SetAMax(const double &amax) { m_amax=amax; }

    inline void SetKT2Max(const double &kt2max) { m_kt2max=kt2max; }

    inline int Mode() const   { return m_mode;   }
    inline int DRMode() const { return m_drmode; }

    inline int Stat() const { return m_stat; }

    inline int Nf() const { return m_nf; }

    inline double AMin() const  { return m_amin;  }
    inline double Kappa() const { return m_kappa; }

    inline double AMax() const { return m_amax; }

    inline double KT2Max() const { return m_kt2max; }

  };// end of class Dipole_Info

}// end of namespace METOOLS

#endif
