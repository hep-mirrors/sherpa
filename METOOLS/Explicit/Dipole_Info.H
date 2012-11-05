#ifndef METOOLS_Explicit_Dipole_Info_H
#define METOOLS_Explicit_Dipole_Info_H

namespace METOOLS {

  class Dipole_Info {
  private:

    int m_mode, m_stat, m_drmode, m_nf, m_msv;

    double m_amin, m_kappa, m_amax[4], m_kt2max, m_mu2;

  public:

    inline void SetMode(const int mode)   { m_mode=mode;   }
    inline void SetDRMode(const int mode) { m_drmode=mode; }

    inline void SetStat(const int stat) { m_stat=stat; }

    inline void SetNf(const int nf) { m_nf=nf; }

    inline void SetMassive(const int msv) { m_msv=msv; }

    inline void SetAMin(const double &amin) { m_amin=amin; }
    inline void SetKappa(const double &kap) { m_kappa=kap; }

    inline void SetAMax(const int t,const double &amax) { m_amax[t]=amax; }

    inline void SetKT2Max(const double &kt2max) { m_kt2max=kt2max; }

    inline void SetMu2(const double &mu2) { m_mu2=mu2; }

    inline int Mode() const   { return m_mode;   }
    inline int DRMode() const { return m_drmode; }

    inline int Stat() const { return m_stat; }

    inline int Nf() const { return m_nf; }

    inline int Massive() const { return m_msv; }

    inline double AMin() const  { return m_amin;  }
    inline double Kappa() const { return m_kappa; }

    inline double AMax(const int t) const { return m_amax[t]; }

    inline double KT2Max() const { return m_kt2max; }

    inline double Mu2() const { return m_mu2; }

  };// end of class Dipole_Info

}// end of namespace METOOLS

#endif
