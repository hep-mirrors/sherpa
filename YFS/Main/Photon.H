#ifndef YFS_Tools_Photon_H
#define YFS_Tools_Photon_H

#include "YFS/Main/Dipole.H"

#include <vector>

namespace YFS {

  class Photon {
  public:
    Photon() = default;
    Photon(const ATOOLS::Vec4D &k, const Dipole *dip) : m_k(k), p_dipole(dip) {}

    // -- kinematics ---------------------------------------------------------
    const ATOOLS::Vec4D &K() const { return m_k; }
    ATOOLS::Vec4D       &K()       { return m_k; }
    double E() const { return m_k[0]; }

    // -- origin -------------------------------------------------------------
    const Dipole *Dip() const { return p_dipole; }
    dipoletype::code Origin() const { return p_dipole->m_type; }
    bool IsISR() const { return p_dipole->m_type == dipoletype::initial; }
    bool IsFSR() const { return p_dipole->m_type == dipoletype::final;   }
    bool IsIF()  const { return p_dipole->m_type == dipoletype::ifi;     }

    // -- generation record --------------------------------------------------
    double m_crude   = 0.;
    double m_scale   = 1.;

    // m_y / m_z (Sudakov) and m_masswgt were here, mirroring ISR::m_yini,
    // m_zini and FSR::m_MassWls. Removed: those never leave their generator,
    // so there is no consumer to serve. Add them back when one appears.

    double m_beta10    = 0.;
    double m_beta11    = 0.;
    double m_ifi_ratio = 1.;

    double beta10() const { return m_beta10; }
    double beta11() const { return m_beta11; }

    // The soft / unstable / zero-amplitude flags belong here, but the
    // counters they would feed (m_softRV, m_rvUnstable, m_zero_real_amp)
    // still live on NLO_Base. Added when those move, not before.

  private:
    ATOOLS::Vec4D  m_k;
    const Dipole  *p_dipole = nullptr;
  };

  typedef std::vector<Photon> Photon_Vector;

  class PhotonPair {
  public:
    PhotonPair(const Photon &k1, const Photon &k2) : p_k1(&k1), p_k2(&k2) {}

    const Photon &First()  const { return *p_k1; }
    const Photon &Second() const { return *p_k2; }

    double m_beta20 = 0.;    //!< this pair's real-real contribution
    double beta20() const { return m_beta20; }

  private:
    const Photon *p_k1 = nullptr;
    const Photon *p_k2 = nullptr;
  };

}// end of namespace YFS

#endif
