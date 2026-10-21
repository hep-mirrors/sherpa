#ifndef METOOLS_HadronCurrents_Photon_Emission_H
#define METOOLS_HadronCurrents_Photon_Emission_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

namespace METOOLS {

  //////////////////////////////////////////////////////////////////////////
  /// Conventions for a real photon emitted from INSIDE a current.
  ///
  /// Kept in one place because interference between photons emitted by
  /// different currents -- a lepton line and, later, a charged hadron --
  /// is only right if every current agrees on all of: the polarisation
  /// vector, the coupling, the vertex phase and the eikonal factor.
  ///
  /// For an emitting external leg i with momentum p_i, charge Q_i (the leg's
  /// OWN flavour, as Flavour::Charge() gives it, so -Q_f for an
  /// antiparticle) and eta_i = +1 outgoing / -1 incoming, the leg contributes
  ///
  ///   e Q_i [ eta_i (p_i.eps*)/(p_i.k) J0  +  S_i / (2 p_i.k) ]
  ///
  /// where J0 is the photon-less current and S_i is the current's spin term:
  /// the Dirac equation on the external spinor has already removed the
  /// propagator's p-slash + m, leaving one k-slash. This single rule covers
  /// all four leg types (outgoing/incoming, particle/antiparticle): working
  /// each case through with ubar(p) p-slash = m ubar(p), p-slash u = m u,
  /// vbar p-slash = -m vbar and p-slash v = -m v, the four sign patterns
  /// collapse onto the one formula above. VA_0_FF.H spells out S_i. The vertex
  /// -ieQ times the propagator's i gives a real e Q, so there is no extra
  /// phase relative to the photon-less diagram.
  ///
  /// Default coupling: NONE. A current that is asked to radiate before a
  /// coupling has been set throws, rather than quietly radiating with e=0.
  /// Exchange_ME sets it from the model's alpha_QED, or from Photon_Alpha
  /// on the Currents: block.
  //////////////////////////////////////////////////////////////////////////
  class Photon_Emission {
    double m_e;
    bool   m_gaugetest;
  public:
    Photon_Emission(): m_e(-1.), m_gaugetest(false) {}

    void SetAlpha(const double alpha);
    /// Test scaffold: EpsStar returns k itself, so a gauge-invariant
    /// amplitude must vanish. Never on in production.
    inline void SetGaugeTest(const bool on) { m_gaugetest = on; }
    inline bool GaugeTest() const           { return m_gaugetest; }
    double E() const;

    /// k with its energy reset to |k|. A photon from the phase-space
    /// generator is light-like only to rounding; XYZFunc takes a pseudo-leg's
    /// mass from sqrt(p^2), and the spin sum needs it to be exactly zero.
    ATOOLS::Vec4D LightLike(const ATOOLS::Vec4D & k) const;

    /// eps*_lambda(k) for an outgoing photon, lambda = 0, 1 (the two
    /// transverse states of METOOLS::Polarization_Vector) -- or k under the
    /// gauge test.
    ATOOLS::Vec4C EpsStar(const ATOOLS::Vec4D & k,const int lambda) const;

    /// One leg's contribution, as in the class comment. Zero for a neutral
    /// leg.
    ATOOLS::Vec4C Emit(const ATOOLS::Flavour & f,const ATOOLS::Vec4D & p,
                       const bool incoming,const ATOOLS::Vec4D & k,
                       const ATOOLS::Vec4C & eps,const ATOOLS::Vec4C & J0,
                       const ATOOLS::Vec4C & S) const;
  };
}

#endif
