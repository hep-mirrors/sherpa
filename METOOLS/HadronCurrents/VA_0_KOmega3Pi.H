#ifndef METOOLS_Current_Library_VA_0_KOmega3Pi_H
#define METOOLS_Current_Library_VA_0_KOmega3Pi_H

#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include <list>

namespace METOOLS { class Total_Width_Base; }

namespace METOOLS {
  // K-omega -> K pi+ pi- pi0 (tau_two_meson_currents_KS_RChiT.tex,
  // Sec."K omega -> K pi+pi-pi0: four-pseudoscalar channel"). A
  // genuinely NEW, self-contained Current class - NOT built on the
  // FF_0_PPP scalar-form-factor (F1,F2,F3,FS) pattern used by every
  // other 3-meson Current class in this directory, because the
  // physics itself is different in kind: a quasi-two-body K-omega
  // production current (Eq. PV-current) with the omega polarization
  // eliminated via its own propagator, feeding into the anomalous
  // omega->3pi vertex - a genuine intermediate massive-vector
  // resonance sitting between the weak vertex and the final state,
  // not a simple sum of scalar propagators multiplying fixed
  // kinematic structures.
  //
  // Order: 0 = K, 1 = pi^+, 2 = pi^-, 3 = pi^0. Unlike every other
  // channel in this directory, there is NO Bose symmetrization
  // needed here at all - all four daughters are physically distinct
  // particles (K^-pi^+pi^-pi^0), so this exact order is required but
  // not because of any special convention; it is simply unambiguous.
  //
  // Physics content (see class-level comment in VA_0_KOmega3Pi.C for
  // the full derivation): a generic PV production current
  //   H^{mu,rho}(Q^2) = -g(Q^2) eps^{mu,rho,alpha,kappa} p_omega_alpha p_K_kappa
  //                      - i f(Q^2) g^{mu,rho}
  //                      - i{a1(Q^2) p_omega^mu + a2(Q^2) p_K^mu} Q^rho
  // contracted through the omega propagator into the anomalous
  // omega->3pi vertex. This reduces (using the fact that the
  // omega->3pi vector is exactly transverse to the omega's own
  // momentum) to ordinary 4-vector/Levi-Civita arithmetic - no
  // genuine rank-2 tensor code is needed.
  class VA_0_KOmega3Pi : public Current_Base {
  private:
    // --- f,g,a1,a2 resonance content (Eqs. PV-f..PV-a2, PV-DR) ---
    // P_i = K (single state, stable: Gamma_K=0 exactly).
    double  m_MK;
    // A_i = K1(1270), K1(1400).
    double  m_MK1a, m_GK1a, m_MK1b, m_GK1b;
    // V_i = K*(892), K*(1680).
    double  m_MKsta, m_GKsta, m_MKstb, m_GKstb;
    // omega mass (fixed pole mass used in the k^rho k^alpha/M_omega^2
    // term and as the running-width propagator's pole) - width comes
    // from Sherpa's own registered omega(782) running width via
    // p_omega_width below, not a fixed number, per Eq.(PV-vector-prop)
    // explicitly writing Gamma_V(k^2).
    double  m_Momega;
    Total_Width_Base * p_omega_width;
    // rho(770) running width, reused 3x (at s_+-,s_+0,s_-0) for the
    // omega->3pi vertex's VMD realisation, Eq.(omega3pi-vmd).
    Total_Width_Base * p_rho_width;

    // h_R*t_{RomegaK} coupling PRODUCTS [GeV], exposed as their own
    // parameter block. Defaults are the meson-dominance values of
    // arXiv:1408.0086 Tab.I, which that paper applies to K omega as
    // well as K rho - see the extended note at the assignment site in
    // the .C for the source quotation and the units question.
    Complex m_hKtKwK;             // P_i=K
    Complex m_hA1tA1wK, m_hA2tA2wK; // A_i=K1(1270),K1(1400)
    Complex m_hV1tV1wK, m_hV2tV2wK; // V_i=K*(892),K*(1680)
    // N = g_omegarhopi*g_rhopipi for G_omega3pi [GeV^-1]. Not given a
    // number in the note; fixed by calibrating this vertex against the
    // measured Gamma(omega->pi+pi-pi0) - see the assignment site in the
    // .C. The contact term g_omega3pi^ct stays 0 (pure rho-exchange
    // implementation, as the note's own default).
    Complex m_Nomega3pi;

    double  m_norm; // overall Vud/Vus-type CKM + normalisation prefactor

    Complex DRfixed(const double & Q2,const double & M,const double & G) const;
    Complex FF_f(const double & Q2)  const;
    Complex FF_g(const double & Q2)  const;
    Complex FF_a1(const double & Q2) const;
    Complex FF_a2(const double & Q2) const;
    Complex Gomega3pi(const double & spm,const double & sp0,
		      const double & sm0) const;

  public:
    VA_0_KOmega3Pi(const ATOOLS::Flavour_Vector& flavs,
		   const std::vector<int>& indices,
		   const std::string& name);
    ~VA_0_KOmega3Pi();
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
