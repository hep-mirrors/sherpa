#ifndef YFS_Main_YFS_Base_H
#define YFS_Main_YFS_Base_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "YFS/Main/YFS_Options.H"
#include <memory>


using namespace ATOOLS;

namespace YFS {
struct submode {
  enum code {
    off        = 0,
    local       = 1,
    global       = 2
  };
};
struct yfsmode {
  enum code {
    off       = 0,
    isr       = 1,
    isrfsr    = 2,
    fsr       = 3,
  };
};

//! How a W+W- -> 4f final state is exponentiated (YFS: WW_Scheme).
struct wwscheme {
  enum code {
    flat = 0,   // the four external fermions are the emitters
    pole = 1    // production stage (W's emit) + one decay stage per W
  };
};

/*!
  Where ACRAIC's O(alpha) virtual comes from.
    acraic    ACRAIC's own vertex and boxes, per helicity (2 -> 2 only)
    external  the provider's helicity-summed virtual applied to |A|^2 (any
              process; box helicity dependence is averaged)
    helicity  the provider's virtual per helicity (Loop_Generator: Recola),
              delta^h = A1^h/A0^h - alpha B_4; v/2 where unresolved
              (NLO_LoopHelicity.C, Acraic_LoopHelicity.C)
    none      no virtual in the ACRAIC column: Born + real only
    auto      (default) external if a loop provider supplies the virtual, else none
*/
struct acraicvirt {
  enum code {
    acraic     = 0,
    external = 1,
    helicity = 2,
    none     = 3,
    automatic = 4
  };
};

struct wgt {
  enum code {
    off       = 0,
    full      = 1,
    mass      = 2,
    hide      = 3,
    jacob     = 4,
  };
};




std::istream &operator>>(std::istream &str, wwscheme::code &sm);
std::istream &operator>>(std::istream &str, acraicvirt::code &sm);
std::ostream &operator<<(std::ostream &str, const acraicvirt::code &sm);
std::ostream &operator<<(std::ostream &str, const wwscheme::code &sm);
std::istream &operator>>(std::istream &str, submode::code &sm);
std::istream &operator>>(std::istream &str, yfsmode::code &sm);
std::istream &operator>>(std::istream &str, wgt::code &sm);
std::ostream &operator<<(std::ostream &str, const yfsmode::code &sm);
std::ostream &operator<<(std::ostream &str, const wgt::code &sm);

class YFS_Base {
protected:

  double m_exponent, m_xmax;
  bool   m_on;

public:
  YFS_Base();
  virtual ~YFS_Base();

  bool On() const   { return m_on;   }

  //! |F((sum p_pi)^2)|^2 dressing an external (pointlike) ME so it matches the
  //! dressed Born; 1 if no form factor is configured.
  double ExternalFormFactor(const ATOOLS::Vec4D_Vector &p,
                            const ATOOLS::Flavour_Vector &fl,
                            const size_t &nin=2);

  //! Set this event's s = (p_beam1 + p_beam2)^2 (YFS_Handler::SetEventS, for
  //! beam spectra). Cuts fixed at construction (m_isrcut, m_vmax) stay; owners
  //! of further YFS objects forward.
  virtual void SetS(const double s) {
    m_s = s;
    if (m_ifiomega_is_emin) m_ifiomega = 0.5*sqrt(m_s)*m_isrcut;
  }

  double m_s,m_irscale;
  double m_v, m_vmax, m_isrcut, m_sp;
  //! YFS: EXCHANGE_T_MIN, GeV^2
  double m_exchange_tmin;
  //! YFS: EXCHANGE_T_MIN_FINAL, GeV^2: the same bound on the physical event
  double m_exchange_tmin_final = 0.;
  //! YFS: BORN_PHOTON_REDUCED_ANGLE, degrees
  double m_bornphoton_reduced_angle = 0.;
  //! YFS: RESOLVED_PHOTON_E [GeV] and RESOLVED_PHOTON_ANGLE [deg]
  double m_resolved_photon_e = 0., m_resolved_photon_angle = 0.;
  double m_alpha, m_alpi, m_gp, m_g;
  bool m_isr_debug, m_fsr_debug, m_coulomb;
  //! YFS: COULOMB_HIGHER_ORDER (flat four-fermion W pairs, see YFS_Options.H)
  coulombho::code m_coulombho;
  double m_photonMass, m_oneloop, m_born, m_isrWeight, m_fsrWeight, m_rescale_alpha;
  double m_deltacut, m_betatWW, m_eps, m_fixed_weight, m_resonace_max;
  double m_rv_test_x, m_rv_test_theta, m_rv_test_phi;
  double m_rv_soft_cut; // soft-photon energy fraction below which RV is skipped
  //! YFS: RV_COLLINEAR_CUT, min_i (p_i.k)/s below which RV is skipped
  double m_rv_collinear_cut = 0.;
  // Soft-photon gate (E/sqrt(s)) for the double-real beta_2^2, per photon; <=0 disables.
  double m_rr_soft_cut;
  // Skip an RV photon when |d1| < m_rv_cancel_eps*max(|m_rv|,|m_rvsub|)
  // (roundoff-dominated cancellation); <=0 disables.
  double m_rv_cancel_eps;
  // Skip an RV photon when |rv| > m_rv_me_max_ratio*max(|rvsub|,|aB|), a guard
  // against unstable soft/collinear provider MEs; O(1e6) recommended, <=0 disables.
  double m_rv_me_max_ratio;

  int  m_order, m_betaorder, m_fullform;
  bool m_formWW, m_ifisub;
  //! YFS: IFI - initial-final interference (false: none, IFI_Sub and IFI_Real too)
  bool m_ifi = true;
  wwscheme::code m_wwscheme;
  bool m_wwonshell, m_wwpoleemission;
  // Emission-side IFI: reweight photons by the IF radiation function and set
  // the IF form factor's cutoff to m_ifiomega. Both or neither.
  bool   m_ifireal;
  double m_ifiomega;
  //! YFS: IF_VIRTUAL - the initial-final dipoles' virtual B
  ifvirtual::code m_ifvirtual = ifvirtual::tbvirt_exact_ir;
  //! IFI_Omega at its default, FSR's m_Emin = sqrt(s)/2 * m_isrcut (follows s via SetS)
  bool   m_ifiomega_is_emin = false;
  // Per-photon IF reweight factors outside [m_ifi_rclip, 1/m_ifi_rclip] are
  // dropped and counted in m_ifi_clipped.
  double m_ifi_rclip;
  long   m_ifi_clipped;
  // One-shot beta_0/beta_1 dump at a fixed point (NLO_Base::ACRAICComparePoint()).
  int  m_acraic_compare;
  // Fractional envelope on the EEX estimate of the missing double-virtual
  // (YFS_Handler::GenerateWeight): the band is (1 +/- this) * VV_EEX.
  double m_vv_approx_unc;
  // Fill the in-run RV plateau-scan histograms (RV vs cancellation ratio and photon energy).
  int  m_rv_cancel_hist;
  int  m_check_mass_reg, m_check_poles, m_check_real, m_check_rv, m_rv_hard_photon;
  //! YFS: CHECK_REAL_REAL - dump the double real's momenta and ME
  int  m_check_rr = 0;
  int m_rmode, m_qedmode, m_kkmcAngles;
  bool m_virt_sub;
  // m_no_born_setting: card value (YFS: No_Born), a global override.
  // m_no_born: effective flag, the override OR'd with the per-process SetNoBorn().
  bool m_no_born_setting;
  // YFS: NLO_Weight_Breakdown: the named YFS.* weights; off by default (nominal
  // weight only), and a stream switched on without it is a fatal error.
  bool m_nlo_weight_breakdown;
  // YFS: NLO_Columns: the YFS.NLO and YFS.NLO_EW named weights, off by default
  bool m_nlo_columns;
  bool m_ladder_weights;
  int m_dump_dipoles;
  int m_check_real_sub, m_check_rr_sub, m_check_virt_born;
  bool m_no_born, m_no_subtraction, m_noflux, m_massless_sub;
  //! Values of YFS: Check_Real_Sub (m_check_real_sub); see NLO_Checks.C.
  enum CheckRealSubMode {
    CHECK_REAL_SUB_OFF       = 0,
    CHECK_REAL_SUB_SOFT      = 1, //!< CheckRealSub:          E_gamma -> 0 at fixed angle
    CHECK_REAL_SUB_SCATTER   = 2, //!< RecordSubScatter:      accumulating (angle, energy, residual) scatter
    CHECK_REAL_SUB_COLLINEAR = 3  //!< CheckRealCollinearSub: theta -> 0 at fixed energy
  };
  
  bool m_looptool, m_virtual_only, m_real_only, m_use_model_alpha,m_failcut,m_nlocuts;
  //! ACRAIC_WEIGHT, resolved: take the O(alpha) factor from ACRAIC rather than
  //! the EEX beta expansion.
  bool m_acraic_weight = false;
  //! NLO_FSR_FROM_EVENT: take the fixed-order real's final-state photons from the
  //! event record (Dipole::GetPhotons()) rather than Dipole::GetMEPhotons(),
  //! which is captured before Dipole::Boost() and keeps hidden photons.
  bool m_nlo_fsr_from_event = false;
  bool m_isrinital, m_realtool, m_useacraic, m_coll_real, m_dim_reg;
  /*
    Which NNLO tools exist. Only NLO_Base assigns these; other YFS_Base objects
    should ask p_nlo (NLO_Base::HasVV and friends).
  */
  bool m_realvirt = false, m_rrtool = false, m_vvtool = false;
  bool m_skipNegWeights;
  bool m_nlo_fsr_photons;
  fluxmode::code m_flux_mode;
  bool m_hidephotons, m_eex_virt;
  int m_photon_split, m_int_nlo;
  std::string m_name, m_type, m_debugDIR_ISR, m_debugDIR_FSR, m_debugDIR_NLO;
  Vec4D m_beam1, m_beam2;

  acraicvirt::code m_acraicvirtsrc = acraicvirt::acraic;
  //! MODE unset and implied by ACRAIC: 1 (YFS_Handler::SetFlavours picks ISRFSR, or ISR)
  bool m_mode_from_acraic = false;
  submode::code m_submode;
  yfsmode::code m_mode;


  void RegisterSettings();
  void RegisterDefaults();
  double Eikonal(const Vec4D &k, const Vec4D &p1, const Vec4D &p2);
  double EikonalMassless(const Vec4D &k, const Vec4D &p1, const Vec4D &p2);
  inline bool  HasISR()   const {return (m_mode != yfsmode::off && m_mode != yfsmode::fsr); }
  inline bool  HasFSR()   const {return (m_mode != yfsmode::off && m_mode != yfsmode::isr); }

};
}

#endif
