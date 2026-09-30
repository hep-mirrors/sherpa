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

// How a W+W- -> 4f final state is exponentiated. See the WW_Scheme comment in
// YFS_Base::RegisterDefaults.
struct wwscheme {
  enum code {
    flat = 0,   // the four external fermions are the emitters
    pole = 1    // production stage (W's emit) + one decay stage per W
  };
};

/*!
  Where CEEX's O(alpha) VIRTUAL comes from. It used to be decided implicitly -
  CEEX supplied its own if and only if no loop provider happened to exist -
  which is how e+e- -> H l+ l- silently took CEEX's 2 -> 2 vertex and boxes and
  reported a +2200% correction with nothing to say so.

  The three are genuinely different physics, not implementations:

    ceex      CEEX's own hand-coded vertex and boxes, per helicity. Amplitude
              level, which is what CEEX is for - and 2 -> 2 ONLY, because the
              expressions are.
    external  the loop provider's helicity-SUMMED virtual through
              Virtual_ME2_Base, applied to |A|^2. Works for any process today.
              The helicity dependence of the boxes is then averaged over, which
              is exactly the coherence CEEX exists to keep, so this is a
              deliberate trade and not a free one.
    helicity  the provider's virtual resolved BY HELICITY, so amplitude-level
              matching survives at any multiplicity (ACRAIC). Needs
              Loop_Generator: Recola (get_amplitude_r1_rcl, no extra loop
              call); delta^h = A1^h/A0^h - alpha B_4 per helicity at the Born
              point, sum_h |A_1 + delta^h A_0|^2 in the CEEX weight. Falls
              back to v/2 (= external) where the provider cannot resolve a
              helicity. YFS/NLO/NLO_LoopHelicity.C, Ceex_LoopHelicity.C.
    none      no virtual in the CEEX column at all: Born + real only, whatever
              NLO_Part says. This is what `external` silently was until
              2026-09-25 (the value was read and never consumed), kept as a
              named choice so Born+real comparisons stay reproducible.
*/
/*
    auto      (default) external when a loop provider supplies the virtual
              (NLO_Part with V and a Loop_Generator), no virtual otherwise;
              CEEX's own 2 -> 2 virtual only when asked for with ceex.
*/
struct ceexvirt {
  enum code {
    ceex     = 0,
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

struct fixed_order {
  enum code {
    off       = 0,
    full      = 1,
    lo        = 2,
    nlo       = 3,
    nnlo      = 4,
  };
};



std::istream &operator>>(std::istream &str, wwscheme::code &sm);
std::istream &operator>>(std::istream &str, ceexvirt::code &sm);
std::ostream &operator<<(std::ostream &str, const ceexvirt::code &sm);
std::ostream &operator<<(std::ostream &str, const wwscheme::code &sm);
std::istream &operator>>(std::istream &str, submode::code &sm);
std::istream &operator>>(std::istream &str, yfsmode::code &sm);
std::istream &operator>>(std::istream &str, wgt::code &sm);
std::istream &operator>>(std::istream &str, fixed_order::code &sm);
std::ostream &operator<<(std::ostream &str, const yfsmode::code &sm);
std::ostream &operator<<(std::ostream &str, const wgt::code &sm);
std::ostream &operator<<(std::ostream &str, const fixed_order::code &sm);

class YFS_Base {
protected:

  double m_exponent, m_xmax;
  bool   m_on;

  // Matrix elements from external providers (OpenLoops, Recola, ...) are
  // pointlike sQED: they know nothing about the pion form factor, while the
  // Born YFS subtracts against is dressed - by the Comix vertex in Vertex mode,
  // by Single_Process in Factored mode. Subtracting a dressed S~*B from an
  // undressed R leaves a residual (1-|F|^2)*S~*B that diverges as k->0 and
  // drives beta_1 negative wherever |F|^2>1. So dress the external ME here,
  // at the same Q^2 the Born is dressed at. Created on first use, since most
  // YFS_Base instances never need it.

public:
  YFS_Base();
  virtual ~YFS_Base();

  bool On() const   { return m_on;   }

  // |F((sum p_pi)^2)|^2, or 1 if no form factor is configured.
  double ExternalFormFactor(const ATOOLS::Vec4D_Vector &p,
                            const ATOOLS::Flavour_Vector &fl,
                            const size_t &nin=2);

  double m_s,m_irscale;
  double m_v, m_vmax, m_isrcut, m_sp;
  double m_alpha, m_alpi, m_hardmin, m_gp, m_g;
  bool m_CalForm, m_fillblob, m_isr_debug, m_fsr_debug, m_coulomb;
  double m_photonMass, m_oneloop, m_born, m_isrWeight, m_fsrWeight, m_rescale_alpha;
  // m_nbar is deliberately NOT here. ISR and FSR each have their own photon
  // multiplicity and each declares it; a shared base copy meant ISR silently
  // used the base member while FSR shadowed it, which is the arrangement that
  // produced the m_betatWW/m_nbar class of bug. Likewise m_setparticles, which
  // is YFS_Handler's flavour latch and was never touched here, and m_mass2[2],
  // which nothing ever read or wrote.
  double m_deltacut, m_betatWW, m_eps, m_fixed_weight, m_resonace_max;
  double m_rv_test_x, m_rv_test_theta, m_rv_test_phi;
  double m_rv_soft_cut; // soft-photon energy fraction below which RV is skipped
  // Same coarse soft-photon gate for the double-real beta_2^2, applied to each
  // of the two photons: the RR residual vanishes when either photon goes soft,
  // while the ME and its subtraction both diverge like S~(k1)*S~(k2), so the
  // pair is roundoff-dominated long before it is physically relevant. Energy
  // fraction E/sqrt(s); <=0 disables.
  double m_rr_soft_cut;
  // Relative-cancellation threshold for the RV beta_1^1. The RV matrix element
  // m_rv and its subtraction m_rvsub each diverge like S(k)~1/k^2 and cancel to
  // form d1; when |d1| drops below m_rv_cancel_eps*max(|m_rv|,|m_rvsub|) the
  // result is roundoff-dominated (and, being that small, physically
  // negligible), so the photon is skipped. <=0 disables. This is the
  // process-independent replacement for the flat RV_SOFT_CUT energy gate.
  double m_rv_cancel_eps;
  // Matrix-element stability guard for the RV. The one-loop provider is
  // numerically unstable for soft, beam-collinear (ISR) photons, returning a
  // real-virtual ME rv up to ~1e30x larger than the scale its own subtraction
  // predicts (max(|rvsub|,|aB|)); the failures are confined to soft/collinear
  // (verified: hard wide-angle photons are always stable), so skipping them
  // loses negligible physics. Skip when |rv| > m_rv_me_max_ratio*subscale.
  // Recommended O(1e6) - the healthy/unstable ratio gap spans 1e4..1e25.
  // <=0 disables.
  double m_rv_me_max_ratio;

  int  m_order, m_betaorder, m_fullform;
  bool m_formWW, m_ifisub;
  wwscheme::code m_wwscheme;
  bool m_wwonshell, m_wwpoleemission;
  // Emission-side IFI: reweight generated photons by the IF radiation function
  // and drop the IF form factor's cutoff to m_ifiomega. Both or neither.
  bool   m_ifireal;
  double m_ifiomega;
  // Per-photon IF reweight factors outside [m_ifi_rclip, 1/m_ifi_rclip] are
  // dropped, and counted here so the rate is reportable rather than silent.
  double m_ifi_rclip;
  long   m_ifi_clipped;
  // One-shot beta_0/beta_1 dump at a fixed phase-space point, for the
  // number-for-number comparison against KKMC's CEEX. See
  // NLO_Base::CEEXComparePoint().
  int  m_ceex_compare;
  // Scale at which Delta_alpha is evaluated for the vacuum-polarisation
  // double-counting subtraction in Virtual/RealVirtual/VirtualVirtual.
  // 0 = photon virtuality (s' or t) - the scale the loop ME itself uses for the
  //     photon self-energy, and therefore the only choice that actually cancels
  //     it. 1 = the runcard mu_R (SCALES), correct only if SCALES is set to that
  //     same virtuality.
  // Fractional envelope on the EEX estimate of the missing double-virtual
  // (YFS_Handler::GenerateWeight): the band is (1 +/- this) * VV_EEX. Default 1
  // = from "no VV at all" to "twice the estimate". Values > 1 make the down
  // variation flip the sign of the term rather than remove it, which is a claim
  // rather than an envelope.
  double m_vv_approx_unc;
  // When set, fill the in-run RV plateau-scan histograms (RV contribution
  // binned by log10 cancellation ratio and by photon energy fraction) so a
  // single production run yields RV_total(cut) and the cut can be validated.
  int  m_rv_cancel_hist;
  int  m_check_mass_reg, m_check_poles, m_check_real, m_check_rv, m_rv_hard_photon;
  //! YFS: CHECK_REAL_REAL - dump the double real's momenta and ME
  //! (Real_Correction with two photons). Read in RegisterSettings; it was
  //! registered but never read, so the dump ran on uninitialised memory.
  int  m_check_rr = 0;
  int m_rmode, m_qedmode, m_kkmcAngles, m_mingammaN;
  bool m_virt_sub;
  // m_no_born_setting is the card value (YFS: No_Born), a global override.
  // m_no_born is the effective flag: the override OR'd with the per-process
  // request from SetNoBorn(), since one handler is shared by all processes.
  bool m_no_born_setting;
  // The YFS.* named weights carry a full per-order breakdown (Real, Virtual,
  // RealReal, the *_1g/_2g multiplicity truncations, the VV band, ...), which
  // is a lot of columns to write for every event. Off by default: only the
  // headline YFS.LO / YFS.NLO / YFS.NNLO are emitted. Set the card option
  // YFS: NLO_Weight_Breakdown: true to get the individual contributions back.
  bool m_nlo_weight_breakdown;
  bool m_ladder_weights;
  int m_dump_dipoles;
  //! YFS: CHECK_INVARIANTS - run YFS_Handler::CheckInvariants each event.
  bool m_check_invariants;
  int m_check_real_sub, m_check_rr_sub, m_tchannel, m_check_virt_born;
  bool m_no_born, m_no_subtraction, m_noflux, m_massless_sub;
  //! Meaning of the YFS: Check_Real_Sub runcard integer (m_check_real_sub).
  //! See NLO_Checks.C for what each check actually does.
  enum CheckRealSubMode {
    CHECK_REAL_SUB_OFF       = 0,
    CHECK_REAL_SUB_SOFT      = 1, //!< CheckRealSub:          E_gamma -> 0 at fixed angle
    CHECK_REAL_SUB_SCATTER   = 2, //!< RecordSubScatter:      accumulating (angle, energy, residual) scatter
    CHECK_REAL_SUB_COLLINEAR = 3  //!< CheckRealCollinearSub: theta -> 0 at fixed energy
  };
  
  bool m_looptool, m_virtual_only, m_real_only, m_use_model_alpha,m_failcut,m_nlocuts;
  /*!
    CEEX_WEIGHT, resolved: take the O(alpha) correction factor from CEEX
    rather than from the EEX beta expansion (the card's auto follows CEEX).
  */
  bool m_ceex_weight = false;
  /*!
    NLO_FSR_FROM_EVENT: take the final-state photons for the fixed-order real
    correction from the EVENT RECORD (Dipole::GetPhotons(), what CEEX is
    given) rather than from Dipole::GetMEPhotons().

    The two are not the same photons. GetMEPhotons() is captured in
    Dipole::MakeFSR before Dipole::Boost(), GetPhotons() after it, so they sit
    in different frames - measured on one event, the same three FSR photons
    came out as 0.0313/0.00255/11.48 GeV against 0.0151/0.00029/1.166.
    GetMEPhotons() additionally keeps the photons HidePhotons() drops, which
    is its stated purpose (see the note in Dipole.C), so switching source also
    changes which photons are corrected at all - the two effects are not
    separable from here.

    Card default true (YFS_Base::RegisterDefaults).
  */
  bool m_nlo_fsr_from_event = false;
  bool m_isrinital, m_realtool, m_useceex, m_coll_real, m_dim_reg;
  /*
    Which NNLO tools exist. Only NLO_Base assigns these (from its p_realvirt,
    p_realreal and p_vv); every other YFS_Base - YFS_Handler, Real_Correction,
    RealVirtual, ... - is a separate object that never does, and used to read
    whatever the heap held. YFS_Handler's m_vvtool gated the VV_EEX/NNLO_VV*
    columns on that garbage (a bool whose low bit was set turned them off:
    NNLO_VV 767 pb in one run, the nominal 1413 pb in the next). Ask p_nlo
    (NLO_Base::HasVV and friends) instead of reading these outside NLO_Base.
  */
  bool m_realvirt = false, m_rrtool = false, m_vvtool = false;
  bool m_skipNegWeights;
  bool m_nlo_fsr_photons;
  bool m_fb_analysis;
  int m_fb_kf;
  fluxmode::code m_flux_mode;
  bool m_hidephotons, m_eex_virt;
  int m_photon_split, m_int_nlo;
  std::string m_name, m_type, m_debugDIR_ISR, m_debugDIR_FSR, m_debugDIR_NLO;
  Vec4D m_beam1, m_beam2;

  ceexvirt::code m_ceexvirtsrc = ceexvirt::ceex;
  //! MODE was not set and CEEX: 1 implied it; YFS_Handler::SetFlavours
  //! then picks ISRFSR or, without charged final-state particles, ISR
  bool m_mode_from_ceex = false;
  submode::code m_submode;
  yfsmode::code m_mode;
  fixed_order::code m_fixedOrder;


  void RegisterSettings();
  void RegisterDefaults();
  double Eikonal(const Vec4D &k, const Vec4D &p1, const Vec4D &p2);
  double EikonalMassless(const Vec4D &k, const Vec4D &p1, const Vec4D &p2);
  inline bool  HasISR()   const {return (m_mode != yfsmode::off && m_mode != yfsmode::fsr); }
  inline bool  HasFSR()   const {return (m_mode != yfsmode::off && m_mode != yfsmode::isr); }
  inline fixed_order::code FixedOrder()  const {return m_fixedOrder;}

};
}

#endif
