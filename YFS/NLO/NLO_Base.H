#ifndef YFS_NLO_BASE_H
#define YFS_NLO_BASE_H

#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/FSR.H"
#include "YFS/Main/Leg_Set.H"
#include "YFS/NLO/Virtual.H"
#include "YFS/NLO/Real_Correction.H"
#include "YFS/NLO/RealVirtual.H"
#include "YFS/NLO/VirtualVirtual.H"
#include "YFS/Main/YFS_Form_Factor.H"
#include "YFS/Main/Photon.H"
#include "YFS/Main/Define_Dipoles.H"
#include "YFS/Main/Beta0_2.H"
#include "YFS/Main/Beta0_2_Boxes.H"
#include "YFS/Main/Dyson_VP.H"
#include "ATOOLS/Math/Histogram_2D.H"
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <memory>
#include "ATOOLS/Math/Histogram.H"

#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "PHASIC++/Channels/Multi_Channel.H"




namespace YFS {

typedef std::map<std::string, std::unique_ptr<ATOOLS::Histogram_2D> > histograms2d;
typedef std::map<std::string, std::unique_ptr<ATOOLS::Histogram> > histograms1d;

/*!
  ACRAIC: REAL_VIRTUAL helicity - one RV photon's per-helicity remainder,
  from the (n+1)-point loop call of RV_MODE remainder
  (NLO_Base::BuildRealVirtualHelicityFactors), read by ACRAIC
  (Acraic_Base::RealVirtualHelicityToAcraic).
*/
struct RVHelicityFactors {
  ATOOLS::Vec4D k;                          //!< the lab photon
  double dv = 0.;                           //!< v_{n+1} + c - v_B (squared level)
  bool ok = false;                          //!< rem/a0/hel are filled
  std::vector<std::vector<int> > hel;       //!< provider labels, photon last
  std::vector<Complex> a0;                  //!< provider tree amplitudes
  std::vector<Complex> rem;                 //!< per-helicity remainder
  ATOOLS::Vec4D_Vector plab;                //!< the (n+1)-point in the lab
};

class NLO_Base : YFS_Base {
public:
  //! This event's s, forwarded to the dipoles and form factor NLO_Base owns.
  void SetS(const double s) override;
  // subloc/subb seen by the IFI_Real restoration, accumulated so the tail that
  // drives the MC error is reportable instead of inferred from sigma.
  long   m_ifi_n = 0;
  double m_ifi_sum = 0., m_ifi_min = 1e300, m_ifi_max = -1e300, m_ifi_sum2 = 0.;
  // prod_i subloc(k_i)/subb(k_i) over one event's photons.
  double m_ifi_prod = 1.;
  // B (m_ifi_prod - 1), left by the last CalculateReal(): what the IF ratios
  // add to the real sum (0 without IFI_Real). The breakdown columns that are
  // not built from the full real sum add it back, see YFS_Handler.
  double m_real_ifi_soft = 0.;
  // Photons of the real sum above IFIOmega(), and those that never reached
  // IFIRealRatio, so their IF factor stayed 1: below PhotonEminNLO(), failed
  // momentum conservation, a zero or bad real ME, or a ratio IFIRealRatio
  // refused.
  long m_ifi_eligible = 0, m_ifi_skip_emin = 0, m_ifi_skip_momcons = 0,
       m_ifi_skip_unusable = 0, m_ifi_skip_other = 0;
  // Profiled in photon energy fraction x = E_gamma/sqrt(s), 5 bins of 0.1.
  // m_ifi_x_r  : the reweight actually applied, 1 + S_IF/S_crude at BORN momenta
  // m_ifi_x_e  : the "exact-cancellation" ratio subloc/subb, at MAPPED momenta
  // Their difference per bin IS the residue left by using the consistent pair
  // instead of the cancelling one, and it should grow with x if that is what
  // distorts the E_gamma spectrum above the cut cliff.
  long   m_ifi_x_n[5] = {0,0,0,0,0};
  double m_ifi_x_r[5] = {0.,0.,0.,0.,0.};
  double m_ifi_x_e[5] = {0.,0.,0.,0.,0.};


  NLO_Base();
  ~NLO_Base();

  /*!
    Not owned: the fixed-order providers belong to the YFS_Process they were
    built for (see SetProviders), since one NLO_Base is shared by every
    process in the run card.
  */
  YFS::Virtual     *p_virt     = nullptr;
  YFS::Real        *p_real     = nullptr;
  YFS::RealReal    *p_realreal = nullptr;
  YFS::RealVirtual *p_realvirt = nullptr;
  YFS::VirtualVirtual *p_vv    = nullptr;
  std::unique_ptr<YFS::Define_Dipoles> p_nlodipoles;
  std::unique_ptr<YFS::YFS_Form_Factor> p_yfsFormFact;
  //! Not owned: both belong to the YFS_Handler.
  YFS::FSR  *p_fsr = nullptr;
  YFS::Define_Dipoles *p_dipoles = nullptr;
  //! Fixed-order real providers indexed by photon multiplicity; [0] unused,
  //! [1] aliases p_real and [2] p_realreal. Anything above that exists only if
  //! the runcard asked for it via NLO_MAX_PHOTONS.
  std::vector<YFS::Real_Correction*> m_realprov;
  Vec4D_Vector m_plab, m_reallab, m_FSRPhotons, m_ISRPhotons, m_bornMomenta;
  // m_ISRPhotons ++ m_FSRPhotons, each carrying the dipole it came from.
  YFS::Photon_Vector m_photons;

  //! GeV threshold below which a photon gets no fixed-order real correction,
  //! or -1 when the guard is disabled. See NLO_PHOTON_EMIN.
  double PhotonEminNLO() const;
  ATOOLS::Flavour_Vector m_flavs;
  std::ofstream m_subscatter;  // accumulating sub-check angle/energy scatter
  histograms2d m_histograms2d;
  histograms1d m_histograms1d;
  ATOOLS::Momenta_Stretcher m_stretcher;

  double m_evts, m_recola_evts, m_ranTheta, m_ranPhi, m_real,m_ravg, m_rv, m_rvsub, m_eikeex, m_subloc;
  // YFS: REAL_TRACE. Per-photon record of the Born+real assembly for the
  // current event, printed by CalculateReal() when |1 + real/Born| exceeds
  // the threshold. Filled by CalculateReal(k), so it sees every term.
  std::ostringstream m_realtrace;
  // Per-nominal-call contributions of the hardest photon(s), captured
  // inline while CalculateReal/CalculateRealVirtual/CalculateRealReal sum
  // over the full photon multiplicity - no separate recomputation needed.
  // *_hard1 = contribution of the single hardest photon; *_hard2 = sum over
  // the two hardest photons (m_rr_hard2 is the single pair formed by them).
  // Used to report the matching truncated to 1 or 2 photons.
  double m_real_hard1, m_rv_hard1, m_rr_hard2;
  /*
    The single real's final per-photon conventions, left by the last
    CalculateReal(k, raw) call: flux, subtraction eikonal and the denominator
    (the crude) of tot = (r flux - subloc B/resc)/denominator. The double real
    reads them so that its beta_2 subtraction uses exactly the beta_1 the
    single real uses (REAL_SUB_EIK, REAL_FSR_FLUX, REAL_FSR_MAP included).
  */
  double m_lastflux = 1., m_lastsubloc = 0., m_lastdenom = 0.;
  //! the same two before any REAL_SUB_EIK / REAL_FSR_FLUX modification
  double m_lastsubloc0 = 0., m_lastdenom0 = 0.;
  double m_real_hard2, m_rv_hard2;
  // Product of the two per-photon crude eikonal factors S~(k1)*S~(k2) by which
  // the last CalculateRealReal(k1,k2) divided its result (line "tot = .../sub1/sub2").
  // Recorded alongside the RR residual so the 1/S~ blow-up for hard wide-angle
  // photons can be separated from the physical (eikonal-weighted) contribution.
  double m_rr_eik;
  int m_zeroV, m_zeroRV, m_zeroRR, m_nonZeroRR, m_nonZeroRV;
  // Counts photons whose RV beta_1^1 was skipped as numerically unreliable:
  // either softer than the coarse RV_SOFT_CUT energy pre-filter, or (primary
  // gate) failing the RV_CANCEL_EPS relative-cancellation test. Their
  // contribution is negligible but roundoff-dominated.
  int m_softRV;
  // Counts double-real photon pairs skipped by the RR_SOFT_CUT energy
  // pre-filter, i.e. where at least one of the two photons was too soft for
  // the beta_2^2 residual to be computed reliably.
  int m_softRR;
  // Counts RV photons skipped by the matrix-element stability guard
  // (RV_ME_MAX_RATIO): |rv| grossly exceeds the scale its subtraction predicts,
  // i.e. the one-loop provider is unstable. Always active in production when
  // RV_ME_MAX_RATIO > 0 (independent of the RV_CANCEL_HIST diagnostics).
  int m_rvUnstable;
  // RV blow-up mechanism diagnostics, filled only when RV_CANCEL_HIST is set.
  // m_rvHiC: photons with cancellation C>=1 (subtraction not cancelling).
  // m_rvBlowup: photons with grossly unphysical |tot| (> 1e3*|Born|).
  // m_rvBlowupRtree0/HiC/Soft: of those, how many had a zero real-emission ME
  // (rtree==0), C>=1, or were soft (E/sqrt(s)<0.01) - to test whether the
  // blow-ups are exactly the rtree==0 soft photons.
  // m_rvBlowupHardWide: of the blow-ups, how many were hard AND wide-angle
  // (E/sqrt(s)>0.1, no charged leg within ~26deg). If this stays 0, the loop-ME
  // instability is confined to soft/collinear photons and skipping them is safe.
  int m_rvHiC, m_rvBlowup, m_rvBlowupRtree0, m_rvBlowupHiC, m_rvBlowupSoft,
      m_rvBlowupHardWide;
  // Counts CalculateReal(k,raw,...) calls where the raw real-emission
  // amplitude p_real->Calc_R(p) itself came back zero, before any cut or
  // subtraction logic runs - covers both the standalone NLO real correction
  // and the real1/real2 sub-calls embedded in CalculateRealReal.
  int m_zero_real_amp;
  bool m_failcut;


  void Init(Flavour_Vector &flav, Vec4D_Vector &plab, Vec4D_Vector &born);
  /*!
   * Point this (process-agnostic) handler at the ME providers of the process
   * that is currently selected. The providers are owned by the YFS_Process
   * they were built for -- a single YFS_Handler, and hence a single NLO_Base,
   * is shared by every process in the run card, so anything built from a
   * Process_Info must not be stored here permanently. Call this whenever the
   * active process changes, alongside YFS_Handler::SetFlavours.
   * Null arguments are allowed and simply mean "this process has no such
   * correction". Ownership is never taken.
   */
  double m_acraicvirt = 1.;
  void SetProviders(YFS::Virtual *virt, YFS::Real *real,
                    YFS::RealVirtual *realvirt, YFS::RealReal *realreal,
                    YFS::VirtualVirtual *vv);
  /*!
    Whether a loop-ME provider has to be built at all.

    Three things can supply the O(alpha) virtual and only one of them needs an
    external loop ME:
      - an external provider (Loop_Generator: OpenLoops, ...),
      - the EEX collinear virtual (YFS: Collinear_Virtual),
      - ACRAIC's beta_0^1 - the two vertex form factors and the IFI boxes -
        which is internal and exact, so with ACRAIC on and no Loop_Generator
        named there is nothing to look up. Without this, NLO_Part: BVR with
        ACRAIC: 1 and no Loop_Generator died on "Couldn't find virtual ME for
        this process" while holding a perfectly good virtual of its own.

    \param ext true if the card named a real external loop generator. When it
    did, the provider is built even with ACRAIC on: both are then evaluated and
    ACRAIC goes out as its own named weight rather than replacing anything.
  */
  inline bool NeedsVirtualProvider(bool ext=true) const
  { return !m_eex_virt && (ext || !m_useacraic); }
  /// whether ACRAIC is the thing supplying the virtual on this run
  /*!
    Whether ACRAIC's OWN virtual is used, rather than the provider's.

    The condition used to be `!m_looptool` - ACRAIC supplied it exactly when no
    loop provider happened to exist. That is an accident of the card, not a
    choice, and it is how e+e- -> H l+ l- came to use ACRAIC's 2 -> 2 vertex and
    boxes: no V in NLO_Part, so no loop tool, so ACRAIC's own, reporting +2200%
    with nothing to flag it. The source is now stated
    (YFS: ACRAIC_Virtual), and `acraic` still needs a provider to be absent,
    since the two would otherwise both be applied.
  */
  inline bool AcraicSuppliesVirtual() const
  { return m_useacraic && !m_eex_virt && !m_looptool
           && m_acraicvirtsrc == acraicvirt::acraic; }
  /*!
    The O(alpha) virtual factor from ACRAIC, rho(Born+virtual)/rho(Born), pushed
    in by YFS_Handler after MakeACRAIC() has run. Only consulted when
    AcraicSuppliesVirtual(); it is 1 until ACRAIC has produced a number, so an
    event where ACRAIC bailed contributes no virtual rather than a garbage one.
  */
  inline void SetAcraicVirtual(double f) { m_acraicvirt = f; }
  /*!
    The two halves of the external virtual, kept so ACRAIC_Compare can test the
    identity that has to hold if the two are in the same scheme:

        virt/Born - sub  ==  2 Re(deltI + deltF)

    virt is the RAW one-loop from the provider, sub the YFS virtual form
    factor B that is subtracted out of it and put into exp(Y) instead. ACRAIC's
    deltI/deltF are the same remainder by construction (KKMC's F1ini/F1fin are
    "IR subtracted"), so a mismatch localises to one side of that equation
    rather than to "the virtual disagrees".
  */
  //! NB not m_virt_sub - that name is already the Virtual_Sub SETTING on
  //! YFS_Base, and shadowing it here silently switched the subtraction off.
  double m_virt_raw = 0., m_virt_subval = 0.;
  //! diagnostic (YFS: ACRAIC_Compare): m_virt_subval split by dipole type,
  //! evaluated in CalculateVirtual on the same dipoles as sub
  double m_virt_sub_ii = 0., m_virt_sub_ff = 0., m_virt_sub_if = 0.;
  //! diagnostic (YFS: ACRAIC_Compare): the last photon's real terms, from
  //! CalculateReal - |M_1|^2 in the Born's normalisation, flux, subtraction
  //! eikonal, crude density, beta_1/S~
  double m_dbg_r = 0., m_dbg_flux = 0., m_dbg_subloc = 0., m_dbg_subb = 0., m_dbg_tot = 0.;
  double m_dbg_dk = 0., m_dbg_dp = 0.;   //!< how far the point moved the photon / the legs
  //! FSR photon (ACRAIC_Compare): X^2 of its M_1 point (final legs + photon) and of the generator (pre-FSR system)
  double m_dbg_x2pt = 0., m_dbg_x2gen = 0.;
  //! ISR photon (ACRAIC_Compare): invariant mass^2 of the final legs at its M_1 point, and its factor
  double m_dbg_x2isr = 0., m_dbg_facisr = 0., m_dbg_facfsr = 0.;
  /*!
    YFS: ACRAIC_Virtual: helicity. The IR-finite, complex virtual factor per
    helicity configuration of the loop provider, at the Born point m_plab
    where CalculateVirtual evaluated the loop:

      delta^h = A1^h/A0^h - run/2 - alpha B_4

    with A0, A1 the provider's (Recola's) tree and one-loop amplitudes at its
    mu_IR (Delta_IR = 0), run = the uniform vacuum-polarisation correction
    Virtual::Calc_V subtracts (run_corr/B), Re 2 alpha B_4 = the YFS virtual
    subtraction YFS.NLO actually subtracted (sub/m_rescale_alpha, whose
    lambda <-> mu_IR matches the loop's), and Im alpha B_4 fixed so that the
    common imaginary part of delta is KKMC's IR-subtracted form factors (see
    LoopHelicityIR). By construction the Born-weighted average of
    2 Re delta^h is v = m_oneloop/m_born, the value ACRAIC_Virtual: external
    uses; with delta^h = v/2 the helicity form IS external.

    Labels are the provider's (physical helicity per leg); Acraic_Base maps
    them to its own index (Acraic_Base::LoopHelicityToAcraic).
  */
  bool m_lhelok = false;
  std::vector<std::vector<int> > m_lhel;
  std::vector<Complex> m_la0, m_ldelta;
  ATOOLS::Vec4D_Vector m_lhelmom;
  double m_lhelv = 0.;
  //! the raw per-helicity loop over tree A1^h/A0^h of the Born point (0 and
  //! m_lratio_ok false where the tree is below VIRT_HEL_MINBORN), and the
  //! uniform running-alpha correction and Im part delta^h was built with:
  //! the helicity real-virtual rebuilds the Born's delta^h on the Born
  //! point's own legs from these (BuildRealVirtualHelicityFactors)
  std::vector<Complex> m_lratio;
  std::vector<char> m_lratio_ok;
  double m_lrun = 0., m_limsh = 0.;
  //! run statistics of the per-event identity check and the fallbacks
  long   m_lhel_n = 0, m_lhel_fail = 0, m_lhel_nfall = 0, m_lhel_ncfg = 0;
  double m_lhel_res_max = 0., m_lhel_res_sum = 0., m_lhel_res_sq = 0.;
  double m_lhel_clo_max = 0., m_lhel_ir_sum = 0., m_lhel_ir_sq = 0.;
  double m_lhel_ir_max = 0.;
  long   m_lhel_reshist[20] = {0};
  void BuildLoopHelicityFactors(double virt, double sub);
  /*!
    The helicity-independent IR coefficient of the one-loop amplitude,
      K = (alpha/2pi) sum_{i<j} (-Z_i Z_j) [CnuA(s_ij, m_i, m_j) - 1],
    Z = theta Q (theta = +1 incoming, -1 outgoing), s_ij = (p_i +- p_j)^2
    (+ for two legs on the same side). d(A1^h/A0^h)/d ln mu_IR^2 = K for
    every h (verified against Recola 2.3.0 to 1e-7). imsh is the lambda-free
    imaginary part of K ln(mu^2/lambda^2) + alpha B_4(lambda),
      imsh = sum_{i<j} Im K_ij ln(mu_IR^2 e^c / |s_ij|),
    c = ACRAIC: VIRT_HEL_IM_C (1/2 gives KKMC's Im F1 for a vertex pair).
    False if a charged leg is massless.
  */
  bool LoopHelicityIR(const ATOOLS::Vec4D_Vector &p, double alpha, double mu2,
                      Complex &K, double &imsh) const;
  //! the same over the legs flavs (the (n+1)-point: its photon is neutral
  //! and drops out)
  bool LoopHelicityIR(const ATOOLS::Flavour_Vector &flavs,
                      const ATOOLS::Vec4D_Vector &p, double alpha, double mu2,
                      Complex &K, double &imsh) const;
  //! photon-mass YFS virtual subtraction (all dipoles), for the IR diagnostic
  double VirtualSubPhotonMass();
  void ReportLoopHelicity();
  /// whether a real-ME provider should be built at all (collinear approx replaces it)
  inline bool NeedsRealProvider() const { return !m_coll_real; }
  double CalculateVirtual();
  double CalculateReal();
  double CalculateRealReal();
  // raw: return the unsubtracted value (real-real subtracts for itself).
  double CalculateReal(Vec4D k, bool raw=false);
  double CalculateRealReal(Vec4D k1, Vec4D k2);

  // ---- general-n real corrections -------------------------------------
  //
  // beta_n is defined by the YFS subset recursion
  //
  //   beta_n(k_1..k_n) = M_n * flux
  //                    - sum over non-empty S subset of {1..n} of
  //                      [ prod_{i in S} Stilde(k_i) ] * beta_{n-|S|}({k}\S)
  //
  // with beta_0 = Born. Writing it once is the only way to reach n=3 without
  // hand-deriving the subtraction again: the n=1 and n=2 forms below are the
  // |mask|=1 and |mask|=2 cases of exactly this.

  //! The real-correction provider for a photon multiplicity, or NULL if the
  //! run was not set up with one. Index 0 is unused.
  YFS::Real_Correction *RealProvider(size_t nphotons) const;
  //! Register the provider for `nphotons`, growing the table as needed.
  void SetRealProvider(size_t nphotons, YFS::Real_Correction *prov);
  //! Highest photon multiplicity a provider exists for (0 if none).
  size_t MaxRealPhotons() const;
  //! Runcard ceiling on the fixed-order photon multiplicity, NLO_MAX_PHOTONS.
  //! Defaults to 2, so a run that does not ask for more behaves as before.
  static size_t RequestedMaxRealPhotons();
  //! Highest multiplicity to BUILD a provider for: NLO_MAX_PHOTONS, or
  //! ACRAIC: TREE_BETA_MAX when ACRAIC is on at ORDER 2 and asks for more. The
  //! fixed-order R_n of YFS.NLO still stops at RequestedMaxRealPhotons().
  static size_t RequestedMaxProviderPhotons();

  /*!
    Matrix element times flux for the photons selected by `mask`, evaluated in
    the mapping that multiplicity induces, together with the LOCAL eikonals
    Stilde in that same frame (written into `subloc`, indexed as `photons`).
    Returns false if the point is unusable.
  */
  bool RealMEForSubset(unsigned mask, double &me, std::vector<double> &subloc);

  /*!
    beta_|mask| by the recursion above, memoised over subsets.

    Memoisation is not an optimisation detail: without it the n=3 residual
    re-evaluates every beta_1 three times and every beta_2 once per superset,
    and the cost of the truncation policy the photon loops need is set by how
    many DISTINCT subsets get evaluated, not how many times they are asked for.
  */
  double RawBeta(unsigned mask);

  /*!
    Raw residuals this event's photons have already had computed, keyed by a
    bitmask over m_photons.

    The key is over m_photons and the cache lives for the event, not the call,
    so beta_1(k) is computed ONCE however many subsets contain k, and every
    beta_2(i,j) is reused by each triple built on it. Keying over a local
    subset instead - as a per-call memo must - gives the same photon a
    different key in every subset and shares nothing.

    This is the raw (unnormalised) residual, which is what the subtraction
    consumes. It is deliberately separate from Photon::m_beta10 and
    PhotonPair::m_beta20, which hold the FINAL normalised contributions for
    HardestBetas and the hard-photon reporting.
  */
  std::unordered_map<unsigned, double> m_rawbeta;

  //! beta_n / (crude eikonal product) for one photon set - the weight the
  //! caller adds to the correction, for any n.
  //! beta_n / (crude eikonal product) for the photons selected by `mask`
  //! out of m_photons - the weight the caller adds, for any n.
  double CalculateRealN(unsigned mask);

  //! Sum beta_n over every n-photon subset of this event's photons. n=2
  //! reproduces the pair loop in CalculateRealReal(); n=3 needs no new loop,
  //! only a three-photon ME provider.
  double CalculateRealMultiplicity(size_t n);
  double CalculateRealVirtual();
  double CalculateRealVirtual(Vec4D k);
  /*
    YFS: RV_MODE 1 (NLO_Base.C, CalculateRealVirtualRemainder): the real-
    virtual as the non-factorisable remainder on top of the real, evaluated
    on the real's own (n+1)-body point with the real's own density, flux and
    coupling. CalculateReal(k) leaves the point and the photon's exact-over-
    crude ratio in m_lastrvinfo; CalculateReal() stores one per photon of
    m_photons in m_rvinfo (only when RV_MODE is on).
  */
  struct RVPointInfo {
    bool ok = false;
    ATOOLS::Vec4D_Vector p;   //!< the (n+1)-body point, photon last
    double rho = 0.;          //!< kappa r flux/(den B): the photon's exact/crude
  };
  std::vector<RVPointInfo> m_rvinfo;
  RVPointInfo m_lastrvinfo;
  //! v = V_sub/B of this event (CalculateVirtual's return over the Born)
  bool m_vborn_ok = false;
  double m_vborn = 0.;
  //! EW_HIGHER_ORDER's part of m_oneloop/m_born (0 without it)
  double m_vewho = 0.;
  //! BETA0_2: the event's photonic factors and the two columns' settings
  YFS::Photonic_Factors m_b02;
  YFS::Beta02_Config m_b02acraic, m_b02nlo;
  long m_b02_n = 0, m_b02_fail = 0;
  std::unique_ptr<YFS::Dyson_VP> p_dysonvp;
  YFS::VP_Factors m_vp;
  long m_vp_n = 0, m_vp_fail = 0;
  YFS::ewhohel::code  m_ewhohel = YFS::ewhohel::none;
  YFS::m1square::code m_m1sq = YFS::m1square::averaged;
  long m_ghel_n = 0, m_ghel_fail = 0;
  mutable long m_m1sq_hel = 0, m_m1sq_avg = 0;
  int m_rvPoleFail = 0, m_rvNoVirt = 0;
  long m_rvLoopChecked = 0, m_rvLoopUnstable = 0;
  //! YFS: RV_BEAM_CONE: photons whose real-virtual was skipped for lying
  //! within the cone of a beam (NLO_RealVirtual.C, RVInBeamCone)
  long m_rvBeamCone = 0;
  bool RVInBeamCone(const ATOOLS::Vec4D_Vector &p) const;
  //! YFS: VIRTUAL_CHECK - Born-level loops checked in a second frame / dropped
  long m_vLoopChecked = 0, m_vLoopUnstable = 0;
  double m_vLoopUnstableBorn = 0., m_vLoopBornSum = 0.;
  //! v with the subtraction on the Born point's own legs (RV_MODE 1), and
  //! that subtraction itself (loop coupling units)
  bool m_vborn_own_ok = false;
  double m_vborn_own = 0., m_vborn_own_sub = 0.;
  bool BornVirtualOnOwnLegs(double &vb);
  //! RV_PHOTON_CT: the external photon's alpha(0) counterterm on v_{n+1}
  bool m_rvct_done = false;
  double m_rvct = 0.;
  void CalibratePhotonCounterterm();
  //! RV_PHOTON_CT 2: the counterterm from the model's parameters
  //! (Photon_Counterterm.C), per external photon, in the real-virtual
  //! provider's convention (OpenLoops' or Recola's); cim receives the
  //! amplitude-level imaginary part (2 Im c_amp/2), used only per helicity
  double AnalyticPhotonCounterterm(double &cim) const;
  double m_rvct_im = 0.;
  //! RV_MODE 1 with a real-virtual provider (YFS_Handler keeps v x real)
  bool RVRemainder() const;
  //! RV_MODE 1: (lab photon, v_{n+1} - v_B) for every photon whose RV was
  //! evaluated this event - read by ACRAIC: REAL_VIRTUAL 2
  std::vector<std::pair<ATOOLS::Vec4D, double> > m_rvdv;
  //! the last RealVirtualFactor's loop over tree, before its subtraction,
  //! that subtraction (B~/kappa on the point's legs), and the point in the
  //! loop frame
  double m_rvLastLt = 0., m_rvLastSub = 0.;
  ATOOLS::Vec4D_Vector m_rvLastLoop;
  /*!
    ACRAIC: REAL_VIRTUAL helicity. The (n+1)-point provider's per-helicity
    amplitudes, read after photon i's RealVirtualFactor call (the same loop
    call, no extra cost), as the per-helicity remainder of RVHelicityFactors:
        delta_{n+1}^c = A1^c/A0^c - (1/2) B~(point legs)/kappa
                        - i imsh(point) + c_photon/2 ,
        delta_B^c'    = A1^c'/A0^c' (Born point) - run/2
                        - (1/2) B~(Born's own legs)/kappa - i imsh(Born) ,
        rem^c = delta_{n+1}^c - delta_B^c' ,
    c' the Born configuration with the same fermion labels (photon last).
    The (n+1)-tree-weighted average of 2 Re delta_{n+1}^c is v_{n+1} + c by
    construction (checked per photon, m_rvhel_res*); the same average of
    2 Re rem^c differs from the squared-level dv by v_B -
    <2 Re delta_B>_{(n+1) tree}, the Born virtual's helicity dependence
    resolved on the real's helicity mix (m_rvhel_dv*, -> 0 in the soft
    limit). Configurations below VIRT_HEL_MINBORN keep rem = dv/2; a Born
    configuration below it keeps delta_B = v_B/2.
  */
  std::vector<RVHelicityFactors> m_rvhel;
  bool RealVirtualHelicityRequested() const;
  bool BuildRealVirtualHelicityFactors(const Vec4D &k, double dv,
                                       const ATOOLS::Vec4D_Vector &plab,
                                       RVHelicityFactors &out);
  long   m_rvhel_n = 0, m_rvhel_fail = 0, m_rvhel_nfall = 0, m_rvhel_ncfg = 0;
  long   m_rvhel_nbornfall = 0;
  double m_rvhel_res_max = 0., m_rvhel_res_sum = 0., m_rvhel_res_sq = 0.;
  double m_rvhel_dv_max = 0., m_rvhel_dv_sum = 0., m_rvhel_dv_sq = 0.;
  void ReportRealVirtualHelicity();
  //! each photon's factor in REAL_COMBINE's product (1 if not evaluated)
  std::vector<double> m_prodfac;

  double CalculateRealVirtualRemainder(size_t i, const Vec4D &k);
  /*
    YFS: RR_MODE 1 (NLO_Base.C, RealRealRemainder): the double real as the
    exact two-photon remainder beyond REAL_COMBINE's product. CalculateReal()
    keeps each photon's product factor, flux and density.
  */
  struct RRPhotonInfo {
    bool ok = false;
    double delta = 0.;   //!< the single-real bracket beta_1/(crude B), as in the product
    double flux = 1.;    //!< the single real's flux on r
    double denom = 0.;   //!< the single real's density (crude)
    double rho = 0.;     //!< the single real's exact over crude, kappa r flux/(denom B)
    ATOOLS::Vec4D_Vector pt;  //!< the single real's (n+1)-body point, photon last
    double subloc = 0.;  //!< the single real's subtraction over its density
    double factor = 1.;  //!< this photon's factor in REAL_COMBINE's product
  };
  std::vector<RRPhotonInfo> m_rrphot;
  RRPhotonInfo m_lastrrinfo;
  double RealRealRemainder(size_t i, size_t j);
  //! YFS: RR_MODE singles_product (v x RR is then left out of the weight)
  bool RRSinglesProduct() const;
  //! RR_MODE 1, two-photon events: rho_ij, the exact two-photon weight the
  //! product plus RR must equal (read by YFS_Handler's RR_PROBE identity line)
  double m_rrExactPair = 0.;
  //! the event's v = V/B - B(generation legs), as VIRTUAL_COMBINE uses it
  double VBornEvent() const { return m_vborn; }
  //! the pair's (n+2)-body point from the event, the other photons reduced
  bool EventPairPoint(size_t i, size_t j, ATOOLS::Vec4D_Vector &P);
  //! the incoherent crude S~_II + sum S~_FF of photon k on the legs of pt
  double CrudeOnLegs(const ATOOLS::Vec4D_Vector &pt, const Vec4D &k);
  //! v_{n+1} = V_fin/T - B~(point legs) at the (n+1)-body point p (photon
  //! last): the IR-subtracted one-loop over the tree, loop coupling units.
  bool RealVirtualFactor(const ATOOLS::Vec4D_Vector &p, double &v1);
  void ProbeRealVirtual(const Vec4D &kk, const ATOOLS::Vec4D_Vector &p,
                        double r, double rtree, double subloc, double eikloc,
                        double flux, double subb, double tot);
  double CalculateVV();
  /*!
    YFS: VV_POLE_CHECK (default true). With a double-virtual provider that
    reports its massive Laurent series (VV_Generator: Internal for
    e+e- -> gamma gamma, EXTRA_XS/NNLO/ee_aa_QED_VV.H), check per event that
    YFS's own subtraction cancels the two-loop soft poles:
      m1m_-1 = sigma_-1,  m2m_-2 = sigma_-1^2/2,  m2m_-1 = sigma_-1 m1m_0,
    sigma = 2 B/(alpha/pi) of this event's dipoles (BVV_full_eps, shifted to
    mu^2 = s), m1m/m2m the provider's massified one- and two-loop series per
    Born (alpha/2pi units, its own one loop). Diagnostic only: no weight
    changes. YFS: VV_POLE_PROBE n prints the first n events (@@@ VVPOLE).
  */
  void CheckVVPoles();
  void ReportVVPoles();
  //! YFS: BETA0_2_DIPHOTON: Born points with the provider, and those skipped
  //! for BETA0_2_DIPHOTON_POWER_CUT
  long m_b02aa_n = 0, m_b02aa_skip = 0;
  long   m_vvp_n = 0, m_vvp_fail = 0;
  double m_vvp_max[3] = {0., 0., 0.};
  long   m_vvp_digits[3][17] = {{0}};
  //! Assert that the mapped hard system balances the total photon momentum.
  //! Multiplicity-agnostic: pass the summed photon momentum.
  void CheckMappingRecoil(const ATOOLS::Vec4D_Vector &p,
                          const ATOOLS::Vec4D &ksum);
  //! Map the Born configuration onto one accommodating any number of photons.
  //! The photons are boosted in place. The fixed-multiplicity overloads below
  //! are wrappers around this.
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D_Vector &k);
  /*
    The (n+1)-body point for beta_1(k_j) when the event has other photons,
    for photons radiated from the initial state (YFS: REAL_MAP). All three
    keep k_j's direction and the beam directions in the lab, so the eikonal
    of the point is the event's up to a scale, and boost the final state
    rigidly; they differ in how the other photons are taken out. Each returns
    false when it does not apply (a final-state photon among k, unphysical
    kinematics), and MapMomenta then uses the rest-frame construction
    (REAL_MAP: 0, the pre-2026-09-24 behaviour), see the comments in
    NLO_Base.C for the measurements behind each.
      1: beams reduced by the others' longitudinal momentum, photon as is.
      2: beams and photon scaled by one x (KKMC's beta_1 convention). Default.
      3: as 2, but with (final state + k_j)^2 kept as in the event.
  */
  bool MapMomentaBeamAxis(ATOOLS::Vec4D_Vector &p, const ATOOLS::Vec4D_Vector &k);
  bool MapMomentaScaled(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D_Vector &k,
                        bool anylabel = false);
  /*
    The Born process (Comix), for the Born at momenta other than the event's.
    Needed by REAL_SUB_EIK 6: the crude a photon could have been generated
    with carries the Born of EACH assignment, S~_II B(s' with the photon
    removed) + S~_FF B(pre-emission pair), as ACRAIC's rho_crude does.
  */
  PHASIC::Process_Base *p_bornproc = nullptr;
  double m_me2born = -1.;     // Comix |M_0|^2 at m_bornMomenta, cached per event
  double BornME2At(const ATOOLS::Vec4D_Vector &p);
  //! (N+n)!/N! for n extra photons over N photons in the Born; 1 without
  //! Born photons or with YFS: REAL_BORN_PHOTON_SYM: 0.
  double BornPhotonSym(size_t nextra) const;
  /*
    Born-photon multichannel (YFS: REAL_BORN_PHOTON_MULTICHANNEL, NLO_Base.C).
    With photons in the Born final state an initial-state photon that passes
    the generation cuts is itself Born-eligible, and the same final state is
    reached through every assignment of photons to "Born pair" and "ISR"
    that passes them. BornPhotonChannelSum returns
        G = sum_b h_b / h_gen ,   h = prod_{ISR photons} S~_II(k) B(Born pair) s/s'
    over the assignments b that pass the Born selector and the generator's
    v_max (h_gen, the generated assignment, included, so G >= 1 and G = 1
    exactly when no other assignment passes). B is the Born |M|^2 at the
    generator's reduced point of that pair: beams back to back along the lab
    axis in the pair's rest frame (pure boost), the pair as in the event.
      lab   : beams (unused) + the Born final legs in the lab, pre-FSR
      isr   : the ISR photons of the event
      jswap : >= 0 only the assignments that exchange ISR photon jswap with
              one Born photon (the channels of that photon's one-photon real);
              < 0 every choice of Born photons among Born and ISR photons.
    nalt, if given, receives the number of other assignments that passed.
  */
  double BornPhotonChannelSum(const ATOOLS::Vec4D_Vector &lab,
                              const std::vector<ATOOLS::Vec4D> &isr,
                              YFS::Dipole &dII, int jswap,
                              int *nalt = nullptr);
  //! 1/G for the ISR photon k of this event (1 when the switch is off or k
  //! is not an ISR photon of the event); used on the one-photon real.
  double BornPhotonChannelWeight(const ATOOLS::Vec4D &k);
  //! the generator's v_max and s for the channel check, set by YFS_Handler
  double m_bpmc_vmax = 1., m_bpmc_s = -1.;
  //! diagnostics of the last BornPhotonChannelWeight call (WEIGHT_PROBE)
  double m_bpmc_lastG = 1.;  int m_bpmc_lastn = 0;
  double m_bpmc_hardG = 1., m_bpmc_hardR = 0., m_bpmc_hardx = 0.,
         m_bpmc_hardsub = 0.;
  bool   ReducedBornISR(const ATOOLS::Vec4D &k, double &bratio);
  /*
    The Born of the event with photon k moved to the other channel, relative
    to the event's Born: sign = -1 takes k out of the Born beams (the ISR
    alternative of an FSR-labelled photon, invariant (P - k)^2), sign = +1
    puts it back into beams and pair (the FSR alternative of an
    ISR-labelled photon, invariant (P + k)^2). s2ratio = new invariant / P^2.
  */
  bool   ReducedBornAlt(const ATOOLS::Vec4D &k, int sign, double &bratio,
                        double &s2ratio);
  /*
    REAL_SUB_EIK 8. The Born at invariant R^2 relative to the Born of the
    event, both built by the same construction from the event's PRE-FSR
    final state Q_pre (NLO m_plab, the system m_born is evaluated at): beams
    back to back along the lab axis in R's rest frame (the generator's
    reduced-beam convention), final legs along their directions in Q_pre's
    rest frame, rescaled to R^2. With R = Q_pre this is the event's Born
    point (m_me2pre; checked against m_born by SUB8_TRACE).
  */
  double m_me2pre = -1.;
  double m_sub8[4] = {-1., -1., -1., -1.};   // SUB8_TRACE scratch
  bool   PreFSRBornRatio(const ATOOLS::Vec4D &R, double &bratio);
  bool   GeneratorReaches(const ATOOLS::Vec4D_Vector &p);
  bool   PreFSRBornPoint(const ATOOLS::Vec4D &R, ATOOLS::Vec4D_Vector &p) const;
  ATOOLS::Vec4D PreFSRSystem() const;
  //! Final-state photons: Q fixed, the final legs recoil (YFS: REAL_FSR_MAP).
  //! 1: every final leg rescaled in the common rest frame (one resonant pair).
  bool MapMomentaFSR(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D_Vector &k);
  //! 2 (default): only the legs of the photon's OWN dipole recoil, the other
  //! pairs sit at their pre-emission momenta. Reduces to 1 for one pair.
  bool MapMomentaFSRDipole(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D_Vector &k);
  //! the same point for photons the pole scheme's W-pair dipole radiated
  bool MapMomentaFSRPole(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D_Vector &k);
  //! the beams of an FSR-mapped point, at the pre-emission final state's s'
  bool RebuildBeamsAtPreFSR(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D_Vector &k,
                            const char *tag);
  //! The event photon with momentum k (nullptr if none), matched by momentum
  //! as PhotonIsFSR does; carries the dipole it was radiated from.
  const YFS::Photon *FindPhoton(const ATOOLS::Vec4D &k) const;
  //! The event's legs AFTER final-state emission (m_plab holds the legs before).
  ATOOLS::Vec4D_Vector m_postlab;
  /*!
    Per photon of the last CalculateReal(): the ratio of the coherent eikonal
    beta_1 subtracted to the crude it was divided by. YFS_Handler builds the
    IFI_Real weight from THESE, so that the Born term and the subtraction
    telescope to r flux/(S~ B) per photon; with the ratio taken from the
    event dipoles instead, the two coherent eikonals (mapped point vs event
    legs) differ at O(1) for hard final-state photons and the difference
    survived as a heavy tail.
  */
  std::vector<double> m_wifterms;
  inline const std::vector<double> &WIFTerms() const { return m_wifterms; }
  inline void SetPostEmissionMomenta(const ATOOLS::Vec4D_Vector &p) { m_postlab = p; }
  bool MapMomentaInvariant(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D_Vector &k);
  bool PhotonIsFSR(const ATOOLS::Vec4D &k) const;
  // Set by MapMomenta: true when one of the new reductions ran AND the event
  // had other photons to remove, i.e. the mapped point is not the event.
  // CalculateReal then divides beta_1 by the CRUDE eikonal of the mapped
  // point (see there); at n = 1 it is false and the assembly is the legacy
  // one, bit for bit.
  bool m_map_reduced = false;
  /*
    Same construction with the Born reference passed in rather than read from
    m_bornMomenta. ACRAIC needs this: MakeACRAIC runs BEFORE CalculateNLO, so
    NLO's member still holds the previous event and the rotation would be
    taken against the wrong axis.
  */
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k,
                  const ATOOLS::Vec4D_Vector &ref);
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k);
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k1, ATOOLS::Vec4D &k2);
  // Defined in NLO_Checks.C, along with the other *Sub checks and the
  // debug histograms.
  void BookHistograms();
  void WriteHistograms();
  void CheckMassReg();
  //! YFS: VIRTUAL_POLE_PROBE: the loop's and the subtraction's 1/eps per event
  void ProbeVirtualPoles(Define_Dipoles *vsub);
  //! YFS: VIRTUAL_PROBE_ABOVE: large |V/B| re-evaluated in two more frames
  void ProbeLargeVirtual(const double virt, const double sub);
  //! YFS: VIRTUAL_CHECK: the loop in a second frame; unstable -> no virtual
  void CheckVirtualStability(double &virt, const double sub, const bool vcanon);
  //! YFS: RR_LOG_CHECK: the double real's soft x hard log against B~(pt) - B~(gen)
  void ProbeRRLog(const Vec4D &k, const Vec4D_Vector &legs, const double Keps);
  //! YFS: PHOTON_DUMP: a photon of the real loop and its contribution
  void DumpPhoton(const YFS::Photon &g);
  void DumpPhotonContribution(const YFS::Photon &g, const double contrib);
  //! momentum-balance reports of the real's mappings (printed when violated)
  void ReportFSRMapResidual(const Vec4D_Vector &p, const Vec4D_Vector &k,
                            const Vec4D &Kall, const double xi,
                            const double Mp, const double msum);
  void ReportMappingResidual(const Vec4D_Vector &p, const Vec4D_Vector &k,
                             const char *tag, const Vec4D &Q);
  //! YFS: REAL_STAB: the mapped pair invariant of a real point
  void DumpMapQ(const double sqrtsqq, const double eksum, const std::string &extra);
  //! YFS: VIRTUAL_CHECK: end-of-run counts
  void ReportVirtualCheck();
  //! B~(K) and the final-leg form factor on the point's and the generation legs
  void SoftSumsOnBothLegs(const Vec4D_Vector &legs, const double K,
                          double &bt_pt, double &bt_gen,
                          double &dY_pt, double &dY_gen);
  bool CheckPhotonForReal(const Vec4D &k);
  bool CheckPhotonForReal(const Vec4D &k, const Vec4D_Vector &p);
  void CheckMasses(Vec4D_Vector &p, int realmode=0);
  bool CheckMomentumConservation(Vec4D_Vector p);
  void CheckRealSub(Vec4D k, int mode);
  //! Collinear-limit counterpart of CheckRealSub: fixed photon energy, swept
  //! angle, checked against the closed-form eikonal rather than a soft-limit
  //! zero. Triggered by Check_Real_Sub: 3. See NLO_Checks.C for why.
  void CheckRealCollinearSub(Vec4D k, int mode);
  void CheckRealVirtualSub(Vec4D k);
  void CheckRealRealSub(Vec4D k1, Vec4D k2);
  // Accumulating (non-exiting) sub-check scatter: for each real photon record
  // its energy, pT, opening angle to the nearest charged particle, and the
  // subtracted residual, so large residuals for energetic photons can be
  // checked for collinear enhancement. Enabled by m_check_real_sub / m_check_rr_sub == 2.
  void RecordSubScatter(const ATOOLS::Vec4D &k, double residual,
                        const std::string &tag, double eik = 1.0);
  void RandomRotate(Vec4D &k);
  inline void SetBorn(double b) { m_born = b; m_me2born = -1.; m_me2pre = -1.; }
  inline void SetBornProc(PHASIC::Process_Base *p) { p_bornproc = p; }
  inline void SetMomentum(Vec4D_Vector &p) { m_plab = p;}
  inline void SetFSR(YFS::FSR *p) {p_fsr = p;}
  inline bool HasVirtual() {return m_looptool;}
  /*!
    Whether the card asked for the virtual at all (V in NLO_Part), set by
    PHASIC::YFS_Process whether or not a loop provider exists. With
    YFS: ACRAIC_Virtual: acraic and no provider, HasVirtual() is false although
    the virtual was requested, and ACRAIC then added none of its own while
    YFS.NLO took ACRAIC's (virtual-free) factor: no virtual in either column
    (audit 2026-09-23, D3). YFS_Handler tells ACRAIC with this instead.
  */
  bool m_virtrequested = false;
  inline void SetVirtualRequested(bool r) { m_virtrequested = r; }
  inline bool VirtualRequested() const { return m_virtrequested; }
  inline bool HasReal() {return m_realtool;}
  inline bool HasRR() {return m_rrtool;}
  //! An exact double-virtual provider (p_vv); without one YFS_Handler
  //! estimates the VV from the EEX series (VV_EEX, NNLO_VV*).
  inline bool HasVV() {return m_vvtool;}
  //! YFS: EW_HIGHER_ORDER - see NLO_Virtual.C
  void AddEWHigherOrder();
  /*!
    YFS: BETA0_2 / ACRAIC: BETA0_2 - the photonic two-loop beta_0^(2) at
    LL+NLL (YFS/Main/Beta0_2.H for the formulas, NLO_Beta02.C). BuildBeta02
    runs at the end of CalculateVirtual, on the Born point m_plab of the
    loop, and keeps the event's analytic photonic factors; each column
    composes them with its own one-loop delta (ACRAIC in
    YFS_Handler::GenerateWeight, YFS.NLO in YFS_Handler::CalculateNLO).
  */
  void InitBeta02();
  void BuildBeta02();
  void ReportBeta02();
  //! v of the loop provider alone: m_oneloop/m_born less EW_HIGHER_ORDER
  double VOneLoop() const;
  inline bool Beta02Ok() const { return m_b02.ok; }
  inline const YFS::Photonic_Factors &Beta02Factors() const { return m_b02; }
  inline const YFS::Beta02_Config &Beta02Acraic() const { return m_b02acraic; }
  inline const YFS::Beta02_Config &Beta02Nlo()  const { return m_b02nlo; }
  //! ACRAIC with an external virtual: the real delta^(2) of delta = v1/2
  double Beta02External() const;
  //! YFS.NLO, relative to the Born: (1 + v1/2 + d2)^2 - 1 - v1 with
  //! d2 = delta^(2) of delta = v1/2 (the squared-level |1 + delta + d2|^2)
  double Beta02SquaredLevel() const;
  /*!
    YFS: DYSON_VP - the Dyson-resummed photon vacuum polarisation beyond the
    provider's first order (YFS/Main/Dyson_VP.H), both columns. Built with
    beta_0^(2), at the end of CalculateVirtual, at the loop's Born point.
  */
  void BuildDysonVP();
  /*!
    ACRAIC: EW_HIGHER_ORDER_HELICITY - the EW_HIGHER_ORDER term per helicity
    (Griffin_Interface::HigherOrderHelicity), labels and factors delta_G(h)
    in the provider's physical labels, and its Born-weighted 2 Re = v_HO.
  */
  void BuildEWHigherOrderHelicity();
  inline YFS::ewhohel::code EWHOHelicityMode() const { return m_ewhohel; }
  inline double VEWHigherOrder() const { return m_vewho; }
  bool m_ghelok = false;
  std::vector<std::vector<int> > m_ghel;
  std::vector<Complex> m_gdelta;
  inline bool DysonVPOn() const { return p_dysonvp != nullptr; }
  inline const YFS::VP_Factors &DysonVPFactors() const { return m_vp; }
  inline bool HasRealVirtual() {return m_realvirt;}
  inline bool HasNLO()  {return m_realtool || m_looptool;}
  inline bool HasNNLO() {return m_rrtool   || m_realvirt;}
  Vec4D MostEnergeticPhoton() const;
  Vec4D FixedTestPhoton() const;
  // One-shot dump of beta_0 and beta_1 on a single, fully specified phase-space
  // point, for the number-for-number comparison against KKMC's CEEX. Enabled
  // with YFS: ACRAIC_Compare: 1. See Test/SherpaCompare/kkmc_acraic_crosscheck.cxx
  // in the KKMC repo for the other half; the momenta printed here are meant to
  // be fed to it so both sides evaluate at exactly the same point.
  void ACRAICComparePoint();
  bool m_acraic_done;

private:
  /*
    The terms of CalculateReal(k) for one photon, filled stage by stage by
    the helpers below (NLO_Real.C), in the order CalculateReal calls them.
  */
  struct RealTerms {
    ATOOLS::Vec4D kk;           //!< the event's photon
    ATOOLS::Vec4D k;            //!< the photon at the (n+1)-body point
    ATOOLS::Vec4D_Vector p;     //!< the (n+1)-body point, photon last
    ATOOLS::Vec4D_Vector pp;    //!< the point's Born legs
    double r = 0.;              //!< the real, in the Born's normalisation
    double r_nomc = 0.;         //!< r before the Born-photon multichannel share
    double flux = 0.;           //!< the flux r is multiplied by
    double subloc = 0.;         //!< the eikonal beta_1 subtracts
    double subb = 0.;           //!< the event's density (crude)
    double flux0 = 0., subloc0 = 0., subb0 = 0.;   //!< the three before REAL_SUB_EIK
    double subb_loc = 0.;       //!< the density at the point (Sub_Mode global)
    double denom0 = 0.;         //!< subb_loc without REAL_SUB_EIK
    double sub8_bF = -1., sub8_f = -1., sub8_scale = 1.;   //!< SUB8_TRACE only
    double tot;                 //!< beta_1/S~; an unknown Sub_Mode leaves it unset
  };
  void   MapRealPoint(RealTerms &rt);
  void   RealME(RealTerms &rt);
  void   RealFlux(RealTerms &rt);
  void   RealSubtraction(RealTerms &rt);
  void   SubEikPostEmissionFSR(RealTerms &rt);
  void   SubEikCoherentBornFSR(RealTerms &rt);
  void   SubEikAssignmentBornFSR(RealTerms &rt);
  void   SubEikMultichannelFSR(RealTerms &rt);
  void   SubEikPreFSRMultichannelFSR(RealTerms &rt);
  double Sub8FSRChannelDensity(const RealTerms &rt, sub8fsrlegs::code fsrlegs,
                               double sFFgen);
  void   FluxOnFSRCrude(RealTerms &rt);
  void   FillPointHistograms(const RealTerms &rt, const std::string &suffix);
  bool   RealMEUsable(const RealTerms &rt);
  void   PointDenominator(RealTerms &rt);
  void   MultichannelISRPhoton(RealTerms &rt);
  void   PreFSRMultichannelISRPhoton(RealTerms &rt);
  void   Sub8Trace(const RealTerms &rt);
  void   AssembleBeta1(RealTerms &rt);
  void   RecordHardestPhoton(const RealTerms &rt);
  void   RealTrace(const RealTerms &rt);
  void   RealStabilityProbe(const RealTerms &rt);
  void   IFIRealRatio(const RealTerms &rt);
  //! where the last CalculateReal(k) stopped, and whether IFIRealRatio
  //! applied a ratio there; read by CountIFISkip
  enum class real_exit { completed, momentum_conservation, unusable_me };
  real_exit m_real_exit = real_exit::completed;
  bool   m_ifi_applied = false;
  //! IFI_Real bookkeeping for a photon whose IF ratio was not applied
  void   CountIFISkip();
  void   ReportBadReal(const RealTerms &rt);
  void   FillRealHistograms(const RealTerms &rt);
  void   TrackRealAverage(double tot);
  void   StoreRealPointInfo(const RealTerms &rt);
};


}
#endif