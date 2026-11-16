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
#include <functional>
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

/*!
  ACRAIC: LOOP_BETA_MAX 2: one photon pair's two-photon one-loop amplitude per
  helicity, on RVHelicityFactors' remainder convention: rem = delta_S - delta_0,
  delta_S the IR-subtracted (n+2)-point loop over tree, delta_0 the Born's.
*/
struct LoopNHelicityFactors {
  ATOOLS::Vec4D k1, k2;                     //!< the lab photons (provider order)
  std::vector<ATOOLS::Vec4D> ks;            //!< every photon of the set (provider order)
  double v2 = 0.;                           //!< v_{n+2} + 2c (squared level)
  bool ok = false;
  std::vector<std::vector<int> > hel;       //!< provider labels, photons last
  std::vector<Complex> a0, rem;             //!< tree; delta_S - delta_0
};

class NLO_Base : YFS_Base {
public:
  //! This event's s, forwarded to the dipoles and form factor NLO_Base owns.
  void SetS(const double s) override;
  // subloc/subb seen by the IFI_Real restoration (diagnostic sums).
  long   m_ifi_n = 0;
  double m_ifi_sum = 0., m_ifi_min = 1e300, m_ifi_max = -1e300, m_ifi_sum2 = 0.;
  // prod_i subloc(k_i)/subb(k_i) over one event's photons.
  double m_ifi_prod = 1.;
  // B (m_ifi_prod - 1) of the last CalculateReal(): what the IF ratios add to the
  // real sum (0 without IFI_Real); see YFS_Handler.
  double m_real_ifi_soft = 0.;
  // Photons of the real sum above IFIOmega(), and those whose IF factor stayed 1, by reason.
  long m_ifi_eligible = 0, m_ifi_skip_emin = 0, m_ifi_skip_momcons = 0,
       m_ifi_skip_unusable = 0, m_ifi_skip_other = 0;
  // Profiled in x = E_gamma/sqrt(s), 5 bins of 0.1: m_ifi_x_r the reweight applied
  // (1 + S_IF/S_crude, Born momenta), m_ifi_x_e the ratio subloc/subb (mapped).
  long   m_ifi_x_n[5] = {0,0,0,0,0};
  double m_ifi_x_r[5] = {0.,0.,0.,0.,0.};
  double m_ifi_x_e[5] = {0.,0.,0.,0.,0.};


  NLO_Base();
  ~NLO_Base();

  //! Not owned: the providers belong to the YFS_Process they were built for (see SetProviders).
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
  //! Fixed-order real providers by photon multiplicity: [1] is p_real, [2] p_realreal; more only with NLO_MAX_PHOTONS.
  std::vector<YFS::Real_Correction*> m_realprov;
  //! One-loop providers by photon multiplicity: [1] is p_realvirt, [n >= 2] with ACRAIC: LOOP_BETA_MAX.
  std::vector<YFS::RealVirtual*> m_rvprov;
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

  double m_evts, m_recola_evts, m_real,m_ravg, m_rv, m_rvsub, m_eikeex, m_subloc;
  // Contributions of the hardest photon (*_hard1) and of the two hardest (*_hard2),
  // captured while the real, RV and RR sums run, to report the matching truncated
  // to 1 or 2 photons.
  double m_real_hard1, m_rv_hard1, m_rr_hard2;
  /*
    The last CalculateReal(k, raw)'s per-photon conventions: flux, subtraction
    eikonal and denominator (the crude) of tot = (r flux - subloc B/resc)/denominator,
    read by the double real's beta_2 subtraction.
  */
  double m_lastflux = 1., m_lastsubloc = 0., m_lastdenom = 0.;
  //! the same two before any REAL_SUB_EIK / REAL_FSR_FLUX modification
  double m_lastsubloc0 = 0., m_lastdenom0 = 0.;
  double m_real_hard2, m_rv_hard2;
  // Product of the two crude eikonals S~(k1) S~(k2) the last CalculateRealReal(k1,k2) divided by.
  double m_rr_eik;
  int m_zeroV, m_zeroRV, m_zeroRR, m_nonZeroRR, m_nonZeroRV;
  // RV beta_1^1 skipped as unreliable: softer than RV_SOFT_CUT or failing RV_CANCEL_EPS.
  int m_softRV;
  // Double-real pairs skipped by the RR_SOFT_CUT energy pre-filter.
  int m_softRR;
  // RV photons skipped by the RV_ME_MAX_RATIO stability guard (|rv| far above its subtraction's scale).
  int m_rvUnstable;
  // RV_CANCEL_HIST diagnostics: m_rvHiC photons with cancellation C >= 1,
  // m_rvBlowup with |tot| > 1e3 |Born|; m_rvBlowup{Rtree0,HiC,Soft} those of them
  // with rtree == 0, C >= 1 or E/sqrt(s) < 0.01, m_rvBlowupHardWide with
  // E/sqrt(s) > 0.1 and no charged leg within ~26 deg.
  int m_rvHiC, m_rvBlowup, m_rvBlowupRtree0, m_rvBlowupHiC, m_rvBlowupSoft,
      m_rvBlowupHardWide;
  // CalculateReal(k,raw) calls where the raw amplitude p_real->Calc_R(p) was zero.
  int m_zero_real_amp;
  bool m_failcut;


  void Init(Flavour_Vector &flav, Vec4D_Vector &plab, Vec4D_Vector &born);
  /*!
   * Point this process-agnostic handler at the providers of the currently
   * selected process (call alongside YFS_Handler::SetFlavours). Null means no
   * such correction; ownership is never taken.
   */
  double m_acraicvirt = 1.;
  void SetProviders(YFS::Virtual *virt, YFS::Real *real,
                    YFS::RealVirtual *realvirt, YFS::RealReal *realreal,
                    YFS::VirtualVirtual *vv);
  /*!
    Whether a loop-ME provider has to be built. The O(alpha) virtual can come
    from an external provider (Loop_Generator), the EEX collinear virtual
    (YFS: Collinear_Virtual) or ACRAIC's internal beta_0^1, which needs no
    lookup.

    \param ext true if the card named an external loop generator; the provider
    is then built even with ACRAIC on, and ACRAIC goes out as its own weight.
  */
  inline bool NeedsVirtualProvider(bool ext=true) const
  { return !m_eex_virt && (ext || !m_useacraic ||
                           // ACRAIC_Virtual external names the provider even for Loop_Generator Internal
                           m_acraicvirtsrc == acraicvirt::external); }
  /*!
    Whether ACRAIC's own virtual is used rather than the provider's, as chosen
    by YFS: ACRAIC_Virtual; `acraic` also needs no provider to be present.
  */
  inline bool AcraicSuppliesVirtual() const
  { return m_useacraic && !m_eex_virt && !m_looptool
           && (m_acraicvirtsrc == acraicvirt::acraic ||
               // auto: ACRAIC's own when V is requested and no external provider exists
               (m_acraicvirtsrc == acraicvirt::automatic && m_virtrequested)); }
  /*!
    The O(alpha) virtual factor rho(Born+virtual)/rho(Born) from ACRAIC, pushed
    in by YFS_Handler after MakeACRAIC(). Used only when AcraicSuppliesVirtual();
    1 until ACRAIC has produced one.
  */
  inline void SetAcraicVirtual(double f) { m_acraicvirt = f; }
  /*!
    The two halves of the external virtual, for ACRAIC_Compare's scheme check
    virt/Born - sub == 2 Re(deltI + deltF): virt is the raw one-loop, sub the
    YFS virtual form factor B subtracted from it.
  */
  //! Not m_virt_sub, which is the Virtual_Sub setting of YFS_Base.
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
    YFS: ACRAIC_Virtual: helicity. The IR-finite complex virtual factor per
    helicity of the loop provider at the Born point m_plab:

      delta^h = A1^h/A0^h - run/2 - alpha B_4

    A0, A1 the provider's tree and one-loop amplitudes (Delta_IR = 0), run the
    uniform vacuum-polarisation correction Virtual::Calc_V subtracts, Re 2 alpha
    B_4 the YFS virtual subtraction, Im alpha B_4 fixed so the common imaginary
    part is KKMC's (LoopHelicityIR). The Born-weighted average of 2 Re delta^h
    is v = m_oneloop/m_born. Labels are the provider's physical helicities
    (Acraic_Base::LoopHelicityToAcraic).
  */
  bool m_lhelok = false;
  std::vector<std::vector<int> > m_lhel;
  std::vector<Complex> m_la0, m_ldelta;
  ATOOLS::Vec4D_Vector m_lhelmom;
  //! the raw per-helicity A1^h/A0^h of the Born point (m_lratio_ok false where the
  //! tree is below VIRT_HEL_MINBORN), and the running-alpha correction and Im part
  //! delta^h was built with; used by BuildRealVirtualHelicityFactors
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
    Z = theta Q (theta = +1 incoming, -1 outgoing), s_ij = (p_i +- p_j)^2, so
    d(A1^h/A0^h)/d ln mu_IR^2 = K. imsh is the lambda-free imaginary part
      imsh = sum_{i<j} Im K_ij ln(mu_IR^2 e^c / |s_ij|),
    c = ACRAIC: VIRT_HEL_IM_C. False if a charged leg is massless.
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

  // ---- general-n real corrections: beta_n by the YFS subset recursion
  //
  //   beta_n(k_1..k_n) = M_n * flux
  //                    - sum over non-empty S subset of {1..n} of
  //                      [ prod_{i in S} Stilde(k_i) ] * beta_{n-|S|}({k}\S)
  //
  // with beta_0 = Born.

  //! The real-correction provider for a photon multiplicity, or NULL if the
  //! run was not set up with one. Index 0 is unused.
  YFS::Real_Correction *RealProvider(size_t nphotons) const;
  //! Register the provider for `nphotons`, growing the table as needed.
  void SetRealProvider(size_t nphotons, YFS::Real_Correction *prov);
  //! One-loop provider for `nphotons` real photons (1: p_realvirt), or NULL.
  YFS::RealVirtual *RealVirtualProvider(size_t nphotons) const;
  void SetRealVirtualProvider(size_t nphotons, YFS::RealVirtual *prov);
  //! ACRAIC: LOOP_BETA_MAX when ACRAIC is on at ORDER 2, else 1: the highest
  //! multiplicity to BUILD a one-loop provider for.
  static size_t RequestedMaxLoopPhotons();
  //! one record per eligible photon pair (filled by CalculateRealVirtual)
  std::vector<LoopNHelicityFactors> m_loopnhel;
  //! the pair (i, j) of m_photons; false (and !out.ok) if it cannot be built
  bool BuildLoopNHelicityFactors(size_t i, size_t j, LoopNHelicityFactors &out);
  //! the same for any photon set S of m_photons (|S| = 2: the pair record;
  //! |S| = 3: ACRAIC: LOOP_LADDER_PROBE's three-photon record)
  bool BuildLoopNSetFactors(const std::vector<size_t> &S, LoopNHelicityFactors &out);
  void BuildLoopNPairs();
  long m_lnh_n = 0, m_lnh_fail = 0, m_lnh_polefail = 0, m_lnh_nfall = 0;
  long m_lnh_cone = 0, m_lnh_checked = 0, m_lnh_unstable = 0;
  double m_lnh_res_sum = 0., m_lnh_res_max = 0.;
  //! Highest photon multiplicity a provider exists for (0 if none).
  size_t MaxRealPhotons() const;
  //! Runcard ceiling NLO_MAX_PHOTONS on the fixed-order photon multiplicity (default 2).
  static size_t RequestedMaxRealPhotons();
  //! Highest multiplicity to BUILD a provider for: NLO_MAX_PHOTONS, or ACRAIC:
  //! TREE_BETA_MAX at ORDER 2. YFS.NLO's R_n still stops at RequestedMaxRealPhotons().
  static size_t RequestedMaxProviderPhotons();

  /*!
    Matrix element times flux for the photons selected by `mask`, evaluated in
    the mapping that multiplicity induces, together with the LOCAL eikonals
    Stilde in that same frame (written into `subloc`, indexed as `photons`).
    Returns false if the point is unusable.
  */
  bool RealMEForSubset(unsigned mask, double &me, std::vector<double> &subloc);

  /*!
    beta_|mask| by the recursion above, memoised over subsets so that each
    distinct subset is evaluated once.
  */
  double RawBeta(unsigned mask);

  /*!
    Raw (unnormalised) residuals computed this event, keyed by a bitmask over
    m_photons, so beta_1(k) and beta_2(i,j) are shared by every superset. Separate
    from Photon::m_beta10 and PhotonPair::m_beta20, which hold the final
    normalised contributions.
  */
  std::unordered_map<unsigned, double> m_rawbeta;

  //! beta_n / (crude eikonal product) for the photons selected by `mask` out of
  //! m_photons - the weight the caller adds, for any n.
  double CalculateRealN(unsigned mask);

  //! Sum beta_n over every n-photon subset of this event's photons. n=2
  //! reproduces the pair loop in CalculateRealReal(); n=3 needs no new loop,
  //! only a three-photon ME provider.
  double CalculateRealMultiplicity(size_t n);
  double CalculateRealVirtual();
  double CalculateRealVirtual(Vec4D k);
  /*
    YFS: RV_MODE 1 (CalculateRealVirtualRemainder): the real-virtual as the
    non-factorisable remainder on the real's own (n+1)-body point, density, flux
    and coupling. CalculateReal(k) leaves the point and the photon's exact-over-
    crude ratio in m_lastrvinfo; CalculateReal() stores one per photon in m_rvinfo.
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
  // ACRAIC: TNP_IFI: v_IF = (v(p) - v(pbar))/2, pbar = Born point with the final momenta swapped
  bool m_vifok = false;
  double m_vifodd = 0.;
  long m_vifn = 0, m_viffail = 0;
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
  long m_rvBeamCone = 0, m_rvIsrSkipped = 0, m_rvPoleVeto = 0;
  bool RVInBeamCone(const ATOOLS::Vec4D_Vector &p) const;
  //! ACRAIC: HYBRID_ISR_RMIN: does this initial-state photon get V_2?
  bool InHybridIsrCone(const ATOOLS::Vec4D &k) const;
  //! YFS: VIRTUAL_CHECK - Born-level loops checked in a second frame / dropped
  long m_vLoopChecked = 0, m_vLoopUnstable = 0;
  double m_vLoopUnstableBorn = 0., m_vLoopBornSum = 0.;
  //! YFS: VIRTUAL_RMIN counters
  long m_vLoopRmin = 0;
  double m_vLoopRminBorn = 0.;
  //! YFS: VIRTUAL_FRACTION - Born-level loops evaluated / replaced by the
  //! reference v0, and the sum of the evaluated v (the running v0)
  long m_vFracEvaluated = 0, m_vFracSkipped = 0;
  //! ACRAIC: RV_FRACTION - the exact real-virtual loops on a fraction of the
  //! events (REAL_VIRTUAL averaged); m_rvSkipped marks an event without them
  long m_rvFracEvaluated = 0, m_rvFracSkipped = 0;
  bool m_rvSkipped = false;
  double m_vFracSum = 0.;
  //! the reference v0 of YFS: VIRTUAL_FRACTION for the next event
  double VirtualFractionV0() const;
public:
  //! a uniform number in [0, 1) from the bits of the momenta (VIRTUAL_FRACTION,
  //! RV_FRACTION): the selection follows the kinematics, not the random stream
  static double KinematicsHash(const ATOOLS::Vec4D_Vector &p);
  //! ACRAIC: RV_FRACTION in (0, 1], default 1
  static double RealVirtualFraction();
  inline bool RealVirtualSkipped() const { return m_rvSkipped; }
  //! v with the subtraction on the Born point's own legs (RV_MODE 1), and
  //! that subtraction itself (loop coupling units)
  bool m_vborn_own_ok = false;
  double m_vborn_own = 0., m_vborn_own_sub = 0.;
  bool BornVirtualOnOwnLegs(double &vb);
  //! VIRTUAL_SUB_LEGS own_if_exponent: Y_IF(own legs) - Y_IF(event legs)
  double IFExponentLegShift();
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
  //! the beams' invariant mass squared at each m_rvdv entry's (n+1)-point: below s
  //! when other initial-state photons were taken out of the beams
  std::vector<double> m_rvdv_spt;
  //! the last RealVirtualFactor's loop over tree, before its subtraction,
  //! that subtraction (B~/kappa on the point's legs), and the point in the
  //! loop frame
  double m_rvLastLt = 0., m_rvLastSub = 0.;
  ATOOLS::Vec4D_Vector m_rvLastLoop;
  /*!
    ACRAIC: REAL_VIRTUAL helicity. The (n+1)-point provider's per-helicity
    amplitudes, read after photon i's RealVirtualFactor call, as the remainder
    of RVHelicityFactors:
        delta_{n+1}^c = A1^c/A0^c - (1/2) B~(point legs)/kappa
                        - i imsh(point) + c_photon/2 ,
        delta_B^c'    = A1^c'/A0^c' (Born point) - run/2
                        - (1/2) B~(Born's own legs)/kappa - i imsh(Born) ,
        rem^c = delta_{n+1}^c - delta_B^c' ,
    c' the Born configuration with the same fermion labels (photon last).
    Configurations below VIRT_HEL_MINBORN keep rem = dv/2 (Born: delta_B = v_B/2).
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
    YFS: RR_MODE 1 (RealRealRemainder): the double real as the exact two-photon
    remainder beyond REAL_COMBINE's product, from the per-photon product factor,
    flux and density CalculateReal() keeps.
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
  //! the pair's (n+2)-body point from the event, the other photons reduced
  bool EventPairPoint(size_t i, size_t j, ATOOLS::Vec4D_Vector &P);
  //! EventPairPoint for any photon set: Born legs, then the photons of S in S order
  bool EventSetPoint(const std::vector<size_t> &S, ATOOLS::Vec4D_Vector &P);
  //! the incoherent crude S~_II + sum S~_FF of photon k on the legs of pt
  double CrudeOnLegs(const ATOOLS::Vec4D_Vector &pt, const Vec4D &k);
  //! v_{n+1} = V_fin/T - B~(point legs) at the (n+1)-body point p (photon
  //! last): the IR-subtracted one-loop over the tree, loop coupling units.
  bool RealVirtualFactor(const ATOOLS::Vec4D_Vector &p, double &v1);
  void ProbeRealVirtual(const Vec4D &kk, const ATOOLS::Vec4D_Vector &p,
                        double r, double rtree, double subloc, double eikloc,
                        double flux, double subb, double tot);
  /*!
    YFS: VV_POLE_CHECK (default true). With a double-virtual provider that
    reports its massive Laurent series (VV_Generator: Internal, e.g.
    EXTRA_XS/NNLO/ee_aa_QED_VV.H), check per event that YFS's subtraction
    cancels the two-loop soft poles:
      m1m_-1 = sigma_-1,  m2m_-2 = sigma_-1^2/2,  m2m_-1 = sigma_-1 m1m_0,
    sigma = 2 B/(alpha/pi) of the event's dipoles, m1m/m2m the provider's
    massified series (alpha/2pi units). Diagnostic only.
  */
  void CheckVVPoles();
  void ReportVVPoles();
  //! YFS: BETA0_2_DIPHOTON: Born points with the provider, and those skipped
  //! for BETA0_2_DIPHOTON_POWER_CUT
  long m_b02aa_n = 0, m_b02aa_skip = 0;
  long   m_vvp_n = 0, m_vvp_fail = 0;
  double m_vvp_max[3] = {0., 0., 0.};
  long   m_vvp_digits[3][17] = {{0}};
  //! Map the Born configuration onto one accommodating any number of photons,
  //! boosting the photons in place; the fixed-multiplicity overloads wrap this.
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D_Vector &k);
  /*
    The (n+1)-body point for beta_1(k_j) when the event has other photons, for
    initial-state photons (YFS: REAL_MAP). All keep k_j's direction and the
    beam directions in the lab and boost the final state rigidly; they differ in
    how the other photons are removed. Each returns false when it does not
    apply, and MapMomenta then uses the rest-frame construction (REAL_MAP: 0).
      1: beams reduced by the others' longitudinal momentum, photon as is.
      2: beams and photon scaled by one x (KKMC's beta_1 convention). Default.
      3: as 2, but with (final state + k_j)^2 kept as in the event.
  */
  bool MapMomentaBeamAxis(ATOOLS::Vec4D_Vector &p, const ATOOLS::Vec4D_Vector &k);
  bool MapMomentaScaled(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D_Vector &k,
                        bool anylabel = false);
  /*
    The Born process (Comix), for the Born at momenta other than the event's.
    Needed by REAL_SUB_EIK 6, whose crude carries the Born of each photon
    assignment, as ACRAIC's rho_crude does.
  */
  PHASIC::Process_Base *p_bornproc = nullptr;
  double m_me2born = -1.;     // Comix |M_0|^2 at m_bornMomenta, cached per event
  double BornME2At(const ATOOLS::Vec4D_Vector &p);
  //! (N+n)!/N! for n extra photons over N photons in the Born; 1 without
  //! Born photons or with YFS: REAL_BORN_PHOTON_SYM: 0.
  double BornPhotonSym(size_t nextra) const;
  /*
    Born-photon multichannel (YFS: REAL_BORN_PHOTON_MULTICHANNEL). With photons
    in the Born final state, an initial-state photon passing the generation cuts
    is itself Born-eligible, so one final state arises from every assignment of
    photons to "Born pair" and "ISR". BornPhotonChannelSum returns
        G = sum_b h_b / h_gen ,   h = prod_{ISR photons} S~_II(k) B(Born pair) s/s'
    over the assignments b passing the Born selector and the generator's v_max
    (so G >= 1); B is the Born |M|^2 at the pair's reduced point.
      lab   : the Born final legs in the lab, pre-FSR (beams unused)
      isr   : the ISR photons of the event
      jswap : >= 0 only the assignments exchanging ISR photon jswap with one
              Born photon; < 0 every choice of Born photons
    nalt, if given, receives the number of other assignments that passed.
    With a non-default YFS: BORN_PHOTON_CHANNELS stage, isr holds every YFS
    photon and h carries S~_II + sum_FF S~_FF per photon (NLO_Real.C).
    density, if given, replaces h: the generator's density of the assignment
    (slot, pick), slot the Born photon legs and pick[a] the candidate on
    slot[a] (Acraic_Base::CrudeOfPhotonAssignment).
  */
  typedef std::function<double(const std::vector<size_t> &,
                               const std::vector<size_t> &)> ChannelDensity;
  double BornPhotonChannelSum(const ATOOLS::Vec4D_Vector &lab,
                              const std::vector<ATOOLS::Vec4D> &isr,
                              YFS::Dipole &dII, int jswap,
                              int *nalt = nullptr,
                              const ChannelDensity *density = nullptr);
  //! 1/G for the ISR photon k of this event (1 when the switch is off or k
  //! is not an ISR photon of the event); used on the one-photon real.
  double BornPhotonChannelWeight(const ATOOLS::Vec4D &k);
  //! the generator's v_max and s for the channel check, set by YFS_Handler
  double m_bpmc_vmax = 1., m_bpmc_s = -1.;
  //! diagnostics of the last BornPhotonChannelWeight call (WEIGHT_PROBE)
  int m_bpmc_lastn = 0;
  bool   ReducedBornISR(const ATOOLS::Vec4D &k, double &bratio);
  /*
    The Born of the event with photon k moved to the other channel, relative to
    the event's Born: sign = -1 takes k out of the Born beams (invariant
    (P - k)^2), sign = +1 puts it back into beams and pair ((P + k)^2).
    s2ratio = new invariant / P^2.
  */
  bool   ReducedBornAlt(const ATOOLS::Vec4D &k, int sign, double &bratio,
                        double &s2ratio);
  /*
    REAL_SUB_EIK 8. The Born at invariant R^2 relative to the event's Born, both
    built from the event's pre-FSR final state Q_pre (m_plab): beams back to
    back along the lab axis in R's rest frame, final legs along their
    directions in Q_pre's rest frame, rescaled to R^2. R = Q_pre gives the
    event's Born point (m_me2pre).
  */
  double m_me2pre = -1.;
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
    IFI_Real weight from these, so that the Born term and the subtraction
    telescope to r flux/(S~ B) per photon.
  */
  std::vector<double> m_wifterms;
  inline const std::vector<double> &WIFTerms() const { return m_wifterms; }
  inline void SetPostEmissionMomenta(const ATOOLS::Vec4D_Vector &p) { m_postlab = p; }
  bool MapMomentaInvariant(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D_Vector &k);
  bool PhotonIsFSR(const ATOOLS::Vec4D &k) const;
  // Set by MapMomenta: true when a new reduction ran and the event had other
  // photons, so the mapped point is not the event; CalculateReal then divides
  // beta_1 by the mapped point's crude eikonal. False at n = 1.
  bool m_map_reduced = false;
  /*
    As above with the Born reference passed in rather than read from
    m_bornMomenta (MakeACRAIC runs before CalculateNLO, so the member still
    holds the previous event).
  */
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k,
                  const ATOOLS::Vec4D_Vector &ref);
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k);
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k1, ATOOLS::Vec4D &k2);
  // Defined in NLO_Checks.C, along with the other *Sub checks and the
  // debug histograms.
  void WriteHistograms();
  void CheckMassReg();
  //! YFS: VIRTUAL_CHECK: the loop in a second frame; unstable -> no virtual
  void CheckVirtualStability(double &virt, const double sub, const bool vcanon);
  //! is a Born photon within r < rmin of a beam?
  bool BornPhotonNearBeam(const double rmin) const;
  void ProbeLargeVirtual(const double virt, const double sub);
  //! YFS: VIRTUAL_CHECK: end-of-run counts
  void ReportVirtualCheck();
  //! B~(K) and the final-leg form factor on the point's and the generation legs
  void SoftSumsOnBothLegs(const Vec4D_Vector &legs, const double K,
                          double &bt_pt, double &bt_gen,
                          double &dY_pt, double &dY_gen);
  void CheckMasses(Vec4D_Vector &p, int realmode=0);
  bool CheckMomentumConservation(Vec4D_Vector p);
  void CheckRealSub(Vec4D k, int mode);
  //! Collinear-limit counterpart of CheckRealSub: fixed photon energy, swept
  //! angle, checked against the closed-form eikonal (Check_Real_Sub: 3).
  void CheckRealCollinearSub(Vec4D k, int mode);
  void CheckRealVirtualSub(Vec4D k);
  void CheckRealRealSub(Vec4D k1, Vec4D k2);
  // Accumulating sub-check scatter: each real photon's energy, pT, angle to the
  // nearest charged particle and subtracted residual (m_check_real_sub / m_check_rr_sub == 2).
  void RecordSubScatter(const ATOOLS::Vec4D &k, double residual,
                        const std::string &tag, double eik = 1.0);
  inline void SetBorn(double b) { m_born = b; m_me2born = -1.; m_me2pre = -1.; }
  inline void SetBornProc(PHASIC::Process_Base *p) { p_bornproc = p; }
  inline void SetMomentum(Vec4D_Vector &p) { m_plab = p;}
  inline void SetFSR(YFS::FSR *p) {p_fsr = p;}
  inline bool HasVirtual() {return m_looptool;}
  /*!
    Whether the card asked for the virtual at all (V in NLO_Part), set by
    PHASIC::YFS_Process whether or not a loop provider exists. YFS_Handler tells
    ACRAIC with this, since HasVirtual() is false with ACRAIC_Virtual: acraic
    and no provider.
  */
  bool m_virtrequested = false;
  inline void SetVirtualRequested(bool r) { m_virtrequested = r; }
  inline bool VirtualRequested() const { return m_virtrequested; }
  inline bool HasReal() {return m_realtool;}
  //! An exact double-virtual provider (p_vv); without one YFS_Handler
  //! estimates the VV from the EEX series (VV_EEX, NNLO_VV*).
  inline bool HasVV() {return m_vvtool;}
  //! YFS: EW_HIGHER_ORDER - see NLO_Virtual.C
  void AddEWHigherOrder();
  /*!
    YFS: BETA0_2 / ACRAIC: BETA0_2 - the photonic two-loop beta_0^(2) at LL+NLL
    (YFS/Main/Beta0_2.H, NLO_Beta02.C). BuildBeta02 runs at the end of
    CalculateVirtual on the loop's Born point m_plab; each column composes the
    event's factors with its own one-loop delta (YFS_Handler::GenerateWeight,
    YFS_Handler::CalculateNLO).
  */
  void InitBeta02();
  void BuildBeta02();
  void ReportBeta02();
  //! v of the loop provider alone: m_oneloop/m_born less EW_HIGHER_ORDER
  double VOneLoop() const;
  //! ACRAIC: TNP_IFI: v_IF of this event (0 when not built), see m_vifodd
  inline double VirtualIFOdd() const { return m_vifok ? m_vifodd : 0.; }
  inline bool VirtualIFOddOk() const { return m_vifok; }
  //! the one-loop V_sub of the Born point with the two final momenta exchanged
  bool MirroredOneLoop(bool vcanon, YFS::Define_Dipoles *vsub, bool ownlegs, double &onem);
  void ReportVirtualIFOdd() const;
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
  // One-shot dump of beta_0 and beta_1 on a fully specified phase-space point for
  // comparison with KKMC's CEEX (YFS: ACRAIC_Compare: 1); the momenta printed are
  // fed to Test/SherpaCompare/kkmc_acraic_crosscheck.cxx in the KKMC repo.
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
  bool   RealMEUsable(const RealTerms &rt);
  void   PointDenominator(RealTerms &rt);
  void   MultichannelISRPhoton(RealTerms &rt);
  void   PreFSRMultichannelISRPhoton(RealTerms &rt);
  void   AssembleBeta1(RealTerms &rt);
  void   IFIRealRatio(const RealTerms &rt);
  //! where the last CalculateReal(k) stopped, and whether IFIRealRatio
  //! applied a ratio there; read by CountIFISkip
  enum class real_exit { completed, momentum_conservation, unusable_me };
  real_exit m_real_exit = real_exit::completed;
  bool   m_ifi_applied = false;
  //! IFI_Real bookkeeping for a photon whose IF ratio was not applied
  void   CountIFISkip();
  void   ReportBadReal(const RealTerms &rt);
  void   TrackRealAverage(double tot);
  void   StoreRealPointInfo(const RealTerms &rt);
};


}
#endif