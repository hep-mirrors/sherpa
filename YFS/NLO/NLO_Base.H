#ifndef YFS_NLO_BASE_H
#define YFS_NLO_BASE_H

#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/FSR.H"
#include "YFS/NLO/Virtual.H"
#include "YFS/NLO/Real_Correction.H"
#include "YFS/NLO/RealVirtual.H"
#include "YFS/NLO/VirtualVirtual.H"
#include "YFS/Main/YFS_Form_Factor.H"
#include "YFS/Main/Photon.H"
#include "YFS/Main/Define_Dipoles.H"
#include "ATOOLS/Math/Histogram_2D.H"
#include <fstream>
#include <unordered_map>
#include <memory>
#include "ATOOLS/Math/Histogram.H"

#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "PHASIC++/Channels/Multi_Channel.H"




namespace YFS {

typedef std::map<std::string, std::unique_ptr<ATOOLS::Histogram_2D> > histograms2d;
typedef std::map<std::string, std::unique_ptr<ATOOLS::Histogram> > histograms1d;

class NLO_Base : YFS_Base {
public:
  // subloc/subb seen by the IFI_Real restoration, accumulated so the tail that
  // drives the MC error is reportable instead of inferred from sigma.
  long   m_ifi_n = 0;
  double m_ifi_sum = 0., m_ifi_min = 1e300, m_ifi_max = -1e300, m_ifi_sum2 = 0.;
  // prod_i subloc(k_i)/subb(k_i) over one event's photons.
  double m_ifi_prod = 1.;
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
  // Per-nominal-call contributions of the hardest photon(s), captured
  // inline while CalculateReal/CalculateRealVirtual/CalculateRealReal sum
  // over the full photon multiplicity - no separate recomputation needed.
  // *_hard1 = contribution of the single hardest photon; *_hard2 = sum over
  // the two hardest photons (m_rr_hard2 is the single pair formed by them).
  // Used to report the matching truncated to 1 or 2 photons.
  double m_real_hard1, m_rv_hard1, m_rr_hard2;
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
  double m_ceexvirt = 1.;
  void SetProviders(YFS::Virtual *virt, YFS::Real *real,
                    YFS::RealVirtual *realvirt, YFS::RealReal *realreal,
                    YFS::VirtualVirtual *vv);
  /*!
    Whether a loop-ME provider has to be built at all.

    Three things can supply the O(alpha) virtual and only one of them needs an
    external loop ME:
      - an external provider (Loop_Generator: OpenLoops, ...),
      - the EEX collinear virtual (YFS: Collinear_Virtual),
      - CEEX's beta_0^1 - the two vertex form factors and the IFI boxes -
        which is internal and exact, so with CEEX on and no Loop_Generator
        named there is nothing to look up. Without this, NLO_Part: BVR with
        CEEX: 1 and no Loop_Generator died on "Couldn't find virtual ME for
        this process" while holding a perfectly good virtual of its own.

    \param ext true if the card named a real external loop generator. When it
    did, the provider is built even with CEEX on: both are then evaluated and
    CEEX goes out as its own named weight rather than replacing anything.
  */
  inline bool NeedsVirtualProvider(bool ext=true) const
  { return !m_eex_virt && (ext || !m_useceex); }
  /// whether CEEX is the thing supplying the virtual on this run
  inline bool CeexSuppliesVirtual() const
  { return m_useceex && !m_eex_virt && !m_looptool; }
  /*!
    The O(alpha) virtual factor from CEEX, rho(Born+virtual)/rho(Born), pushed
    in by YFS_Handler after MakeCEEX() has run. Only consulted when
    CeexSuppliesVirtual(); it is 1 until CEEX has produced a number, so an
    event where CEEX bailed contributes no virtual rather than a garbage one.
  */
  inline void SetCeexVirtual(double f) { m_ceexvirt = f; }
  /*!
    The two halves of the external virtual, kept so CEEX_Compare can test the
    identity that has to hold if the two are in the same scheme:

        virt/Born - sub  ==  2 Re(deltI + deltF)

    virt is the RAW one-loop from the provider, sub the YFS virtual form
    factor B that is subtracted out of it and put into exp(Y) instead. CEEX's
    deltI/deltF are the same remainder by construction (KKMC's F1ini/F1fin are
    "IR subtracted"), so a mismatch localises to one side of that equation
    rather than to "the virtual disagrees".
  */
  //! NB not m_virt_sub - that name is already the Virtual_Sub SETTING on
  //! YFS_Base, and shadowing it here silently switched the subtraction off.
  double m_virt_raw = 0., m_virt_subval = 0.;
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
  double CalculateVV();
  //! Assert that the mapped hard system balances the total photon momentum.
  //! Multiplicity-agnostic: pass the summed photon momentum.
  void CheckMappingRecoil(const ATOOLS::Vec4D_Vector &p,
                          const ATOOLS::Vec4D &ksum);
  //! Map the Born configuration onto one accommodating any number of photons.
  //! The photons are boosted in place. The fixed-multiplicity overloads below
  //! are wrappers around this.
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D_Vector &k);
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k);
  void MapMomenta(ATOOLS::Vec4D_Vector &p, ATOOLS::Vec4D &k1, ATOOLS::Vec4D &k2);
  // Defined in NLO_Checks.C, along with the other *Sub checks and the
  // debug histograms.
  void BookHistograms();
  void WriteHistograms();
  void CheckMassReg();
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
  inline void SetBorn(double b) { m_born = b;}
  inline void SetMomentum(Vec4D_Vector &p) { m_plab = p;}
  inline void SetFSR(YFS::FSR *p) {p_fsr = p;}
  inline bool HasVirtual() {return m_looptool;}
  inline bool HasReal() {return m_realtool;}
  inline bool HasRR() {return m_rrtool;}
  inline bool HasNLO()  {return m_realtool || m_looptool;}
  inline bool HasNNLO() {return m_rrtool   || m_realvirt;}
  Vec4D MostEnergeticPhoton() const;
  Vec4D FixedTestPhoton() const;
  // One-shot dump of beta_0 and beta_1 on a single, fully specified phase-space
  // point, for the number-for-number comparison against KKMC's CEEX. Enabled
  // with YFS: CEEX_Compare: 1. See Test/SherpaCompare/kkmc_ceex_crosscheck.cxx
  // in the KKMC repo for the other half; the momenta printed here are meant to
  // be fed to it so both sides evaluate at exactly the same point.
  void CEEXComparePoint();
  bool m_ceex_done;
};


}
#endif