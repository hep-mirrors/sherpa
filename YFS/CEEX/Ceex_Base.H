#ifndef YFS_CEEX_Ceex_handler
#define YFS_CEEX_Ceex_handler

#include "ATOOLS/Phys/Flavour.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "EXTAMP/Process.H"
#include "YFS/Main/YFS_Base.H"


using namespace MODEL;
using namespace ATOOLS;
using namespace std;
using namespace METOOLS;

namespace METOOLS { class Spin_Amplitudes; }

namespace YFS {

  class NLO_Base;   // MapMomenta only; not owned


  class Real_Correction;

  typedef std::map<int,double> HelicityMap;
  typedef complex<double> dcmplx;

  /*!
    A helicity amplitude over the external legs.

    Storage is FLAT and indexed through Idx(), so the leg count is a runtime
    property rather than four hard-coded subscripts. That is what lets the same
    container hold a 2 -> 2 amplitude (4 legs, 16 entries) and a 2 -> N one;
    the previous Complex m_A[2][2][2][2] could not, and no amount of work
    elsewhere could get past it.

    Idx packs leg 0 into the least significant bit, matching the convention
    METOOLS::Spin_Amplitudes uses and therefore the one Comix hands over.

    s_nlegs is set ONCE per run from the process (SetLegs), and the constructor
    zeroes only the active 2^s_nlegs entries - not the whole capacity, which
    would memset 4 KB per partition and is why the capacity can afford to be
    generous.
  */
  class Amplitude{
  public:
    static const int s_maxlegs = 8;              //!< 2 -> 4 plus a photon is 7
    static const int s_maxhel  = 1 << s_maxlegs;
    static int s_nlegs;            //!< legs carrying TWO helicity states
    static inline int NHel() { return 1 << s_nlegs; }
    static void SetLegs(int n);
    /*!
      Spin-aware version: counts only the legs with more than one spin state
      and records where each leg sits in the packed index.

      METOOLS indexes a Spin_Amplitudes as a MIXED-RADIX product over the
      per-leg spin multiplicities (Spin_Structure::GetNumber:
      num += mult*spins[i]; mult *= m_spins[i]), and a scalar has
      multiplicity ONE - it contributes index 0 and does not advance mult. So
      for e+e- -> H mu+mu- the 2x2x1x2x2 mixed-radix index collapses to
      h0 + 2h1 + 4h3 + 8h4, which is exactly this binary packing over the four
      two-valued legs. The two conventions already agree; what was wrong was
      counting the scalar as a leg, which made CEEX ask for 2^5 entries
      against Comix's 16.

      s_bit[i] is leg i's bit position, or -1 for a leg with one spin state.
    */
    static std::vector<int> s_bit;
    static void SetLegs(const ATOOLS::Flavour_Vector &flavs);
    static inline int Bit(size_t i)
    { return i < s_bit.size() ? s_bit[i] : -1; }
    static inline int Idx(int h0, int h1, int h2, int h3)
    { return h0 | (h1 << 1) | (h2 << 2) | (h3 << 3); }
    Amplitude();
    ~Amplitude(){}
    Complex m_A[s_maxhel];
    Complex m_U[2][2];
    Complex m_V[2][2];
  };

  class Ceex_Base : public YFS::YFS_Base {  
  public:

    Ceex_Base(const Flavour_Vector &flavs);
    virtual ~Ceex_Base() { ReportBornShape(); ReportComixReal(); ReportBeta2(); ReportLoopHelicity(); }

    Flavour_Vector m_flavs;
    /*!
      Sherpa processes supplying the spin amplitudes, when CEEX_COMIX_AMPS is
      on: the Born and the one-photon real. Comix's FillAmplitudes() returns
      helicity-indexed amplitudes that reproduce the matrix element exactly -
      measured, Differential/SumSquare = 1/(4 n!) to the last bit over 5864
      points - so they can replace the hand-coded spinor algebra below.
      Not owned.
    */
    PHASIC::Process_Base *p_bornproc = nullptr;
    PHASIC::Process_Base *p_realproc = nullptr;
    inline void SetBornProc(PHASIC::Process_Base *p) { p_bornproc = p; }
    inline void SetRealProc(PHASIC::Process_Base *p) { p_realproc = p; }
    //! Born spin amplitudes from Comix; false if unavailable.
    bool ComixBornAmplitude(const ATOOLS::Vec4D_Vector &p, Amplitude &A,
                            double *me2 = NULL,
                            const double propscale = -1.0,
                            const double decscale = -1.0);
    /*
      The COMIX leg bitmask of the outgoing fermion pair - the DECAY line. At
      2 -> 2 this is the whole final state, so it names the same current the
      s-channel branch already pins and the second scale is a no-op.
    */
    //! Per-event quantities the weight probe reports; see CEEX: WEIGHT_PROBE.
    inline size_t NPhot()  const { return m_allphotons.size(); }
    inline double SvarQ()  const { return m_svarQ; }
    inline double SvarY()  const { return m_svarY; }
    //! Legs and photons of the current event (CEEX frame), for the
    //! heavy-event trace; see CEEX: TRACE_FACTOR_ABOVE.
    void PrintKinematics(std::ostream &o) const;
    inline void SetLabMomenta(const ATOOLS::Vec4D_Vector &p) { m_plabmom = p; }
    //! the event's final legs BEFORE final-state emission (YFS m_reallab),
    //! where the generator's m_born sits; see GeneratorBornAt
    inline void SetPreFSRMomenta(const ATOOLS::Vec4D_Vector &p) { m_prefsr = p; }
    inline size_t DecayCId() const
    { return (((size_t)1)<<m_if1) | (((size_t)1)<<m_if2); }

    /*!
      CEEX: COMIX_REAL. When set, the O(alpha) REAL correction of a
      one-photon event is taken from Comix's helicity amplitudes for the
      e+e- -> f fbar gamma process instead of the hand-coded spinor algebra
      in Ceex_Real.C. Default 0 - nothing changes unless it is asked for.

      Only the real is affected: the virtual stays where it is, and events
      with a photon multiplicity other than one keep the hand-coded path
      (see ApplyComixReal for why).
    */
    int m_comixreal = 0;
    /*!
      Real matrix-element providers, indexed by photon multiplicity. One per
      multiplicity because the n-photon amplitude comes from a DIFFERENT
      process (this one plus n photons), not from a recursion over the
      one-photon one. Slot 0 is unused.
    */
    /*!
      Whether a virtual was ASKED FOR - i.e. whether NLO_Part contains V.
      YFS_Base::m_looptool is the LOOP_TOOL setting and says nothing about it;
      NLO_Base shadows that name with (p_virt != NULL), which is the real
      answer, so the handler pushes it in.
    */
    bool m_hasvirtual = false;
    inline void SetHasVirtual(bool h) { m_hasvirtual = h; }
    /*!
      Does CEEX build its OWN virtual into its amplitude?

      Only when it was asked for one AND it is the nominated source. Without
      the first condition CEEX added a virtual to a Born+Real run: measured on
      e+e- -> H mu+ mu- with NLO_Part: BR, the NLO column had no virtual at all
      (YFS.Virtual = 0) while the CEEX column still carried its own, and the
      two differed by 36% for that reason alone. Without the second, CEEX would
      add its own on top of a provider's.
    */
    inline bool CeexOwnVirtual() const
    { return m_hasvirtual && m_ceexvirtsrc == ceexvirt::ceex; }

    std::vector<YFS::Real_Correction*> m_realprov;
    /*
      Not owned. Only for MapMomenta: beta_1(k_j) belongs at the REDUCED
      kinematics with the OTHER photons removed, and that construction - boost
      to the Q frame, rebuild the beams - already exists there. Re-deriving it
      here would repeat a mapping that has carried a sign bug before.
    */
    inline void SetNLO(YFS::NLO_Base *p) { p_nlo = p; }

    inline void SetRealProvider(size_t n, YFS::Real_Correction *p)
    { if (m_realprov.size() <= n) m_realprov.resize(n+1, NULL);
      m_realprov[n] = p; }
    inline YFS::Real_Correction *RealProvider(size_t n) const
    { return n < m_realprov.size() ? m_realprov[n] : NULL; }
    Amplitude m_comixM1;
    Amplitude m_comixraw[2];
    int m_comixdrawnhel = 0;
    /*!
      CEEX: COMIX_REAL_FLIP - which external legs Comix labels with the
      OPPOSITE helicity index to CEEX. Bit i flips leg i: 0,1 the beams, 2,3
      the outgoing pair, 4 the photon. See the note in FetchComixReal; the
      @@@ CEEXFLIP diagnostic is what determines it.
    */
    int m_comixflip = -1;
    /*!
      CEEX: COMIX_REAL_PHOTON_FLIP - the one bit the Born cannot determine.
      The Born has no photon, so the calibration below fixes the four FERMION
      bits and is blind to the fifth. The soft limit is blind to it too, and
      for a reason: the soft factor is a scalar in helicity space with
      |s_+| = |s_-|, so the wrong photon helicity is an overall constant there.
      It comes from the hard-photon comparison against the hand-coded real,
      which is decisive (0.046 with the bit set against 0.84 without, on 563
      of 570 events).
    */
    int m_comixphoflip = 1;
    /*!
      CEEX: COMIX_REAL_NORM - Comix's amplitude in CEEX's normalisation.

      MEASURED, not derived: the soft probe (SoftProbe) finds
      |M1_comix|^2 / |e^2 (S_ini+S_fin) Born_hand|^2 -> 4.0000000 as the
      photon softens, on every point tried and to seven digits, so the two
      conventions differ by exactly a factor 2 in the amplitude. The factor
      is invisible in the hand-coded CEEX because its weight is
      rho1/RhoCrude and the same soft factors sit in both; it is NOT
      invisible once the numerator comes from somewhere else.
    */
    double m_comixnorm = -1.;
    //! Has the Born calibration been attempted yet? One attempt per run.
    bool m_comixcalibrated = false;
    /*!
      What a Born-level calibration returns. The Born is where the
      Comix-to-CEEX map can be DERIVED rather than fitted: it has no photon,
      so the leg permutation and the overall normalisation are the only two
      things that can differ, and both are visible directly - no soft limit,
      no reference eikonal, no cancellation to lean on.
    */
    struct ComixCalib {
      int    mask = -1;     //!< winning fermion flip mask
      double met = -1.;     //!< its shape metric
      double next = -1.;    //!< the runner-up's, i.e. how decisive the win is
      double N = 0.;        //!< sqrt(sum|A_comix|^2 / sum|e^2 A_hand|^2)
      double sh = 0., sc = 0., me2 = 0.;
      double rmin = 0., rmax = 0.;   //!< per-helicity |C/H| spread
      int    nlive = 0;     //!< helicity slots carrying real weight
    };
    /*!
      Per-photon, per-helicity factor turning the hand-coded one-photon
      amplitude into Comix's. Built once per event per photon and applied
      inside the partition loop, where the soft weights sProd/Sactu live.
      A RATIO rather than a difference, because it must survive being applied
      to each stage's share of the photon separately: r is a constant in the
      partition index, so sum_stages w_s B_s r = r sum_stages w_s B_s.
    */
    std::vector<Amplitude> m_cxratio;
    //! Comix's beta_1 per photon, kept so it can be compared with what the
    //! partition loop actually accumulates (m_realphot), rather than with a
    //! standalone reconstruction of it.
    std::vector<Amplitude> m_cxcsub;
    std::vector<char>      m_cxratiook;
    int m_perphoton = 0;
    YFS::NLO_Base *p_nlo = nullptr;
    //! Per-PHOTON accounting, kept apart from the per-EVENT counters above:
    //! one run has as many of these as it has photons, not as it has events,
    //! and averaging the two together produces a number that is neither.
    long   m_cxppn = 0, m_cxppfail = 0;
    double m_cxppdev = 0.;
    void BuildComixPhotonRatios();
    //! |beta_1|/beta_0 against E_gamma, the real-subtraction validation
    //! figure, for the hand-coded and the Comix beta_1 on one normalisation.
    void Beta1Scan();
    /*!
      Per-helicity factor turning Comix's Born into CEEX's, built once per
      event. Its modulus is COMIX_REAL_NORM; its phase is the per-leg spinor
      convention difference, which the Born calibration showed is NOT constant
      across helicities and so cannot be one number. Scale independent if the
      two codes' couplings agree - which is exactly what using it at every
      partition scale tests.
    */
    //! Per-helicity spread of the virtual factor over one event's partitions.
    std::vector<double> m_vpmin, m_vpmax;
    bool m_vpon = false;
    Amplitude m_cxbalign;
    bool      m_cxbalignok = false;
    int       m_comixborn = 0;
    /*
      The normalisation between Comix's amplitudes and CEEX's is a true
      CONSTANT only when the two differ by the spin-average alone, which is
      the 2 -> 2 case (1/sqrt(initial spin states), exact). Past that the
      hand-coded Born carries the process's own factors - for H l+ l- the
      Higgs coupling and propagator - and the ratio moves with the event.
      Calibrating it ONCE then makes the whole run depend on which event
      happened to calibrate it: measured over three seeds, the derived
      normalisation came out 13.587, 6.869 and 7.273, and the CEEX column with
      COMIX_REAL followed it to -8.7%, +446% and +14.2%.

      So m_normexact says whether the constant is the exact one. When it is
      not, the real takes m_cxrnorm instead - the same ratio rebuilt per event
      from that event's Born, which is legitimate because whatever the
      hand-coded Born is missing multiplies the one-photon real identically
      and cancels in the ratio.
    */
    bool      m_normexact = false;
    /*
      The DECAY-side invariant of the current partition,
      (q_if1 + q_if2 + sum_FSR k)^2. Built in the partition loop as PY, the
      same object the pseudo-flux numerator uses. At 2 -> 2 it equals svarX.
    */
    double    m_svarY = -1.;
    ATOOLS::Vec4D m_PXvec;   //!< this partition's X, the momentum at s'
    /*!
      The generator's own phase-space point at s': the ISR-reduced beams from
      YFS_Handler's m_plab together with the recoiled final legs. Not read by
      the weight any more - beta_0 is built per partition at X_wp by
      BornLegsAt, and this single per-event point is the all-ISR case of that
      - but kept as the handler already supplies it.
    */
    ATOOLS::Vec4D_Vector m_plabmom;
    ATOOLS::Vec4D_Vector m_prefsr;
    /*!
      The generator's Born point at invariant R^2: the PRE-FSR final legs,
      boosted to their own rest frame and rescaled together to invariant mass
      sqrt(R^2), with the beams back to back along the lab axis as a pure
      boost carries it. The same construction as NLO_Base::PreFSRBornPoint
      (which CEEX cannot call: MakeCEEX runs before NLO's momenta are
      refreshed). For any final state, where BornLegsAt rebuilds one pair.
    */
    bool GeneratorBornAt(const ATOOLS::Vec4D &R, ATOOLS::Vec4D_Vector &pb) const;
    /*!
      Project the final legs and every photon onto exact momentum
      conservation against the beams and onto their mass shells, before
      anything is evaluated on them. See RepairMomentumBalance.
    */
    bool RepairMomentumBalance();
    //! the conservation residual found before the repair (FRAME_PROBE)
    double m_repairresid = -1.;
    /*!
      True when the final state is a fermion-antifermion pair of one flavour
      (mu+ mu-, nu nubar, e+ e-, q qbar): the only case the hand-coded 2 -> 2
      spinor structure, its Comix alignment and its virtual describe. Every
      other process, four legs or not (gamma gamma included), takes the
      general Comix branch.
    */
    bool m_ffbar = false;
    //! beta_1 partition spread, see CEEX: BETA1_PARTITION_CHECK
    double    m_b1min = 0., m_b1max = 0., m_b1sum = 0., m_b1sq = 0.;  int m_b1n = 0;
    double    m_cxrnorm = -1.;
    void BuildComixBornAlignment();
    //! Normalisation for Comix REAL amplitudes: exact constant, else per event
    inline double RealNorm() const
    { return (m_normexact || !(m_cxrnorm > 0.)) ? m_comixnorm : m_cxrnorm; }
    /*
      e(alpha_YFS)/e(alpha_model) for each photon of a Comix real amplitude.
      Comix couples every photon with the model's alpha (1/131.9 in the G_mu
      scheme); CEEX's soft factors carry m_e = sqrt(4 pi m_alpha), which is
      alpha(0) under YFS: USE_MODEL_ALPHA 0. Without this factor the real
      exceeds eikonal x Born by alpha_model/alpha(0) in the soft limit: Z-pole
      mu mu, G_mu, CEEX/YFS.NLO = 1.031. YFS.NLO applies the same ratio
      (m_rescale_alpha) to its real. Exactly 1 with USE_MODEL_ALPHA 1.
    */
    inline double ComixPhotonCoupling(size_t nphot = 1) const
    { return std::pow(m_rescale_alpha, 0.5*(double)nphot); }
    bool ComixRealAt(const ATOOLS::Vec4D_Vector &pp, int hel, Amplitude &M1,
                     const double propscale = -1.0);
    /*
      Comix's IR-subtracted one-photon amplitude for photon k, at an imposed
      propagator scale. beta_1 = M1(k) - lambda M1(lambda k): under k -> lam k
      the polarisation vector depends only on the photon DIRECTION and every
      p.k picks up one lambda, so the eikonal obeys s(lam k) = s(k)/lam
      exactly and the soft piece cancels between the two terms. No
      polarisation vector, no relative phase, no projection - both terms are
      Comix's own amplitude at CEEX's own (non-conserving) arguments.
    */
    /*
      The beta-provider form of the weight, built WITHOUT the partition sum:
        A(hel) = beta_0 + sum_i beta_1(k_i)/s_i + ...
        rho    = sum_hel |A|^2
      Run alongside the partition sum so the two can be compared - the
      partition sum is the IFI reference until the beta form is shown to
      reproduce it. Returns false if it cannot handle this multiplicity.
    */
    void Beta1Closure();
    void SoftLimitTest();
    void SoftNormCheck(const ATOOLS::Vec4D &k);
    bool ReducedBeams(const ATOOLS::Vec4D &k,
                      ATOOLS::Vec4D_Vector &pb) const;
    bool BornLegsAt(const ATOOLS::Vec4D &X,
                    ATOOLS::Vec4D_Vector &pb) const;
    //! Beams back-to-back at X, the radiating pair at Y, spectators as they are
    bool LegsAt(const ATOOLS::Vec4D &X, const ATOOLS::Vec4D &Y,
                ATOOLS::Vec4D_Vector &pb) const;
    /*!
      The legs for photon iphot in the CURRENT partition: the other photons
      on a reducing stage come out of the beams, the others on the final
      stage go into the radiating pair, so that {legs + k} balances and the
      two s-channel invariants M_1 sees are X_wp with this photon counted
      initial and X_wp with it counted final.
    */
    bool PartitionLegs(int iphot, ATOOLS::Vec4D_Vector &pb) const;
    Complex TotalEikonal(const ATOOLS::Vec4D_Vector &p,
                         const ATOOLS::Vec4D &k, int hel);
    /*!
      The polarisation vector Comix's external current carries for an
      OUTGOING photon k at the helicity ComixRealAt selects for CEEX's hel:
      the conjugate of METOOLS' EP/EM built from ATOOLS spinors with the
      gauge vector Amplitude::SetGauge(0) installs, (1,0,1,0). The eikonal
      subtracted from a Comix amplitude has to be built from THIS vector, or
      the two differ by an azimuth-dependent phase that no |s|^2 test sees.
    */
    ATOOLS::Vec4C ComixPolarisation(const ATOOLS::Vec4D &k, int hel) const;
    //! The GPS (KKMC) form of the per-leg eikonal, kept for comparison only.
    Complex SfactorLegGPS(const ATOOLS::Vec4D &p, const ATOOLS::Vec4D &k,
                          int hel);
    //! one stage's eikonal current for photon k, on the legs p (m_pceex order)
    Complex StageEikonal(int stage, const ATOOLS::Vec4D_Vector &p,
                         const ATOOLS::Vec4D &k, int hel);
    bool ComixInfraredSubtracted_1_0(const ATOOLS::Vec4D &k, int hel,
                                    int iphot);
    /*!
      Comix at the PHYSICAL legs with the s-channel propagator momenta shifted
      per side (COMIX::Amplitude::SetPropShift): dini added to initial-side
      lines, dfin to final-side ones. This is KKMC's construction - spinors at
      the physical momenta, propagators at the partition's X - and it needs no
      momentum rebuild, so the eikonal M_1 reduces to is the physical one.
    */
    typedef PHASIC::Process_Base::Prop_Shifts PropShifts;
    bool ComixBornShifted(const ATOOLS::Vec4D_Vector &p, Amplitude &A,
                          const PropShifts &shifts);
    bool ComixRealShifted(const ATOOLS::Vec4D_Vector &pp, int hel,
                          Amplitude &M1, const PropShifts &shifts);
    /*!
      One (leg mask, photon sum) per stage for the CURRENT partition, leaving
      photon iphot out (all photons for iphot < 0). Handed to Comix so that
      every internal line sits at its partition invariant.
    */
    PropShifts StageShifts(int iphot) const;
    /*!
      The partition's reduced legs as LEG entries of the shift list, from
      which COMIX::Amplitude::SetPropShifts places every space-like exchange
      line (one initial leg with part of the final state) at the partition's
      reduced invariant. iphot < 0: the Born, BornLegsAt(m_PXvec);
      iphot >= 0: M_1 of that photon, PartitionLegs. Detected by leg content,
      no process is named. CEEX: TCHANNEL_SHIFT = 0 disables it.
    */
    void AddExchangeLineShifts(int iphot, PropShifts &sh) const;
    /*!
      Does the Born have a space-like exchange line? Probed once per Born
      process: a leg-only shift (which COMIX::Amplitude::SetPropShifts
      applies to exchange lines and to nothing else) either changes the
      Born or leaves it bit for bit. Decides CEEX: TCHANNEL_REDUCED_BORN -1.
    */
    bool BornHasExchangeLine();
    std::map<const PHASIC::Process_Base*, int> m_exchline;
    //! this event: reduced-point beta_0 and balanced-point M_1 (see above)
    bool m_redborn = false;
    /*!
      The final-state radiating groups, as leg indices into m_flavs, set by
      the handler from YFS's radiating dipoles (one per resonance pair).
      Empty means one flat final stage. The initial state is always its own
      stage, appended last (m_initstage) so the two-stage case keeps the old
      ordering.
    */
    std::vector<std::vector<int> > m_groups;
    int m_initstage = 1;
    inline void SetStageGroups(const std::vector<std::vector<int> > &g) { m_groups = g; }
    //! Sum of the photons other than iphot (all if iphot < 0) on the reducing
    //! (initial) stages, or on the others, in THIS partition.
    ATOOLS::Vec4D PartitionShift(int iphot, bool reducing) const;
    bool ComixBeta1At(const ATOOLS::Vec4D &k, int hel,
                      const double propscale, Amplitude &B1);
    bool ComixOnePhotonAmplitude(const ATOOLS::Vec4D &k, int hel,
                                 ATOOLS::Vec4D_Vector &bp, Amplitude &M1,
                                 ATOOLS::Vec4D &kmap);
    bool HandOnePhotonAmplitude(const ATOOLS::Vec4D_Vector &bp,
                                const ATOOLS::Vec4D &k, int hel, Amplitude &M1);
    bool CalibrateComixMap(ComixCalib &c);
    void DeriveComixMap();
    bool m_havecomixreal = false;
    //! Fetch and re-index Comix's one-photon amplitude for this event.
    bool FetchComixReal();
    //! Replace the hand-coded real with it. No-op unless m_comixreal.
    void ApplyComixReal();
    //! Soft-limit ladder: beta_1 must vanish as the photon softens.
    void SoftProbe();
    //! End-of-run summary of the Comix real path.
    void ReportComixReal() const;
    //! Run totals of the Comix-vs-hand comparison (SHERPA_CEEX_COMIX).
    double m_cxrnormsum = 0., m_cxrmetsum = 0.;
    long   m_cxrn = 0, m_cxrfail = 0;
    //! Events taken by the Comix path, and refused by it, per photon
    //! multiplicity - so that "it runs" is not mistaken for "it was used".
    std::vector<long> m_cxrnbyn, m_cxrfailbyn;
    //! First reason the Comix path declined at each multiplicity.
    std::vector<std::string> m_cxrwhy;
    //! Agreement with the hand-coded amplitude, per photon multiplicity.
    //! Aggregated over all n it is meaningless: at n = 1 the hand-coded side
    //! is the exact one-photon matching, at n > 1 it is the partition sum,
    //! and those are different objects being compared to different things.
    std::vector<double> m_cxrnormbyn, m_cxrmetbyn;
    void CountComixReal(size_t ng, bool ok);
    Amplitude  m_bornAmp, m_Amp1U, m_Amp1V;
    Amplitude  m_AmpExpo0, m_AmpExpo1;
    Amplitude  m_snapBorn, m_snapVirt, m_snapReal;
    /*!
      THIS partition's Born WITHOUT its soft-factor product: e^2 (flux)
      beta_0(X_wp), the object InfraredSubtractedME_0_0 multiplied by
      prod_i s_i^{wp_i} and added. beta_1 subtracts it once per photon, with
      that photon's own stage eikonal taken on the legs M_1 was evaluated on
      and the other photons' from the partition. Stored rather than
      recomputed so that the addition and the subtractions are the same
      numbers.
    */
    Amplitude  m_partborn0;
    //! this partition's pseudo-flux svarY/svarQ (1 under NO_PSEUDOFLUX)
    double     m_pflux = 1.;
    //! print the beta_1 pieces per partition on the CHECK_XS dump event
    bool       m_b1trace = false;
    Amplitude  m_AmpBornVirt;
    Amplitude  m_AmpBornReal;

    std::vector<Amplitude> m_realphot;
    Vec4D_Vector m_allphotons_lab;
    double     m_resultbv = 0., m_resultbr = 0.;
    double     m_rhocrud = 0.;
    Vec4D_Vector   m_momenta, m_isrphotons, m_bornmomenta;
    Vec4D_Vector   m_fsrphotons, m_allphotons;
    int                  m_nstages = 2;
    std::vector<int>     m_stage;
    std::vector<char>    m_stagereduces;
    /*!
      Photons below SOFT_PARTITION_CUT are not enumerated: their sum over
      stages collapses to the TOTAL eikonal times beta_0 up to a propagator
      shift of order k. They keep a fixed stage (the one with the largest
      eikonal, so the X shift follows their natural emitter) and the
      soft-factor product uses s_tot(k) for them (m_sactu). That takes the
      cost from nstages^n to nstages^(n_hard).

      A collapsed photon is NOT automatically below BETA1_XCUT: the cut here
      is x = 2E/sqrt(s), beta_1's is y = 2k.S/S^2 against the system the
      photon is radiated from, and after a hard ISR photon S^2 << s. Such a
      photon still gets a beta_1 - M_1 added once per partition (it is never
      enumerated, so every assignment of the others is visited exactly once)
      against s_tot(k) beta_0(X_fixed), the same object m_sactu put into
      beta_0. See ComixInfraredSubtracted_1_0.
    */
    std::vector<char>    m_fixedstage;
    std::vector<Complex> m_sactu;
    /*!
      The crude's soft-factor product for this partition: |s^wp_j|^2 for an
      enumerated photon and sum_g |s_g(k_j)|^2 for a collapsed one. The
      generator labels every photon initial or final and produces it with
      |S_I|^2 or |S_F|^2, so its label-blind density carries the INCOHERENT
      sum per photon - never |s_I + s_F|^2, which (i) drops the soft-photon
      initial-final interference from rho_1/rho_crude and (ii) vanishes at
      the zeros of s_I + s_F while beta_1 does not. m_crudefixed says whether
      any photon is collapsed (the product is then not |m_Sprod|^2).
    */
    double m_crudeprod = 1.;
    bool   m_crudefixed = false;
    std::vector<std::vector<Complex> > m_Sfac;
    /*!
      One entry per leg of a stage: the index into m_pceex/m_bornmomenta and
      the weight w = Q * theta (charge in units of e, theta = -1 incoming,
      +1 outgoing). A stage's soft factor is sum_i w_i * (per-leg eikonal), and
      Ceex_Base::Sfactor is exactly that sum for the charge-neutral two-leg
      case, with implicit weights (+1,-1). Built by BuildStages().
    */
    struct StageLeg { int leg; double w; };
    std::vector<std::vector<StageLeg> > m_stagelegs;
    /*!
      Momenta the stage legs index. The first m_flavs.size() entries are the
      external legs - initial from the Born momenta, final from the CEEX ones,
      which is the split the ISR/FSR version made implicitly - and anything
      appended after them is a RECONSTRUCTED RESONANCE.

      That last part is why stage legs index this rather than m_pceex directly:
      a W in e+e- -> W+W- -> 4f radiates and belongs to three stages at once
      (production, and its own decay), but is not an entry in the event's
      flavour vector at all. DipoleSet::BuildPole already reconstructs it for
      the pole scheme; a stage leg has to be able to name it.
    */
    std::vector<ATOOLS::Vec4D> m_stagemom;
    //! (q1+q2)^2, the PHYSICAL final-pair invariant. Fixed for the event,
    //! where m_sp is the propagator scale and changes with the partition.
    double m_svarQ = 0.;
    int    m_nparts = 0;   //!< partitions summed on the last event
    size_t m_maxnphot = 0; //!< highest photon multiplicity seen
    long   m_partbad = 0, m_partn = 0;
    double m_pzmin = 0., m_pzmax = 0., m_pzhardx = 1e30;
    double m_pzspread[16] = {0.}, m_pzsum[16] = {0.};
    long   m_pzn[16] = {0};
    const double m_pzxcut[4] = {1e-5, 1e-4, 1e-3, 1e-2};
    double m_pzsoft[4] = {0.}, m_pzsoftsum[4] = {0.}, m_pzsoftsav[4] = {0.};
    long   m_pzsoftn[4] = {0};
    Vec4D m_eta, m_zeta, m_b, m_pp1, m_pp2, m_kk;
    /*!
      Boost from the lab into the frame Init() puts m_momenta in. Needed
      because Init() boosts m_momenta ONLY, and is called BEFORE
      SetISRPhotons()/SetBornMomenta() - so the photons and the Born momenta
      arrive in the lab and would otherwise be combined with CMS momenta in
      the same spinor product.
    */
    ATOOLS::Poincare m_cms;

    /*!
      The four momenta the CEEX spin amplitudes are built from:
      {p_a, p_b, p_c, p_d} with p_c,p_d the PHYSICAL outgoing fermions.

      m_momenta is {beam, beam, Born pair, lab pair} - MakeCEEX appends the
      lab final state after the four Born momenta - and only the lab pair
      satisfies p_a + p_b = p_c + p_d + sum(k). The Born pair sits at full
      beam energy, so using it puts the final-state spinors at a different
      kinematic point from the propagator, which eq.(one-photon) evaluates at
      X = P - sum(k).
    */
    ATOOLS::Vec4D_Vector m_pceex;
    void BuildCeexMomenta();

    /*!
      Self-validating check on the Born spin amplitude, accumulated over the
      run and reported at destruction.

      For e+e- -> f fbar through gamma/Z the Born must go like
      a(1 + (b/a) cos(theta) + cos^2(theta)): the cos^2 coefficient is exactly
      1 and the linear one carries the forward-backward asymmetry. Fitting
      that needs no external normalisation.


      m_born IS therefore a valid comparator, and fixes the normalisation
      absolutely. Compared at matched kinematics - the Born rebuilt at the
      physical pair's scale, not at a partition's m_sp - the measured relation
      is

          e^4 sum_hel |A_born|^2 / 4  =  m_born * (s/4) / s'

      i.e. m_born = e^4 sum_hel |A_born|^2 * s'/s. The 1/s' is the flux factor
      m_lastxs carries and the helicity sum does not. Verified at two beam
      energies, which is what distinguishes s/4 from a constant: 15625.1
      measured against s/4 = 15625 at 250 GeV, and 2079.06 against 2079.36 at
      91.2 GeV.

    */
    double m_fit[3][3] = {{0.}}, m_fitb[3] = {0.}, m_fity = 0., m_fity2 = 0.;
    long   m_fitn = 0;
    //! Running totals of the O(alpha^0) and O(alpha) CEEX distributions, so
    //! the size of the virtual correction can be reported at end of run.
    double m_rho0sum = 0., m_rho1sum = 0.;
    //! Same running totals split by which beta enters, so the virtual and
    //! the real can be read off separately. In YFS the big collinear logs in
    //! beta_0^1 must be cancelled by beta_1^0; if they are not, these two
    //! will not have opposite signs.
    double m_rhobvsum = 0., m_rhobrsum = 0.;
    //! Incoherent (RhoCrud) vs coherent (rho0) Born. The CEEX event weight
    //! is r1/RhoCrud; these two coincide only for a SINGLE partition, so
    //! their ratio measures how much the denominator choice matters once
    //! FSR opens the partition sum.
    double m_rhocrudsum = 0.;  long m_rhon = 0;
    bool   m_ceexdumped = false, m_justdumped = false;
    //! Worst relative violation of the CEEX soft-factor identity, over the run.
    double m_sfacworst = 0.;
    long   m_sfacn = 0;
    //! Check |s_sigma(k)|^2 against the closed form of hep-ph/0006359
    //! eq.(soft-fac-isr). Self-contained: no second generator needed.
    void   CheckSoftFactor(const ATOOLS::Vec4D &k);
    //! Check U and V against the DIAGONALITY identity of hep-ph/0006359
    //! eq.(diagonality): with both fermion legs at the same momentum p,
    //! U[p,k_sigma,p] = V[p,k_sigma,p] = b_sigma(k,p) delta_{l1 l2}.
    void   CheckUVDiagonality(const ATOOLS::Vec4D &k);
    double m_uvoff = 0., m_uvdiag = 0., m_bsig = 0., m_uvflip = 0.;
    long   m_uvn = 0, m_uvflipn = 0;
    void   AccumulateBornShape(double cth, double val);
    void   ReportBornShape();
    int m_onlyz, m_onlyg;
    bool m_fixedwidth,m_checkxs;
    std::vector<int> m_PhoHel;
    PHASIC::Process_Info m_pi;
    Complex    m_T, m_U, m_sW, m_cW, m_sin2tw, m_gcpl, m_propG, m_propZ, m_cfac;
    Complex m_cL, m_cR, m_prop, m_I_L, m_I_R, m_F_L, m_F_R, m_I, m_Soft;
    /*!
      Overall convention factor of the eikonal against Comix's M_1. Measured
      on the soft limit with ComixPolarisation in place: M_1 -> -(s C)
      exactly, the phase settling on -pi to 1e-6 as lambda -> 0 on every
      photon and configuration tried, so the vertex-sign convention differs
      from Comix's by a sign and nothing else.
    */
    Complex m_sfacphase = Complex(-1., 0.);
    Complex m_Tamp[2][2][2][2], m_Uamp[2][2][2][2],m_ampborn[2][2][2][2], m_Sc, m_Tc, m_Uc;
    Complex m_Sini[2][2], m_TC[2], m_UC[2], m_ampISR1,  m_ae, m_af, m_ve, m_vf, m_norm;
    //! Bhabha t-channel: gamma/Z exchanged between the (p1,p3) electron line
    //! and the (p2,p4) positron line. Same two spinor structures as the
    //! s-channel (Fierz), so only these coupling x propagator factors are new.
    Complex m_propGt, m_propZt, m_TCt[2], m_UCt[2];
    //! Electroweak form factors, KKMC's GPS_EWFFact structure but computed
    //! rather than interpolated from DIZET tables. Tree values are
    //! kappa = rho = VPi = 1, which reproduces KKMC with KeyElw = 0 - i.e.
    //! exactly what this code did before CEEX: WEAK existed.
    Complex m_kapE{1.,0.}, m_kapF{1.,0.}, m_kapEF{1.,0.};
    Complex m_rhoEW{1.,0.}, m_gamVPi{1.,0.}, m_vvcor{1.,0.};
    //! 2*T3 and 4*Q*sin^2(theta_W) kept separately so the kappa factors can
    //! multiply the sin^2 piece alone, as in GPS_EWFFact.
    Complex m_t3e2, m_t3f2, m_qesw, m_qfsw;
    int     m_weak = 0;
    //! t-channel spinor structures: the 2<->3 crossing of the s-channel
    //! pairings, i.e. (2,1)x(3,4) and (2,3)x(1,4).
    Complex m_Tampt[2][2][2][2], m_Uampt[2][2][2][2];
    double  m_tinv = 0.;
    bool    m_bhabha = false;
    Complex m_BoxGGtu, m_BoxGZtu, m_BoxGGut, m_BoxGZut, m_Sprod, m_beta00, m_beta10, m_beta01, m_beta20;
    // Complex m_Amp1U[2][2], m_Amp1V[2][2];
    Complex m_vertexI, m_vertexF, m_b1, m_b2;
    Complex m_rawboxGG, m_rawboxGZ, m_rawboxSub, m_rawcoef;
    double     m_Q1Q2I, m_Q1Q2F, m_QIQF, m_MZ, m_gZ, m_isrformfactor, m_sQ;
    double     m_mass_I, m_mass_F, m_qe, m_qf, m_crude;
    /*!
      Where the outgoing fermion pair sits. At 2 -> 2 this is legs 2 and 3 and
      the whole coupling layer said so literally. With anything else in the
      final state it is not: for e+e- -> H l+ l- leg 2 is the HIGGS, so
      m_qf = flavs[2].Charge() made the final-state Z coupling that of a
      neutral scalar and CouplingZ reported "Z coupling is Zero!".
    */
    size_t m_if1 = 2, m_if2 = 3;
    double     m_result, m_result0, m_e;
    // Re sum_h A_0(h)^* A_1(h) / 4: the Born-level coherent amplitude against
    // the full O(alpha^1) one. What an external virtual multiplies (see
    // YFS_Handler::GenerateWeight, CEEX_Virtual: external).
    double     m_result01 = 0.;
    void RegisterDefaults();
    void Init(const Vec4D_Vector &p);
    virtual void Calculate();
    //! Zero the cumulative spin amplitudes (KKMC ZerAmplit).
    void ZerAmplit();
    //! sum_hel |A(hel)|^2 for each order (KKMC MakeRho).
    void MakeRho();
    void InfraredSubtractedME_0_0();
    void InfraredSubtractedME_0_1();
    /*!
      beta_1^0 for ONE photon attached to the INITIAL state, KKMC's HiniPlus.
      sProd is the product of the S-factors of every photon in the current
      PARTITION and Sactu is this photon's own, so the factor KKMC calls
      Fact0 = sProd/Sactu is passed in rather than rebuilt - with a partition
      sum the caller is the only place that knows either.
    */
    /*!
      beta_1^0 for ONE photon attached to the FINAL state, KKMC's HfinPlus.

      CKine = (q1+q2+k)^2 / (q1+q2)^2. beta_0^0 carries svarX/svarQ (KKMC's
      Cfact0), which is 1 for a purely initial-state partition but not once a
      photon sits on the final line; the sProd*(1-CKine)*Born term here is
      what compensates it, and it vanishes in the infrared as CKine -> 1.
    */
    /*!
      Derive the emission stages from the process and check they are admissible.

      A stage is gauge invariant only if it is separately charge-neutral,
      sum_{i in g} Q_i theta_i = 0. That is the condition this enforces; a decomposition failing it is rejected
      rather than silently producing a gauge-dependent answer.

      Currently returns the 2 -> 2 decomposition (initial pair, final pair),
      which is the KKMC ISR/FSR split. This is the ONE function a 2 -> N core
      has to change: everything downstream reads m_stagelegs and m_nstages.
    */
    //! Helicity index into Amplitude::m_A. Forwarded so the convention lives
    //! in one place (Amplitude) while call sites inside Ceex_Base stay short.
    static inline int Idx(int h0, int h1, int h2, int h3)
    { return Amplitude::Idx(h0, h1, h2, h3); }
    Complex SfactorLeg(const ATOOLS::Vec4D &p, const ATOOLS::Vec4D &k, int hel);
    bool BuildStages();
    void PartitionStart(int &last);
    void PartitionPlus(int &last);
    /*!
      O(alpha^2) tree-level double real, beta_2 (CEEX: ORDER: 2), built from
      Comix's two-photon amplitude; see Ceex_Beta2.C. The beta_1 pass leaves,
      per partition and photon, the aligned M_1 it added (if this partition
      carries it), the stage eikonal it subtracted with, and whether beta_1
      was active at all, so that beta_2 subtracts EXACTLY the beta_1 the
      partition sum contains (the "(1-n) cancellation must use the same
      numbers" rule, one order up).
    */
    int m_order = 1;
    Amplitude m_AmpBeta2;     //!< the beta_2 increments alone (partition summed)
    Amplitude m_AmpExpo2;     //!< m_AmpExpo1 + m_AmpBeta2, built in MakeRho
    /*!
      CEEX: KKMC_FLUX_EMULATION (diagnostic, 2 -> 2, no virtual): KKMC's
      O(alpha^1) amplitude keeps sum_{i<j in F} 2 k_i.k_j/Q^2 sProd B per
      partition (its beta_0 flux minus the (1-CKine) terms), and its O(alpha^2)
      removes it only for pairs above vcut2. These hold the two residuals so
      rho_1/rho_2 can be made KKMC-like for an event-level comparison.
    */
    Amplitude m_AmpFluxAll, m_AmpFluxSoft;
    std::vector<Amplitude> m_b2M1;    //!< aligned M_1(k_j) of THIS partition
    std::vector<char>      m_b2hasM1; //!< this partition carried M_1(k_j)
    std::vector<char>      m_b2on;    //!< beta_1(k_j) active in this partition
    std::vector<Complex>   m_b2sj;    //!< the eikonal beta_1(k_j) subtracted with
    std::vector<char>      m_b2elig;  //!< photon may pair (enumerated, x cut)
    bool m_b2ready = false;           //!< provider, process and Born checked
    long m_b2pairs = 0, m_b2m2 = 0, m_b2fail = 0, m_b2events = 0;
    //! per event: which photons may pair; false = no beta_2 this event
    bool PrepareBeta2();
    void ReportBeta2() const;
    //! zero the per-partition beta_1 record (called at each partition start)
    void ResetBeta1Record();
    //! Comix M_2(k_j,k_l) in CEEX's alignment, photon j in the first photon slot
    bool ComixReal2Shifted(const ATOOLS::Vec4D_Vector &pp, int hj, int hl,
                           Amplitude &M2, const PropShifts &shifts,
                           bool noshift = false);
    //! StageShifts with both photons j and l left out
    PropShifts StageShifts2(int j, int l) const;
    //! PartitionLegs with both photons j and l left out (balanced with k_j+k_l)
    bool PartitionLegs2(int j, int l, ATOOLS::Vec4D_Vector &pb) const;
    void AddExchangeLineShifts2(int j, int l, PropShifts &sh) const;
    //! Is the pair (j,l) eligible for beta_2 in this event (cut, enumeration)?
    bool Beta2PairEligible(int j, int l) const;
    //! eq. (95) of hep-ph/0006359 for one pair in the current partition
    bool ComixInfraredSubtracted_2_0(int j, int l);
    //! n = 2 closure: rho_2 against |M_2|^2 on the event (CEEX: BETA2_CLOSURE)
    void Beta2Closure();
    //! single and double soft limits of beta_2 (CEEX: BETA2_SOFT_TEST)
    void Beta2SoftLimit();
    double m_result1 = 0.;     //!< rho_1 = sum|A_1|^2/4, always
    double m_result01o1 = 0.;  //!< Re<A_0,A_1>/4, always
    double m_result2 = 0.;     //!< rho_2 = sum|A_2|^2/4 (0 unless ORDER 2)
    double m_result12 = 0.;    //!< Re<A_1,A_2>/4
    double m_result02 = 0.;    //!< Re<A_0,A_2>/4
    double m_resultV = 0.;     //!< |A_v|^2/4: what an external v multiplies, squared
    inline int Order() const { return m_order; }
    inline double GetResult1() const { return m_result1; }
    inline double GetResult01Order1() const { return m_result01o1; }
    inline double GetResult2() const { return m_result2; }
    inline double GetResult02() const { return m_result02; }
    inline double GetResult12() const { return m_result12; }
    //! the amplitude an external virtual (v/2) multiplies, |A_v|^2/4: A_0,
    //! or A_1 under CEEX: ORDER2_VIRTUAL_ON_BETA1 (KKMC's (1+v/2) on beta_1)
    inline double GetResultVV() const { return m_resultV; }
    //! Per-photon soft factors, both attachments, at the drawn helicity.
    void CalculateSfactors();
    void MakePhotonHel();
    void MakeBoxMandelstams(const ATOOLS::Vec4D &PX);
    Complex CnuA(double svar, double m1, double m2) const;
    Complex BoxGG();
    Complex BoxGZ();
    Complex BoxSubtract();
    Complex Sfactor(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int hel);
    /*!
      The U / V transition matrices, eq.(222) of hep-ph/0006359.

    */
    Complex UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h2, int i, int j,
                   bool NegMass=0, double M1=-1., double M2=-1.);
    Complex VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h2, int i, int j,
                   double M1=-1., double M2=-1.);
    virtual inline bool ZOnly() {return m_onlyz;}

    Complex S(const Vec4D &p1, const Vec4D &p2, int h1, int h2);
    Complex S(const Vec4D &p1, const Vec4D &p2, double m1,double m2, int h1, int h2);
    Complex T(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1=1, int C2=1,
              double M1=-1., double M2=-1.);
    Complex Tp(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1=1, int C2=1,
               double M1=-1., double M2=-1.);
    Complex U(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1=1, int C2=1,
              double M1=-1., double M2=-1.);
    Complex Up(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1=1, int C2=1,
               double M1=-1., double M2=-1.);
    Complex BornAmplitude(const Vec4D_Vector &k);
    Complex BornAmplitude(const Vec4D_Vector &k, int h0, int h1, int h2, int h3);
    Complex BornAmplitude(Vec4D p1, Vec4D p2, Vec4D p3, Vec4D p4, int h0, int h1, int h2, int h3);
    /*!
      \param Mf3,Mf4 spinor masses for the two FINAL legs; negative means
      "from the flavours". They have to be settable because the beta_1^0 arms
      put a PHOTON in one of those slots, and a photon's spinor mass is zero,
      not sqrt(k^2) of a numerically off-shell momentum.
    */
    void BornAmplitude(const Vec4D_Vector &k, Amplitude &T,
                       double Mf3=-1., double Mf4=-1., int slot=-1);
    /*!
      Spinor structures cached per CALL SITE for the duration of one event.

      m_Tamp/m_Uamp depend only on the momenta and the masses; of everything
      BornAmplitude does, only m_TC/m_UC - two entries, through the two
      propagators - depend on the partition. Without this the partition loop
      recomputed all sixteen spinor products 2^n times to change two numbers.
      Exact, not an approximation: the cached values are bit-for-bit the ones
      that would have been recomputed.

      slot 0 is the Born; photon j owns 1+4j..4+4j, its two ISR arms and its
      two FSR arms (which differ - the ISR arms put the photon in a beam slot,
      the FSR arms in an outgoing one).
    */
    struct SpinorSet { Complex T[2][2][2][2], U[2][2][2][2]; };
    std::vector<SpinorSet> m_spincache;
    std::vector<char>      m_spinvalid;
    void UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpU,
                double M1=-1., double M2=-1.);
    void VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpV,
                double M1=-1., double M2=-1.);

    void SumAmplitude(Complex &sum, const Amplitude &Amp, const Complex fac=1.);
    void AddU(Complex &sum, const Amplitude &Amp, const Amplitude &U, const  Complex fac=1.);
    //! Helicity-resolved forms of AddU/AddV: the same contraction, but kept
    //! per helicity instead of collapsed to a scalar. KKMC's AmpAddI.
    void AddU(Amplitude &out, const Amplitude &Born, const Amplitude &U, const Complex fac=1.);
    void AddV(Amplitude &out, const Amplitude &Born, const Amplitude &V, const Complex fac=1.);
    /*!
      The FINAL-state counterparts, KKMC's AmpAddF. The contraction is the
      mirror image of AmpAddI: U's FIRST index is the output for slot 2 and
      V's SECOND index the output for slot 3, where on the initial side it is
      U's second index for slot 0 and V's first for slot 1.
    */
    void AddUF(Amplitude &out, const Amplitude &U, const Amplitude &Born, const Complex fac=1.);
    void AddVF(Amplitude &out, const Amplitude &Born, const Amplitude &V, const Complex fac=1.);
    void AddV(Complex &sum, const Amplitude &Amp, const Amplitude &V, const  Complex fac=1.);
    void SumAmplitude(Complex &sum, const Amplitude &Amp1, const Amplitude &Amp2, const  Complex fac=1.);

    Complex Splus(const Vec4D &p1, const Vec4D &p2);
    Complex Sminus(const Vec4D &p1, const Vec4D &p2);
    Complex CouplingZ(double i, int mode);
    Complex CouplingG();
    void MakeProp();
    //! t-channel propagators for Bhabha. Spacelike, so no Z width.
    void MakePropT(const ATOOLS::Vec4D_Vector &p);
    //! Electroweak form factors at this partition's scale and angle.
    //! No-op (tree values) unless CEEX: WEAK is set.
    void MakeEWFF(double svar, double costhd);
    void Reset();
    
    inline double GetResult()  {return m_result;}
    //! O(alpha^0) CEEX distribution, i.e. beta_0^0 only.
    inline double GetResult0() {return m_result0;}
    inline double GetResult01() const {return m_result01;}
    /*!
      YFS: CEEX_Virtual: helicity (Ceex_LoopHelicity.C).

      LoopHelicityToCeex puts the loop provider's per-helicity virtual factors
      delta (NLO_Base::m_ldelta, provider labels hel = physical helicity per
      leg, tree amplitudes a0, at the Born point p) into CEEX's helicity index,
      d[f]; entries the provider does not fill keep vhalf. The label map is
      DERIVED, like DeriveComixMap: at the first VIRT_HEL_NCAL events the
      normalised |a0|^2 pattern is matched against Comix's Born at p over all
      leg-flip masks (CEEX: VIRT_HEL_MASK >= 0 overrides), and Comix's index
      goes to CEEX's through the fermion bits of m_comixflip, as in
      InfraredSubtractedME_0_0. Only the labels are needed: the ratio
      A1/A0 cancels every phase convention.

      HelicityVirtualRho(d) = 1/4 sum_f |A_1[f] + d[f] A_0[f]|^2 on this
      event's coherent amplitudes (A_0 = m_AmpExpo0, the partition-summed
      Born, so d[f] multiplies every partition's beta_0 alike). With
      d[f] = v/2 it equals m_result + v m_result01 + v^2/4 m_result0, the
      external form.
    */
    int    m_loopmask = -1, m_loopcal = 0;
    double m_loopmet = -1., m_loopnext = -1., m_looprmin = 0., m_looprmax = 0.;
    long   m_lh_n = 0, m_lh_fail = 0, m_lh_mapbad = 0;
    double m_lh_sumhel = 0., m_lh_sumext = 0., m_lh_sumd = 0., m_lh_sumd2 = 0.;
    double m_lh_regmax = 0., m_lh_dmax = 0.;
    bool   LoopHelicityToCeex(const std::vector<std::vector<int> > &hel,
                              const std::vector<Complex> &a0,
                              const std::vector<Complex> &delta,
                              const ATOOLS::Vec4D_Vector &p, Complex vhalf,
                              std::vector<Complex> &d);
    /*!
      With P1, P2 != 0 the helicity sum is weighted by (1 + l1 P1)(1 + l2 P2),
      l_i the PHYSICAL helicity of beam i in CEEX index f (through the derived
      label map; for f fbar the CEEX labels are physical, KKMC's): the
      longitudinally polarised rho for the unpolarised crude.
    */
    double HelicityVirtualRho(const std::vector<Complex> &d,
                              double P1 = 0., double P2 = 0.) const;
    //! physical helicity (+-1) of leg i in CEEX index f, 0 if unknown
    int    PhysicalHelicity(int f, size_t i) const;
    /*!
      Per event, CEEX_Virtual: helicity: the CEEX correction factor (before
      the IF form-factor cut) with the helicity virtual, m_lhcol[0], and the
      variants the extra named weights report, relative to it:
      [1] external virtual, [2]/[3] beam 1 fully left/right polarised,
      [4]/[5] the same with the external virtual. Reset to 0 per event.
    */
    double m_lhcol[6] = {0., 0., 0., 0., 0., 0.};
    //! run totals: the helicity and external additions to rho_1/rho_crude
    //! and the in-code regression |rho_V(v/2) - external| / rho_1
    void   AccumulateLoopHelicity(bool ok, double addhel, double addext,
                                  double reg, const double split[3]);
    /*!
      rho_V(d) - rho_V(v/2) split, with D = d - v/2, into the part linear in
      D, 1/4 sum 2 Re[(A_1 + v/2 A_0)^* D A_0], and the quadratic ones,
      1/4 sum (Re D)^2 |A_0|^2 and 1/4 sum (Im D)^2 |A_0|^2.
    */
    void   HelicityVirtualSplit(const std::vector<Complex> &d, double vhalf,
                                double split[3]) const;
    double m_lh_split[3] = {0., 0., 0.};
    void   ReportLoopHelicity();
    //! KKMC's RhoCrud: the incoherent partition sum. See m_rhocrud.
    inline double GetRhoCrude() const { return m_rhocrud; }
    /*
      YFS: REAL_BORN_PHOTON_MULTICHANNEL. With photons in the Born the same
      final state is reached through several assignments of photons to Born
      pair and ISR; the generator's density is the SUM of their crudes
      (NLO_Base::BornPhotonChannelSum returns G = sum/generated). Applied by
      YFS_Handler right after Calculate(), before anything reads the crude.
    */
    inline void ScaleRhoCrude(const double g) { m_rhocrud *= g; }
    /*!
      True on the ONE event whose momenta were written to ceex_point.dat.
      CEEX_Compare uses it to mark its line for that event, so the three
      virtuals - CEEX, the external provider and KKMC run on the dumped point -
      are all quoted for the same configuration.
    */
    inline bool JustDumped() const { return m_justdumped; }
    /*!
      rho(Born+virtual)/rho(Born): the O(alpha) virtual correction factor, with
      the real emission excluded. Returns 1 if the Born vanished, so a caller
      can multiply by it unconditionally.
    */
    inline double VirtualFactor() const
    { return m_result0 > 0. ? m_resultbv/m_result0 : 1.; }
    //! rho(Born+real)/rho(Born): the O(alpha) REAL factor, virtual excluded.
    inline double RealFactor() const
    { return m_result0 > 0. ? m_resultbr/m_result0 : 1.; }
    //! Photons in the order CEEX used them: ISR first, then FSR.
    inline const ATOOLS::Vec4D_Vector &AllPhotons() const { return m_allphotons; }
    //! The same photons BEFORE the boost into the CEEX frame. m_allphotons is
    //! boosted, so it cannot be matched against anything the handler holds.
    inline const ATOOLS::Vec4D_Vector &AllPhotonsLab() const { return m_allphotons_lab; }
    //! rho(Born + photon j's real)/rho(Born) - 1, for that photon alone.
    double RealFactorPhoton(size_t j) const;
    inline void SetISRPhotons(Vec4D_Vector &k) {m_isrphotons=k;}
    inline void SetFSRPhotons(Vec4D_Vector &k) {m_fsrphotons=k;}
    inline void SetBornMomenta (Vec4D_Vector &p) {m_bornmomenta =p;}
    inline void SetBorn (const double &born) {m_born=born;}
    inline void SetISRFormFactor (double &form) {m_isrformfactor =form;}
    inline void SetProcessInfo(PHASIC::Process_Info  const &pi) {m_pi = pi;}
    inline bool CheckXs() {return m_checkxs;}
    /*!
      Diagnostics for the t-channel (Bhabha) investigation: the two physical
      t-channel invariants of the CEEX legs and the polar cosine of the
      hardest photon in the CEEX frame. (m_bornmomenta is NOT the event's
      Born point - Phase_Space_Point regenerates the final state after the
      ISR step - so nothing here compares against it.)
    */
    inline double T02() const
    { return m_pceex.size() > 3 ? (m_pceex[0]-m_pceex[2]).Abs2() : 0.; }
    inline double T13() const
    { return m_pceex.size() > 3 ? (m_pceex[1]-m_pceex[3]).Abs2() : 0.; }
    inline double HardestPhotonCosTheta() const
    { double e(-1.), c(0.);
      for (size_t j(0); j < m_allphotons.size(); ++j)
        if (m_allphotons[j][0] > e) { e = m_allphotons[j][0]; c = m_allphotons[j].CosTheta(); }
      return c; }
    double Xi(const Vec4D p, const Vec4D q);
  };
}

#endif
