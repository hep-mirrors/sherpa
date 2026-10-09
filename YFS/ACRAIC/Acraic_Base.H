#ifndef YFS_ACRAIC_Acraic_handler
#define YFS_ACRAIC_Acraic_handler

#include "ATOOLS/Phys/Flavour.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "EXTAMP/Process.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/ACRAIC/Acraic_Amplitude.H"


using namespace MODEL;
using namespace ATOOLS;
using namespace std;
using namespace METOOLS;

namespace METOOLS { class Spin_Amplitudes; }

namespace YFS {

  class NLO_Base;   // MapMomenta only; not owned


  class Real_Correction;
  struct VP_Factors;  // YFS/Main/Dyson_VP.H
  struct RVHelicityFactors;  // YFS/NLO/NLO_Base.H

  typedef std::map<int,double> HelicityMap;
  typedef complex<double> dcmplx;

  class Acraic_Base : public YFS::YFS_Base {  
  public:

    Acraic_Base(const Flavour_Vector &flavs);
    virtual ~Acraic_Base() { ReportBornShape(); ReportComixReal(); ReportBeta2(); ReportLoopHelicity();
                             ReportRealVirtualHelicity(); }

    Flavour_Vector m_flavs;
    //! Sherpa processes supplying the Born and one-photon real spin amplitudes
    PHASIC::Process_Base *p_bornproc = nullptr;
    PHASIC::Process_Base *p_realproc = nullptr;
    inline void SetBornProc(PHASIC::Process_Base *p) { p_bornproc = p; }
    inline void SetRealProc(PHASIC::Process_Base *p) { p_realproc = p; }
    //! Born spin amplitudes from Comix; false if unavailable.
    bool ComixBornAmplitude(const ATOOLS::Vec4D_Vector &p, Amplitude &A,
                            double *me2 = NULL,
                            const double propscale = -1.0,
                            const double decscale = -1.0);
    //! Per-event quantities the weight probe reports; see ACRAIC: WEIGHT_PROBE.
    inline size_t NPhot()  const { return m_allphotons.size(); }
    inline double SvarQ()  const { return m_svarQ; }
    inline double SvarY()  const { return m_svarY; }
    //! Legs and photons of the event (ACRAIC: TRACE_FACTOR_ABOVE)
    void PrintKinematics(std::ostream &o) const;
    inline void SetLabMomenta(const ATOOLS::Vec4D_Vector &p) { m_plabmom = p; }
    //! final legs before FSR (YFS m_reallab), where m_born sits
    inline void SetPreFSRMomenta(const ATOOLS::Vec4D_Vector &p) { m_prefsr = p; }
    //! the Comix leg mask of the outgoing fermion pair
    inline size_t DecayCId() const
    { return LegUnion(LegBit(m_if1), LegBit(m_if2)); }

    //! ACRAIC: COMIX_REAL: the one-photon real from Comix's helicity amplitudes
    bool m_comixreal = false;
  

    bool m_hasvirtual = false;
    inline void SetHasVirtual(bool h) { m_hasvirtual = h; }
   
    inline bool AcraicOwnVirtual() const
    { return m_hasvirtual && m_acraicvirtsrc == acraicvirt::acraic; }

    std::vector<YFS::Real_Correction*> m_realprov;
    //! not owned; for MapMomenta (beta_1(k_j) without the other photons)
    inline void SetNLO(YFS::NLO_Base *p) { p_nlo = p; }

    inline void SetRealProvider(size_t n, YFS::Real_Correction *p)
    { if (m_realprov.size() <= n) m_realprov.resize(n+1, NULL);
      m_realprov[n] = p; }
    inline YFS::Real_Correction *RealProvider(size_t n) const
    { return n < m_realprov.size() ? m_realprov[n] : NULL; }
    Amplitude m_comixM1;
    Amplitude m_comixraw[2];
    int m_comixdrawnhel = 0;
    //! ACRAIC: COMIX_REAL_FLIP: bit i if Comix labels leg i's helicity opposite
    int m_comixflip = -1;
    //! fermion-leg flips for n legs: ACRAIC's f is Comix's FlipLegs(f, this)
    inline int ComixLegFlips(size_t n) const { return LowLegs(m_comixflip, n); }
    //! 1 if m_comixflip reverses the photon helicity (leg n), else 0
    inline int ComixPhotonFlip(size_t n) const { return LegBitOf(m_comixflip, n); }
    //! ACRAIC: COMIX_REAL_PHOTON_FLIP: the photon bit the Born cannot fix
    int m_comixphoflip = 1;
    //! ACRAIC: COMIX_REAL_NORM: Comix's amplitude in ACRAIC's normalisation
    double m_comixnorm = -1.;

    bool m_comixcalibrated = false;
    //! result of the Born-level calibration of the Comix-to-ACRAIC map
    struct ComixCalib {
      int    mask = -1;     //!< winning fermion flip mask
      double met = -1.;     //!< its shape metric
      double next = -1.;    //!< the runner-up's, i.e. how decisive the win is
      double N = 0.;        //!< sqrt(sum|A_comix|^2 / sum|e^2 A_hand|^2)
      double sh = 0., sc = 0., me2 = 0.;
      double rmin = 0., rmax = 0.;   //!< per-helicity |C/H| spread
      int    nlive = 0;     //!< helicity slots carrying real weight
    };
    //! per-photon, per-helicity factor turning the hand-coded M_1 into Comix's
    std::vector<Amplitude> m_cxratio;
    //! Comix's beta_1 per photon, compared with m_realphot
    std::vector<Amplitude> m_cxcsub;
    std::vector<char>      m_cxratiook;
    bool m_perphoton = false;
    YFS::NLO_Base *p_nlo = nullptr;
    //! per-photon accounting, kept apart from the per-event counters
    long   m_cxppn = 0, m_cxppfail = 0;
    double m_cxppdev = 0.;
    void BuildComixPhotonRatios();
    //! |beta_1|/beta_0 against E_gamma, hand-coded and Comix
    void Beta1Scan();
    //! Per-helicity spread of the virtual factor over one event's partitions.
    std::vector<double> m_vpmin, m_vpmax;
    bool m_vpon = false;
    //! per-helicity factor turning Comix's Born into ACRAIC's (once per event)
    Amplitude m_cxbalign;
    bool      m_cxbalignok = false;
    bool      m_comixborn = false;
    //! m_comixnorm is the exact 2 -> 2 constant; else RealNorm uses m_cxrnorm
    bool      m_normexact = false;

    double    m_svarY = -1.;
    ATOOLS::Vec4D m_PXvec;   //!< this partition's X, the momentum at s'
    //! the generator's own phase-space point at s'
    ATOOLS::Vec4D_Vector m_plabmom;
    ATOOLS::Vec4D_Vector m_prefsr;
    //! the generator's Born point at R^2 from the pre-FSR legs (any final state)
    bool GeneratorBornAt(const ATOOLS::Vec4D &R, ATOOLS::Vec4D_Vector &pb) const;
    bool RepairMomentumBalance();
    //! the conservation residual found before the repair (FRAME_PROBE)
    double m_repairresid = -1.;
    //! final state is one f fbar pair, the case the hand-coded 2 -> 2 Born covers
    bool m_ffbar = false;
    //! beta_1 partition spread, see ACRAIC: BETA1_PARTITION_CHECK
    double    m_b1min = 0., m_b1max = 0., m_b1sum = 0., m_b1sq = 0.;  int m_b1n = 0;
    double    m_cxrnorm = -1.;
    void BuildComixBornAlignment();
    //! Normalisation for Comix REAL amplitudes: exact constant, else per event
    inline double RealNorm() const
    { return (m_normexact || !(m_cxrnorm > 0.)) ? m_comixnorm : m_cxrnorm; }
    //! e(alpha_YFS)/e(alpha_model) per photon of a Comix real amplitude
    inline double ComixPhotonCoupling(size_t nphot = 1) const
    { return std::pow(m_rescale_alpha, 0.5*(double)nphot); }
    bool ComixRealAt(const ATOOLS::Vec4D_Vector &pp, int hel, Amplitude &M1,
                     const double propscale = -1.0);
    //! ACRAIC: BETA1_CLOSURE: beta_0 + sum_i beta_1(k_i)/s_i vs partition sum
    void Beta1Closure();
    void SoftLimitTest();
    void SoftNormCheck(const ATOOLS::Vec4D &k);
    //! ACRAIC: BORNNORM_CHECK: the ACRAIC Born against m_born
    void CheckBornNorm();
    //! ACRAIC: FRAME_PROBE: the event in its frames
    void ProbeFrames();
    //! YFS: ME_PROBE - rho_1, sum |M_1|^2 and the crude pieces for one-photon events.
    void ProbeMakeRho(const int nh);
    //! ACRAIC: BORN_ALIGN_CHECK - whether the Comix Born alignment built at sQ still holds at s.
    void CheckBornAlignment(const Amplitude &H, const Amplitude &C,
                            const int nh, const double sp_save);
    //! ACRAIC: PARTITION_BORN_CHECK - this partition's Born against BornLegsAt's and Comix's.
    void CheckPartitionBorn(const Amplitude &AmpBorn);
    //! YFS: ME_PROBE - the momenta ACRAIC hands Comix for its one-photon M_1.
    void ProbeComixMomenta(const ATOOLS::Vec4D_Vector &p,
                           const std::vector<METOOLS::Spin_Amplitudes> *amps,
                           const size_t ng, const size_t nl, const double worst);
    //! ACRAIC: BETA1_CHECK - the Comix beta_1 remainder against the hand-coded one at photon j.
    void Beta1Check(const size_t j, const int nh, const Amplitude &C,
                    const Amplitude &Csub, const Amplitude &H,
                    const Amplitude &Bc, const Complex &Sk);
    //! ACRAIC: BETA1_TRACE / TRACE_FACTOR_ABOVE: whether this event is traced
    void SetBeta1Trace();
    //! ACRAIC: BETA1_SCAN
    void RunBeta1Scan();
    //! ACRAIC: TCHANNEL_REDUCED_BORN: whether this event uses the reduced Born
    void ResolveReducedBorn();
    //! ACRAIC: BETA1_PARTITION_CHECK: one partition's beta_1 norm
    void AccumulateBeta1Partition();
    //! ACRAIC: KKMC_FLUX_EMULATION: KKMC flux streams (diagnostic)
    void EmulateKKMCFlux();
    //! ACRAIC: CLOSURE_1PHOT: |A_1|^2 against the exact one-photon real
    void CheckClosureOnePhoton();
    //! ACRAIC: BETA1_PARTITION_CHECK: the spread of beta_1 over partitions
    void ReportBeta1Partitions();
    //! ACRAIC: BETA1_CHECK: |Csub|/|realphot| per photon
    void CheckBeta1Ratio();
    //! ACRAIC: VIRT_PARTITION_CHECK: the virtual factor across partitions
    void CheckVirtualPartitions(const int nparts);
    //! the partition count against m_nstages^(enumerated photons)
    void CheckPartitionCount(const long nparts);
    //! ACRAIC: BETA1_CLOSURE
    void RunBeta1Closure();
    //! ACRAIC: SOFT_LIMIT_TEST
    void RunSoftLimitTest();
    //! ACRAIC: SOFT_NORM_CHECK on the softest photon
    void RunSoftNormCheck();
    bool ReducedBeams(const ATOOLS::Vec4D &k,
                      ATOOLS::Vec4D_Vector &pb) const;
    bool BornLegsAt(const ATOOLS::Vec4D &X,
                    ATOOLS::Vec4D_Vector &pb) const;
    //! Beams back-to-back at X, the radiating pair at Y, spectators as they are
    bool LegsAt(const ATOOLS::Vec4D &X, const ATOOLS::Vec4D &Y,
                ATOOLS::Vec4D_Vector &pb) const;
    //! legs for photon iphot in this partition, balanced with k
    bool PartitionLegs(int iphot, ATOOLS::Vec4D_Vector &pb) const;
    Complex TotalEikonal(const ATOOLS::Vec4D_Vector &p,
                         const ATOOLS::Vec4D &k, int hel);
    //! the polarisation vector Comix's current carries for photon k at hel
    ATOOLS::Vec4C ComixPolarisation(const ATOOLS::Vec4D &k, int hel) const;
    //! the GPS (KKMC) form of the per-leg eikonal, for comparison only
    Complex SfactorLegGPS(const ATOOLS::Vec4D &p, const ATOOLS::Vec4D &k,
                          int hel);
    //! one stage's eikonal current for photon k (iphot) on the legs p
    Complex StageEikonal(int stage, const ATOOLS::Vec4D_Vector &p,
                         const ATOOLS::Vec4D &k, int hel, int iphot = -1);
    bool ComixInfraredSubtracted_1_0(const ATOOLS::Vec4D &k, int hel,
                                    int iphot);
    //! Comix at the physical legs with shifted propagator momenta (as KKMC)
    typedef PHASIC::Process_Base::Prop_Shifts PropShifts;
    bool ComixBornShifted(const ATOOLS::Vec4D_Vector &p, Amplitude &A,
                          const PropShifts &shifts);
    bool ComixRealShifted(const ATOOLS::Vec4D_Vector &pp, int hel,
                          Amplitude &M1, const PropShifts &shifts);
    //! one (leg mask, photon sum) per stage, photon iphot left out
    PropShifts StageShifts(int iphot) const;
    //! YFS: IFI false: veto diagrams emitting a photon off its partition stage
    void AddPlacementVetoes(const std::vector<int> &photons,
                            PropShifts &sh) const;
    //! YFS: IFI false: why this process cannot have it (empty = it can)
    std::string NoIFIRefusal();
    /*!
      YFS: MODE ISR with charged final-state particles: ISR only. Every photon
      sits on the production stage, and the final-stage emission diagrams are
      vetoed in M_1 and M_2 (AddPlacementVetoes), so beta_1, beta_2 are those
      of the initial-state current alone; no IF boxes, no greco. Set per
      process in BuildStages. (A neutral final state is ISR only anyway.)
    */
    bool m_isronly = false;
    inline bool IsrOnly() const { return m_isronly; }
    //! ACRAIC: RESONANCE_EXPONENTIATION (greco, Acraic_Resonance.C) is on
    inline bool ResonanceExponentiation() const
    { return m_resexpmode != acraicresexp::none; }
    //! a_g of stage g on the legs p (external indices of p)
    double StageResonanceCoefficient(int g, const ATOOLS::Vec4D_Vector &p) const;
    //! the (stage leg mask, a_g) entries of this event on the physical legs
    void BuildResonanceExponents();
    //! append this event's resonance entries to a shift list
    void AddResonanceEntries(PropShifts &sh) const;
    //! the resonant log's share of the provider's virtual at the Born point
    void BuildResonanceBornDelta();
    inline double ResonanceBornDeltaV() const { return m_resexpdv; }
    inline const std::vector<Complex> &ResonanceBornDeltaHelicity() const
    { return m_resexpdh; }
    //! ACRAIC: RESONANCE_EXPONENTIATION_PROBE: stage invariants per partition
    void ResonanceProbe(int partition) const;
    //! the reduced legs as leg shift entries for the exchange lines (TCHANNEL_SHIFT)
    void AddExchangeLineShifts(int iphot, PropShifts &sh) const;
    //! ACRAIC: TCHANNEL_SHIFT: 1 on, 0 off, -1 (default) on for 2 -> 2 only
    bool ExchangeLineShiftsOn() const;
    //! whether the Born has a space-like exchange line (TCHANNEL_REDUCED_BORN -1)
    bool BornHasExchangeLine();
    std::map<const PHASIC::Process_Base*, int> m_exchline;
    //! this event: reduced-point beta_0 and balanced-point M_1 (see above)
    bool m_redborn = false;
    //! ACRAIC: RESONANCE_EXPONENTIATION and this event's state for it
    acraicresexp::code m_resexpmode = acraicresexp::none;
    struct ResonanceStage {
      int    stage;   //!< index into m_stagelegs
      size_t legs;    //!< its Comix leg mask (StageShiftMask)
      double a;       //!< a_g on the physical legs m_pacraic
    };
    std::vector<ResonanceStage> m_resexp;
    double m_resexpdv = 0.;
    std::vector<Complex> m_resexpdh;
    //! ACRAIC: TCHANNEL_MULTIPHOTON (with m_redborn and >= 2 photons)
    tchmultiphoton::code TchannelMultiphotonMode() const;
    bool TchannelMultiphotonActive() const
    { return m_redborn && m_allphotons.size() >= 2
             && TchannelMultiphotonMode() != tchmultiphoton::partition_legs; }
    //! photon iphot's one-photon point (TCHANNEL_MULTIPHOTON >= 1)
    bool OnePhotonScaledLegs(int iphot, ATOOLS::Vec4D_Vector &legs,
                             ATOOLS::Vec4D &kpt) const;
    //! TCHANNEL_MULTIPHOTON 2: this partition's beta_1 terms v_j
    std::vector<Amplitude> m_b1terms;
    //! the factorised higher-order remainder of this partition's beta_1 sum
    void AddFactorisedRemainder();
    //! photons whose one-photon point could not be built (HEAD legs used)
    long m_tmpfail = 0;
    //! final-state radiating groups from YFS's dipoles (empty: one final stage)
    std::vector<std::vector<int> > m_groups;
    int m_initstage = 1;
    inline void SetStageGroups(const std::vector<std::vector<int> > &g)
    { m_groups = g; m_wgroups.clear(); }
    //! a final-state radiating group, with its W for a decay stage (W_STAGES)
    struct StageGroup {
      std::vector<int> externalLegs;
      std::vector<int> resonanceDaughters;
      ATOOLS::Flavour  resonance = ATOOLS::Flavour(kf_none);
      bool HasResonance() const { return !resonanceDaughters.empty(); }
    };
    std::vector<StageGroup> m_wgroups;
    inline void SetStageGroups(const std::vector<StageGroup> &g)
    { m_wgroups = g; m_groups.clear(); }
    //! a reconstructed resonance: its daughters plus its decay stage's photons
    struct Resonance {
      ATOOLS::Flavour  flav;
      std::vector<int> daughters;       //!< external legs, e.g. {2,4}
      int              stagemomIndex;   //!< its entry in m_stagemom
      int              decayStage = -1; //!< the stage it is incoming in
    };
    std::vector<Resonance> m_resonances;
    //! ACRAIC: W_STAGES (off, on, auto: on inside the pole window)
    tristate::code m_wstages = tristate::off;
    //! ACRAIC: W_EIKONAL: W pole at its daughters, or plus decay photons
    weikonal::code m_weikonal = weikonal::daughters;
    inline bool WStagesRequested() const { return m_wstages != tristate::off; }
    inline tristate::code WStagesMode() const { return m_wstages; }
    //! true once BuildStages has put a resonance into the stage table
    inline bool WStagesActive() const { return !m_resonances.empty(); }
    bool IsResonanceLeg(int leg) const;
    const Resonance &ResonanceOf(int leg) const;
    //! the external legs a stage covers, a resonance replaced by its daughters
    std::vector<int> StageExternalLegs(int stage) const;
    //! the Comix leg mask of a stage's shift entry (StageShifts)
    size_t StageShiftMask(int stage) const;
    //! momentum of a final stage's radiating system (daughters or legs)
    ATOOLS::Vec4D StageSystemMomentum(int stage) const;
    //! photons this partition assigns to stage g, leaving iphot and lphot out
    ATOOLS::Vec4D StagePhotonSum(int stage, int iphot, int lphot = -1) const;
    //! numerator and pole momentum of a stage leg's eikonal for photon iphot
    struct LegMomentum { ATOOLS::Vec4D num, pole; };
    struct StageLeg;   // defined below with m_stagelegs
    LegMomentum StageLegMomentum(const StageLeg &l, int iphot,
                                 const ATOOLS::Vec4D_Vector &p) const;
    //! recompute m_Sfac of resonance stages for this partition (W_EIKONAL 1)
    void RecomputeResonanceSfactors();
    //! the eikonal of a whole stage for photon iphot on legs p
    Complex StageCurrent(int stage, int iphot, const ATOOLS::Vec4D &k, int hel,
                         const ATOOLS::Vec4D_Vector &p);
    //! this partition's photons except iphot on the reducing stages, or the rest
    ATOOLS::Vec4D PartitionShift(int iphot, bool reducing) const;
    //! Comix's beta_1 = M_1(k) - S(k) B_0 at propagator scale propscale
    bool ComixBeta1At(const ATOOLS::Vec4D &k, int hel,
                      const double propscale, Amplitude &B1);
    bool ComixOnePhotonAmplitude(const ATOOLS::Vec4D &k, int hel,
                                 ATOOLS::Vec4D_Vector &bp, Amplitude &M1,
                                 ATOOLS::Vec4D &kmap);
    bool HandOnePhotonAmplitude(const ATOOLS::Vec4D_Vector &bp,
                                const ATOOLS::Vec4D &k, int hel, Amplitude &M1);
    void AddInitialStateArms(const ATOOLS::Vec4D_Vector &bp,
                             const ATOOLS::Vec4D &k, int hel, Amplitude &M1);
    void AddFinalStateArms(const ATOOLS::Vec4D_Vector &bp,
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
    //! run totals of the Comix-vs-hand comparison
    double m_cxrnormsum = 0., m_cxrmetsum = 0.;
    long   m_cxrn = 0, m_cxrfail = 0;
    //! events taken and refused by the Comix path, per photon multiplicity
    std::vector<long> m_cxrnbyn, m_cxrfailbyn;
    //! First reason the Comix path declined at each multiplicity.
    std::vector<std::string> m_cxrwhy;
    //! agreement with the hand-coded amplitude, per photon multiplicity
    std::vector<double> m_cxrnormbyn, m_cxrmetbyn;
    void CountComixReal(size_t ng, bool ok);
    Amplitude  m_bornAmp, m_Amp1U, m_Amp1V;
    Amplitude  m_AmpExpo0, m_AmpExpo1;
    Amplitude  m_snapBorn, m_snapVirt, m_snapReal;
    //! this partition's Born without its soft-factor product (beta_1 subtracts it)
    Amplitude  m_partborn0;
    //! this partition's pseudo-flux svarY/svarQ (1 under NO_PSEUDOFLUX)
    double     m_pflux = 1.;
    //! print the beta_1 pieces per partition on the CHECK_XS dump event
    bool       m_b1trace = false;
    Amplitude  m_AmpBornVirt;
    //! Born terms weighted as m_AmpBornVirt (not rho_0), for VirtualFactor()
    Amplitude  m_AmpBornOnly;
    double     m_resultbo = 0., m_incbo = 0.;
    Amplitude  m_AmpBornReal;

    std::vector<Amplitude> m_realphot;
    //! ACRAIC: REAL_VIRTUAL 2: photon j's M_1 terms, partition summed
    std::vector<Amplitude> m_realphotM1;
    Vec4D_Vector m_allphotons_lab;
    double     m_resultbv = 0., m_resultbr = 0.;
    double     m_rhocrud = 0.;
    Vec4D_Vector   m_momenta, m_isrphotons, m_bornmomenta;
    Vec4D_Vector   m_fsrphotons, m_allphotons;
    int                  m_nstages = 2;
    std::vector<int>     m_stage;
    std::vector<char>    m_stagereduces;
    //! photons on a fixed stage (SOFT_PARTITION_CUT); this partition's soft factors
    std::vector<char>    m_fixedstage;
    std::vector<Complex> m_sactu;
    //! the crude's (incoherent) soft-factor product; whether a collapsed photon enters it
    double m_crudeprod = 1.;
    bool   m_crudefixed = false;
    //! ACRAIC: IFI; false sums the partitions incoherently (KKMC KeyINT = 0)
    bool   m_ifi_coherent = true;
    double m_inc00 = 0., m_inc11 = 0., m_inc01 = 0.;
    //! the same at ORDER 2 and for Born+virtual (YFS: IFI false)
    double m_inc22 = 0., m_inc12 = 0., m_inc02 = 0., m_incbv = 0.;
    std::vector<std::vector<Complex> > m_Sfac;
    //! stage leg: m_pacraic index, w = Q theta (-1 in, +1 out); S = sum w_i S_i
    struct StageLeg { int leg; double w; };
    std::vector<std::vector<StageLeg> > m_stagelegs;
    //! momenta the stage legs index: external legs, then resonances
    std::vector<ATOOLS::Vec4D> m_stagemom;
    //! (q1+q2)^2, the physical final-pair invariant (m_sp is the partition's)
    double m_svarQ = 0.;
    int    m_nparts = 0;   //!< partitions summed on the last event
    size_t m_maxnphot = 0; //!< highest photon multiplicity seen
    long   m_partbad = 0, m_partn = 0;
    double m_pzmin = 0., m_pzmax = 0.;
    double m_pzspread[16] = {0.}, m_pzsum[16] = {0.};
    long   m_pzn[16] = {0};
    const double m_pzxcut[4] = {1e-5, 1e-4, 1e-3, 1e-2};
    double m_pzsoft[4] = {0.}, m_pzsoftsum[4] = {0.}, m_pzsoftsav[4] = {0.};
    long   m_pzsoftn[4] = {0};
    Vec4D m_eta, m_zeta, m_b;
    //! lab -> Init() frame (Init() boosts only m_momenta)
    ATOOLS::Poincare m_cms;

    //! legs of the spin amplitudes: the beams and the physical outgoing legs
    ATOOLS::Vec4D_Vector m_pacraic;
    void BuildAcraicMomenta();
    //! m_if1/m_if2 = the first two outgoing fermions (2/3 if none)
    void FindOutgoingFermionPair(const ATOOLS::Flavour_Vector &flavs);

    //! the Born-shape fit a(1 + (b/a) cos + cos^2), accumulated over the run
    double m_fit[3][3] = {{0.}}, m_fitb[3] = {0.}, m_fity = 0., m_fity2 = 0.;
    long   m_fitn = 0;
    //! run totals of the O(alpha^0) and O(alpha) ACRAIC distributions
    double m_rho0sum = 0., m_rho1sum = 0.;
    //! the same split into Born+virtual and Born+real
    double m_rhobvsum = 0., m_rhobrsum = 0.;
    //! run total of the incoherent Born (RhoCrud)
    double m_rhocrudsum = 0.;  long m_rhon = 0;
    bool   m_acraicdumped = false, m_justdumped = false;
    //! Worst relative violation of the ACRAIC soft-factor identity, over the run.
    double m_sfacworst = 0.;
    long   m_sfacn = 0;
    //! |s_sigma(k)|^2 against the closed form, hep-ph/0006359 eq.(soft-fac-isr)
    void   CheckSoftFactor(const ATOOLS::Vec4D &k);
    //! U[p,k_sigma,p] = V[p,k_sigma,p] = b_sigma(k,p) delta_{l1 l2},
    //! hep-ph/0006359 eq.(diagonality)
    void   CheckUVDiagonality(const ATOOLS::Vec4D &k);
    double m_uvoff = 0., m_uvdiag = 0., m_bsig = 0., m_uvflip = 0.;
    long   m_uvn = 0, m_uvflipn = 0;
    void   AccumulateBornShape(double cth, double val);
    void   ReportBornShape();
    bool m_onlyz, m_onlyg;
    bool m_fixedwidth,m_checkxs;
    std::vector<int> m_PhoHel;
    Complex    m_T, m_U, m_sW, m_cW, m_sin2tw, m_gcpl, m_propG, m_propZ, m_cfac;
    Complex m_cL, m_cR, m_prop, m_I_L, m_I_R, m_F_L, m_F_R, m_I;
    //! sign of the eikonal against Comix's M_1 (soft limit M_1 -> -(s C))
    Complex m_sfacphase = Complex(-1., 0.);
    Complex m_Tamp[2][2][2][2], m_Uamp[2][2][2][2],m_ampborn[2][2][2][2], m_Sc, m_Tc, m_Uc;
    Complex m_TC[2], m_UC[2], m_ae, m_af, m_ve, m_vf, m_norm;
    //! Bhabha t-channel coupling x propagator factors
    Complex m_propGt, m_propZt, m_TCt[2], m_UCt[2];
    //! electroweak form factors (KKMC GPS_EWFFact); tree values 1
    Complex m_kapE{1.,0.}, m_kapF{1.,0.}, m_kapEF{1.,0.};
    Complex m_rhoEW{1.,0.}, m_gamVPi{1.,0.}, m_vvcor{1.,0.};
    //! 2 T3 and 4 Q sin^2(theta_W) apart, so kappa multiplies the sin^2 piece
    Complex m_t3e2, m_t3f2, m_qesw, m_qfsw;
    bool    m_weak = false;
    //! t-channel spinor structures, (2,1)x(3,4) and (2,3)x(1,4)
    Complex m_Tampt[2][2][2][2], m_Uampt[2][2][2][2];
    double  m_tinv = 0.;
    bool    m_bhabha = false;
    Complex m_BoxGGtu, m_BoxGZtu, m_BoxGGut, m_BoxGZut, m_Sprod, m_beta00, m_beta10, m_beta01;
    Complex m_vertexI, m_vertexF, m_b1, m_b2;
    Complex m_rawboxGG, m_rawboxGZ, m_rawboxSub, m_rawcoef;
    double     m_Q1Q2I, m_Q1Q2F, m_QIQF, m_MZ, m_gZ, m_isrformfactor, m_sQ;
    double     m_mass_I, m_mass_F, m_qe, m_qf, m_crude;
    //! the outgoing fermion pair's legs (FindOutgoingFermionPair)
    size_t m_if1 = 2, m_if2 = 3;
    double     m_result, m_result0, m_e;
    // Re sum_h A_0(h)^* A_1(h)/4, what an external virtual multiplies
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
    //! helicity index into Amplitude::m_A
    static inline int Idx(int h0, int h1, int h2, int h3)
    { return Amplitude::Idx(h0, h1, h2, h3); }
    Complex SfactorLeg(const ATOOLS::Vec4D &p, const ATOOLS::Vec4D &k, int hel);
    //! per-leg eikonal e (pnum.eps)/(ppole.k), as a shifted Comix propagator gives
    Complex SfactorLeg(const ATOOLS::Vec4D &pnum, const ATOOLS::Vec4D &ppole,
                       const ATOOLS::Vec4D &k, int hel);
    //! derive the charge-neutral emission stages (m_stagelegs, m_nstages)
    bool BuildStages();
    void PartitionStart(int &last);
    void PartitionPlus(int &last);
    //! PartitionStart with the softest photons collapsed (MAX_PARTITION_PHOTONS)
    void PartitionStartCollapsed(int &last);
    //! ACRAIC: CRUDE_FROM_GENERATOR: the crude on the generator's stages (< 0: none)
    double CrudeFromGenerator();
    //! no beta_1 for photon j in any partition (ACRAIC: NEXT_TO_SOFT_FLUX)
    bool Beta1NeverEvaluated(size_t j) const;
    //! photon j's next-to-soft flux: its stage invariant with over without it
    double NextToSoftFlux(size_t j) const;
    std::vector<std::vector<int> > m_crudegroups;
    inline void SetCrudeStageGroups(const std::vector<std::vector<int> > &g)
    { m_crudegroups = g; }
    //! ACRAIC: CRUDE_FROM_GENERATOR: off, on, compare, auto (on with W stages)
    crudegen::code m_crudegen = crudegen::off;
    //! saves the stage tables, restores them on destruction
    struct StageTableSwap {
      Acraic_Base &c;
      std::vector<std::vector<StageLeg> > legs;
      std::vector<ATOOLS::Vec4D> mom;
      std::vector<Resonance> res;
      std::vector<std::vector<Complex> > sfac;
      std::vector<int> stage;
      std::vector<char> fixed, reduces;
      std::vector<StageGroup> wgroups;
      std::vector<std::vector<int> > groups;
      int initstage, nstages;
      StageTableSwap(Acraic_Base &cb);
      ~StageTableSwap();
    };
    //! ACRAIC: ORDER; at 2 beta_2 from Comix's two-photon amplitude
    int m_order = 1;
    Amplitude m_AmpBeta2;     //!< the beta_2 increments alone (partition summed)
    Amplitude m_AmpExpo2;     //!< m_AmpExpo1 + m_AmpBeta2, built in MakeRho
    //! ACRAIC: KKMC_FLUX_EMULATION: sum_{i<j in F} 2 k_i.k_j/Q^2 sProd B (all, < vcut2)
    Amplitude m_AmpFluxAll, m_AmpFluxSoft;
    std::vector<Amplitude> m_b2M1;    //!< aligned M_1(k_j) of this partition
    std::vector<char>      m_b2hasM1; //!< this partition carried M_1(k_j)
    std::vector<char>      m_b2on;    //!< beta_1(k_j) active in this partition
    std::vector<Complex>   m_b2sj;    //!< the eikonal beta_1(k_j) subtracted with
    std::vector<char>      m_b2elig;  //!< photon may pair (enumerated, x cut)
    long m_b2events = 0;
    //! sorted photon indices (tree residuals beta-hat_S)
    typedef std::vector<int> PhotonSet;
    int m_treebetamax = 1;            //!< ACRAIC: TREE_BETA_MAX resolved (1 or ORDER)
    int m_treebetacap = 0;            //!< this event's cap (providers, eligibility)
    std::map<PhotonSet, Amplitude> m_treebeta;   //!< this partition's beta_S, |S| >= 2
    std::vector<long> m_bnterms, m_bnamps;       //!< per |S|: terms, M_|S| evaluations
    long m_b2fail = 0;
    //! ACRAIC: TREE_BETA_PROBE: per order m, sum_{|S| = m} W_S beta_S
    std::vector<Amplitude> m_treeorder;
    bool m_treeprobe = false;
    //! per event: which photons may enter beta_m (m >= 2), and the cap
    bool PrepareBeta2();
    void ReportBeta2() const;
    //! zero the per-partition beta_1 record (called at each partition start)
    void ResetBeta1Record();
    //! Comix M_m in ACRAIC's alignment, hel[q] the q-th photon leg's helicity
    bool ComixRealNHel(const ATOOLS::Vec4D_Vector &pp, const std::vector<int> &hel,
                       Amplitude &M, const PropShifts &shifts, bool noshift);
    //! Comix M_m(S) in ACRAIC's alignment: photon S[q] in the q-th photon slot
    bool ComixRealNShifted(const ATOOLS::Vec4D_Vector &pp, const PhotonSet &S,
                           Amplitude &M, const PropShifts &shifts,
                           bool noshift = false);
    //! the two-photon case of ComixRealNShifted (soft-limit diagnostics)
    bool ComixReal2Shifted(const ATOOLS::Vec4D_Vector &pp, int hj, int hl,
                           Amplitude &M2, const PropShifts &shifts,
                           bool noshift = false);
    //! the photons the current partition puts on stage, except those in S
    ATOOLS::Vec4D StagePhotonSumExcept(int stage, const PhotonSet &S) const;
    //! StageShifts with every photon of S left out
    PropShifts SubsetStageShifts(const PhotonSet &S) const;
    PropShifts StageShifts2(int j, int l) const;
    //! PartitionLegs with every photon of S left out (balanced with sum_S k)
    bool SubsetPartitionLegs(const PhotonSet &S, ATOOLS::Vec4D_Vector &pb) const;
    void AddSubsetExchangeLineShifts(const PhotonSet &S, PropShifts &sh) const;
    //! Is the pair (j,l) eligible for beta_2 in this event (cut, enumeration)?
    bool Beta2PairEligible(int j, int l) const;
    //! the eligible photons of this event, ascending
    PhotonSet TreeBetaPhotons() const;
    //! beta_1 of photon j as this partition's beta_1 pass added it
    Complex RecordedBeta1(int j, int f) const;
    //! the product of this partition's soft factors m_sactu over S
    Complex SoftProduct(const PhotonSet &S) const;
    //! beta_S for one S in the current partition, added with its W_S
    bool ComixInfraredSubtracted_N_0(const PhotonSet &S);
    //! every beta_S of the current partition, 2 <= |S| <= cap, by size
    void AddTreeBetas(bool dobetan);
    //! the m = 0 and m = 1 rungs of the probe ladder for this partition
    void AddProbeLowOrders();
    //! @@@ BNCLOS: |A_N|^2 for N = 0..cap against |M_n|^2 (TREE_BETA_PROBE)
    void TreeBetaClosure();
    //! n = 2 closure: rho_2 against |M_2|^2 on the event (ACRAIC: BETA2_CLOSURE)
    void Beta2Closure();
    //! single and double soft limits of beta_2 (ACRAIC: BETA2_SOFT_TEST)
    void Beta2SoftLimit();
    double m_result1 = 0.;     //!< rho_1 = sum|A_1|^2/4, always
    double m_result01o1 = 0.;  //!< Re<A_0,A_1>/4, always
    double m_result2 = 0.;     //!< rho_2 = sum|A_2|^2/4 (0 unless ORDER 2)
    double m_result12 = 0.;    //!< Re<A_1,A_2>/4
    double m_result02 = 0.;    //!< Re<A_0,A_2>/4
    double m_resultV = 0.;     //!< |A_v|^2/4: what an external v multiplies, squared
    //! ACRAIC: BETA0_2 (delta^(2) on A_0 alone) needs only overlaps with A_0
    double m_result0x = 0.;    //!< Re<A_0,A_X>/4
    double m_result0v = 0.;    //!< Re<A_0,A_v>/4
    inline int Order() const { return m_order; }
    inline double GetResult1() const { return m_result1; }
    inline double GetResult01Order1() const { return m_result01o1; }
    inline double GetResult02() const { return m_result02; }
    inline double GetResult12() const { return m_result12; }
    //! |A_v|^2/4, A_v = A_0, or A_1 under ACRAIC: ORDER2_VIRTUAL_ON_BETA1
    inline double GetResultVV() const { return m_resultV; }
    //! Re<A_0, A_X>/4 and Re<A_0, A_v>/4 (ACRAIC: BETA0_2)
    inline double GetResult0X() const { return m_result0x; }
    inline double GetResult0V() const { return m_result0v; }
    //! ACRAIC: REAL_VIRTUAL (Acraic_Base.C, RealVirtualMode)
    static acraicrv::code RealVirtualMode();
    //! ACRAIC: REAL_VIRTUAL 2: sum_j dv_j Re<A, M_1j>/4, A = A_1 (A_2 at
    //! ORDER 2), dv_j = v_{n+1,j} - v_B per photon (0 where unknown)
    double RealVirtualRemainderRho(const std::vector<double> &dv) const;
    /*!
      ACRAIC: REAL_VIRTUAL helicity (Acraic_RealVirtualHelicity.C): the
      same linear term resolved by helicity,
        sum_j sum_h 2 Re[conj(A^h) r_j^h M_1j^h]/4 ,
      r_j^h the remainder of photon j in ACRAIC's index (complex). With
      r_j^h = dv_j/2 on every h it is RealVirtualRemainderRho(dv).
    */
    double HelicityRealVirtualRemainderRho
      (const std::vector<std::vector<Complex> > &r) const;
    //! photon j's per-helicity remainder in ACRAIC's index at the photon's
    //! drawn helicity phohel; configurations the provider does not give
    //! keep dvhalf. False (r = dvhalf everywhere) without a label map.
    bool RealVirtualHelicityToAcraic(const RVHelicityFactors &rv, int phohel,
                                     Complex dvhalf, std::vector<Complex> &r);
    //! the photon's drawn helicity (m_PhoHel), 0 if j is out of range
    inline int PhotonHelicity(size_t j) const
    { return j < m_PhoHel.size() ? m_PhoHel[j] : 0; }
    //! the (n+1)-point label map and run totals of REAL_VIRTUAL helicity
    int    m_rvmask = -1, m_rvcal = 0;
    double m_rvmet = -1., m_rvnext = -1.;
    long   m_rvh_n = 0, m_rvh_avg = 0, m_rvh_mapbad = 0;
    double m_rvh_sumhel = 0., m_rvh_sumavg = 0., m_rvh_sumd = 0., m_rvh_sumd2 = 0.;
    double m_rvh_dmax = 0.;
    void AccumulateRealVirtualHelicity(bool anyhel, double addhel, double addavg);
    void ReportRealVirtualHelicity();
    //! true when photon j carried an M_1 term in some partition
    bool PhotonHasM1(size_t j) const;
    //! Per-photon soft factors, both attachments, at the drawn helicity.
    void CalculateSfactors();
    void MakePhotonHel();
    void MakeBoxMandelstams(const ATOOLS::Vec4D &PX);
    Complex CnuA(double svar, double m1, double m2) const;
    Complex BoxGG();
    Complex BoxGZ();
    Complex BoxSubtract();
    Complex Sfactor(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int hel);
    //! U / V transition matrices, hep-ph/0006359 eq. (222)
    Complex UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h2, int i, int j,
                   bool NegMass=0, double M1=-1., double M2=-1.);
    Complex VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h2, int i, int j,
                   double M1=-1., double M2=-1.);

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
    //! \param Mf3,Mf4 final-leg spinor masses (< 0: from the flavours; 0 for a photon)
    void BornAmplitude(const Vec4D_Vector &k, Amplitude &T,
                       double Mf3=-1., double Mf4=-1., int slot=-1);
    //! spinor structures cached per call site (slot 0 Born, photon j 1+4j..4+4j)
    struct SpinorSet { Complex T[2][2][2][2], U[2][2][2][2]; };
    std::vector<SpinorSet> m_spincache;
    std::vector<char>      m_spinvalid;
    void UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpU,
                double M1=-1., double M2=-1.);
    void VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpV,
                double M1=-1., double M2=-1.);

    void SumAmplitude(Complex &sum, const Amplitude &Amp, const Complex fac=1.);
    void AddU(Complex &sum, const Amplitude &Amp, const Amplitude &U, const  Complex fac=1.);
    //! helicity-resolved AddU/AddV (KKMC's AmpAddI)
    void AddU(Amplitude &out, const Amplitude &Born, const Amplitude &U, const Complex fac=1.);
    void AddV(Amplitude &out, const Amplitude &Born, const Amplitude &V, const Complex fac=1.);
    //! final-state counterparts (KKMC's AmpAddF): U's first index is the
    //! output for slot 2, V's second for slot 3 (mirror of AmpAddI)
    void AddUF(Amplitude &out, const Amplitude &U, const Amplitude &Born, const Complex fac=1.);
    void AddVF(Amplitude &out, const Amplitude &Born, const Amplitude &V, const Complex fac=1.);
    void AddV(Complex &sum, const Amplitude &Amp, const Amplitude &V, const  Complex fac=1.);
    void SumAmplitude(Complex &sum, const Amplitude &Amp1, const Amplitude &Amp2, const  Complex fac=1.);

    Complex Splus(const Vec4D &p1, const Vec4D &p2);
    Complex Sminus(const Vec4D &p1, const Vec4D &p2);
    Complex CouplingZ(double i, int mode);
    Complex CouplingG();
    //! resonance: ACRAIC: RESONANCE_EXPONENTIATION on the hand-coded Born
    void MakeProp(bool resonance = true);
    //! t-channel propagators for Bhabha. Spacelike, so no Z width.
    void MakePropT(const ATOOLS::Vec4D_Vector &p);
    //! electroweak form factors at this partition's scale (ACRAIC: WEAK)
    void MakeEWFF(double svar, double costhd);
    void Reset();
    
    inline double GetResult()  {return m_result;}
    //! O(alpha^0) ACRAIC distribution, i.e. beta_0^0 only.
    inline double GetResult0() {return m_result0;}
    inline double GetResult01() const {return m_result01;}
    //! YFS: ACRAIC_Virtual: helicity: provider-to-ACRAIC label map, run totals
    int    m_loopmask = -1, m_loopcal = 0;
    double m_loopmet = -1., m_loopnext = -1., m_looprmin = 0., m_looprmax = 0.;
    long   m_lh_n = 0, m_lh_fail = 0, m_lh_mapbad = 0;
    double m_lh_sumhel = 0., m_lh_sumext = 0., m_lh_sumd = 0., m_lh_sumd2 = 0.;
    double m_lh_regmax = 0., m_lh_dmax = 0.;
    bool   LoopHelicityToAcraic(const std::vector<std::vector<int> > &hel,
                              const std::vector<Complex> &a0,
                              const std::vector<Complex> &delta,
                              const ATOOLS::Vec4D_Vector &p, Complex vhalf,
                              std::vector<Complex> &d);
    //! helicity sum weighted by (1 + l1 P1)(1 + l2 P2), l_i physical beam helicities
    double HelicityVirtualRho(const std::vector<Complex> &d,
                              double P1 = 0., double P2 = 0.) const;
    //! physical helicity (+-1) of leg i in ACRAIC index f, 0 if unknown
    int    PhysicalHelicity(int f, size_t i) const;
    //! 1/4 sum_f |A_X[f] + vhalf A_v[f] + e[f] A_0[f]|^2 (YFS: DYSON_VP),
    //! A_X = A_2 at ORDER 2 else A_1, A_v = A_1 if vonbeta1 else A_0
    double A0InsertionRho(const std::vector<Complex> &e, double vhalf,
                          bool vonbeta1) const;
    //! the Dyson VP factor per ACRAIC helicity (Born-weighted mean if unknown)
    std::vector<Complex> HelicityShares(const VP_Factors &vp) const;
    //! values in physical labels on ACRAIC's helicity index (else `fallback`)
    std::vector<Complex> PhysicalLabelsToAcraic(const std::vector<std::vector<int> > &hel,
                                              const std::vector<Complex> &val,
                                              const Complex &fallback) const;
    //! sum_f 2 Re(e[f]) |A_0[f]|^2 / sum_f |A_0[f]|^2 (YFS: DYSON_VP_PROBE)
    double A0WeightedTwoRe(const std::vector<Complex> &e) const;
    //! where the nominal's external v/2 sits: on A_1 (true) or on A_0
    inline bool VirtualOnBeta1() const { return m_vonbeta1; }
    bool m_vonbeta1 = false;
    //! per event, ACRAIC_Virtual: helicity: [0] the correction factor with the
    //! helicity virtual; [1] external virtual, [2]/[3] beam 1 left/right, [4]/[5] both
    double m_lhcol[6] = {0., 0., 0., 0., 0., 0.};
    //! run totals of the helicity and external additions and their regression
    void   AccumulateLoopHelicity(bool ok, double addhel, double addext,
                                  double reg, const double split[3]);
    //! rho_V(d) - rho_V(v/2), D = d - v/2, split into 1/4 sum 2 Re[(A_1 + v/2
    //! A_0)^* D A_0], 1/4 sum (Re D)^2 |A_0|^2 and 1/4 sum (Im D)^2 |A_0|^2
    void   HelicityVirtualSplit(const std::vector<Complex> &d, double vhalf,
                                double split[3]) const;
    double m_lh_split[3] = {0., 0., 0.};
    void   ReportLoopHelicity();
    //! KKMC's RhoCrud: the incoherent partition sum. See m_rhocrud.
    inline double GetRhoCrude() const { return m_rhocrud; }
    //! YFS: REAL_BORN_PHOTON_MULTICHANNEL: G = sum of channel crudes/generated
    inline void ScaleRhoCrude(const double g) { m_rhocrud *= g; }
    //! true on the event whose momenta were written to acraic_point.dat
    inline bool JustDumped() const { return m_justdumped; }
    //! rho(Born+virtual)/rho(Born), 1 if the Born vanished
    inline double VirtualFactor() const
    { return m_resultbo > 0. ? m_resultbv/m_resultbo : 1.; }
    //! rho(Born+real)/rho(Born): the O(alpha) real factor, virtual excluded.
    inline double RealFactor() const
    { return m_result0 > 0. ? m_resultbr/m_result0 : 1.; }
    //! the photons in the lab, before the boost into the ACRAIC frame
    inline const ATOOLS::Vec4D_Vector &AllPhotonsLab() const { return m_allphotons_lab; }
    //! rho(Born + photon j's real)/rho(Born) - 1, for that photon alone.
    double RealFactorPhoton(size_t j) const;
    inline void SetISRPhotons(Vec4D_Vector &k) {m_isrphotons=k;}
    inline void SetFSRPhotons(Vec4D_Vector &k) {m_fsrphotons=k;}
    inline void SetBornMomenta (Vec4D_Vector &p) {m_bornmomenta =p;}
    inline void SetBorn (const double &born) {m_born=born;}
    inline void SetISRFormFactor (double &form) {m_isrformfactor =form;}
    //! t-channel diagnostics: t invariants of the ACRAIC legs, hardest photon's cos
    inline double T02() const
    { return m_pacraic.size() > 3 ? (m_pacraic[0]-m_pacraic[2]).Abs2() : 0.; }
    inline double T13() const
    { return m_pacraic.size() > 3 ? (m_pacraic[1]-m_pacraic[3]).Abs2() : 0.; }
    inline double HardestPhotonCosTheta() const
    { double e(-1.), c(0.);
      for (size_t j(0); j < m_allphotons.size(); ++j)
        if (m_allphotons[j][0] > e) { e = m_allphotons[j][0]; c = m_allphotons[j].CosTheta(); }
      return c; }
    double Xi(const Vec4D p, const Vec4D q);
  };
}

#endif
