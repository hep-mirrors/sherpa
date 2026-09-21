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

namespace YFS {

  typedef std::map<int,double> HelicityMap;
  typedef complex<double> dcmplx;

  class Amplitude{
  public:
    Amplitude();
    ~Amplitude(){}
    Complex m_A[2][2][2][2];
    Complex m_U[2][2];
    Complex m_V[2][2];
  };

  class Ceex_Base : public YFS::YFS_Base {  
  public:

    Ceex_Base(const Flavour_Vector &flavs);
    virtual ~Ceex_Base() { ReportBornShape(); }

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
    bool ComixBornAmplitude(const ATOOLS::Vec4D_Vector &p, Amplitude &A);
    Amplitude  m_bornAmp, m_Amp1U, m_Amp1V;
    /*!
      Cumulative spin amplitudes, one per order in alpha, exactly as KKMC's
      m_AmpExpo0/1 (KKceex.cxx:1185). The betas are ADDED INTO these, helicity
      indices intact, and only at the very end is the distribution formed as
      sum_hel |A(hel)|^2 (MakeRho, KKceex.cxx:1234).

      This replaced a scalar accumulation in which each beta was collapsed over
      helicities by SumAmplitude() and the total then squared - i.e.
      |sum_hel sum_beta|^2 instead of sum_hel |sum_beta|^2. Helicity
      configurations are distinguishable, so that summed interference between
      states that cannot interfere.
    */
    Amplitude  m_AmpExpo0, m_AmpExpo1;
    //! Per-helicity snapshots of the three additive pieces of m_AmpExpo1,
    //! taken as each is added: the Born term, the O(alpha) virtual increment
    //! (beta_0^1) and the real increment (beta_1^0). Lets the comparison
    //! against KKMC be made piece by piece instead of only on the total.
    Amplitude  m_snapBorn, m_snapVirt, m_snapReal;
    /*!
      Born + virtual, with no real emission, accumulated over the SAME
      partition sum as m_AmpExpo0/1. Its ratio to m_AmpExpo0 is the O(alpha)
      VIRTUAL factor on its own, which is what NLO_Base wants when CEEX is
      supplying the virtual in place of a loop ME.

      Kept separate from m_snapVirt, which is a per-partition ASSIGNMENT for
      the diagnostics and therefore only ever holds the last partition.
    */
    Amplitude  m_AmpBornVirt;
    //! Born + real, no virtual. Same partition sum; the mirror of m_AmpBornVirt.
    Amplitude  m_AmpBornReal;
    /*!
      Per-photon real emission: m_realphot[j] is photon j's beta_1^0, summed
      over the partitions but NOT over the other photons. Lets the real be
      compared against the external provider photon by photon - the external
      stores its own in YFS::Photon::m_beta10 - instead of only in the total,
      where a disagreement cannot be attributed to a particular emission.

      Their sum is NOT the total real: rho is |B + sum_j R_j|^2, so the cross
      terms between photons belong to no single j.
    */
    std::vector<Amplitude> m_realphot;
    Vec4D_Vector m_allphotons_lab;
    double     m_resultbv = 0., m_resultbr = 0.;
    /*!
      The INCOHERENT partition sum of the Born, KKMC's RhoCrud: each
      partition's beta_0^0 squared and added, rather than the amplitudes added
      and then squared.

      This is the denominator the CEEX weight needs. rho0 is the COHERENT sum,
      so the partitions can cancel and it can approach zero - while the crude
      weight it divides is a factorised eikonal times Born and does not. At
      250 GeV that turned one event into a 1107 pb spike on a 4.86 pb cross
      section. A sum of positive terms cannot do that.
    */
    double     m_rhocrud = 0.;
    Vec4D_Vector   m_momenta, m_isrphotons, m_bornmomenta;
    //! Final-state photons, and the two lists concatenated (ISR first). The
    //! partition flag, the helicity and the two S-factors are all indexed
    //! against m_allphotons.
    Vec4D_Vector   m_fsrphotons, m_allphotons;
    std::vector<int>     m_isrflag;
    std::vector<Complex> m_Sfac_ini, m_Sfac_fin;
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
      that needs no external normalisation, which matters because there is no
      valid external comparator here - Sherpa's m_born is a differential that
      turns out to carry NO angular dependence at all (R^2 = 8e-05 against
      cos(theta)), so an earlier version of this check measured the comparator
      rather than the amplitude.
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
    double     m_result, m_result0, m_e;
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
    void InfraredSubtractedME_1_0(const Vec4D &k, int hel,
                                  const Complex &sProd, const Complex &Sactu,
                                  int slot=-1, int iphot=-1);
    /*!
      beta_1^0 for ONE photon attached to the FINAL state, KKMC's HfinPlus.

      CKine = (q1+q2+k)^2 / (q1+q2)^2. beta_0^0 carries svarX/svarQ (KKMC's
      Cfact0), which is 1 for a purely initial-state partition but not once a
      photon sits on the final line; the sProd*(1-CKine)*Born term here is
      what compensates it, and it vanishes in the infrared as CKine -> 1.
    */
    void InfraredSubtractedME_1_0_FSR(const Vec4D &k, int hel,
                                      const Complex &sProd, const Complex &Sactu,
                                      double CKine, int slot=-1, int iphot=-1);
    /*!
      The ISR/FSR partition of the photons, KKMC's GPS_PartitionStart and
      PartitionPlus. m_isrflag[j] is 1 if photon j hangs off the initial state
      in the CURRENT partition and 0 if it hangs off the final state; the
      amplitude is summed coherently over all 2^n of them, and THAT coherent
      sum is what real initial-final interference is in CEEX.

      last: 0 = keep going, 1 = this is the last one, 2 = stop immediately.
    */
    void PartitionStart(int &last);
    void PartitionPlus(int &last);
    void InfraredSubtractedME_2_0();
    Complex BetaDouble_2_0(Vec4D &k1, Vec4D &k2, int h1, int h2);
    Complex BetaSingle_2_0(Vec4D &k1, Vec4D &k2, int h1, int h2);
    Complex BetaRest_2_0(Vec4D &k1, Vec4D &k2, int h1, int h2);
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

      M1/M2 are the SPINOR masses of p1/p2; negative means "take it from the
      momentum". They matter only in the helicity-FLIPPING element, which is
      the one proportional to the masses - and one of these legs is normally
      the PHOTON, whose Vec4D::Mass() is the square root of whatever numerical
      off-shellness the momentum happens to carry (measured: p^2 = -1.8e-15,
      so Mass() = 4.2e-08) rather than zero. KKMC passes p1.M/p2.M, so its
      photon leg carries MasPhot and its fermion leg the nominal mass.
    */
    Complex UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h2, int i, int j,
                   bool NegMass=0, double M1=-1., double M2=-1.);
    Complex VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h2, int i, int j,
                   double M1=-1., double M2=-1.);
    Complex bsigma(Vec4D k, Vec4D p, int hel);
    virtual inline bool ZOnly() {return m_onlyz;}

    Complex S(const Vec4D &p1, const Vec4D &p2, int h1, int h2);
    Complex S(const Vec4D &p1, const Vec4D &p2, double m1,double m2, int h1, int h2);
    /*!
      GPS spinor structures. C1/C2 are the spinor-type signs of the two legs.

      M1/M2 are the SPINOR masses. A negative value means "take it from the
      momentum", i.e. Vec4D::Mass() = sqrt(p^2), which is what every call site
      used to do implicitly. Pass them explicitly to use the nominal flavour
      mass instead: the momenta are only on shell to about 3e-7 here, and the
      mass enters the helicity-conserving (mass-suppressed) branch linearly.
    */
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
    //! KKMC's RhoCrud: the incoherent partition sum. See m_rhocrud.
    inline double GetRhoCrude() const { return m_rhocrud; }
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
    double Xi(const Vec4D p, const Vec4D q);
  };
}

#endif
