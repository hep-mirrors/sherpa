/*!
  \file Ceex_Comix.C

  The O(alpha) REAL correction taken from Comix instead of from hand-coded
  spinor products.

*/

#include "YFS/CEEX/Ceex_Base.H"
#include "YFS/NLO/Real_Correction.H"

#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Flavour.H"

#include <cstdlib>
#include <iostream>

using namespace YFS;
using namespace ATOOLS;

namespace {

  /*!
    CEEX and Comix both index a helicity as 0 -> +1, 1 -> -1, and with every
    external leg two-valued the flat Spin_Amplitudes index is just the bits:
    leg i in bit i, the photon in the bit above the last fermion leg. That is
    also how Amplitude::Idx packs them, so the fermion part of a Comix index
    and a CEEX index are THE SAME NUMBER and no per-leg unpacking is needed.

    The consequence used throughout below: flipping the helicity LABEL of a
    set of legs is an XOR of the packed index with that set's bit mask, so the
    relabelling is one operation rather than a loop over legs.
  */
  inline size_t FlatIndex(size_t ferm, int hg, int nlegs)
  { return ferm | ((size_t)hg << nlegs); }

}


bool Ceex_Base::FetchComixReal()
{
  m_havecomixreal = false;
  const size_t nl(m_flavs.size());
  const size_t ng(m_allphotons.size());
  /*
    Every refusal below says why, once per multiplicity. A bare count of
    refusals is not a diagnosis: "no process was built at this multiplicity"
    and "the momenta do not balance" are different problems with different
    fixes, and without the reason a run that never once used the Comix
    amplitude looks the same as one that used it and disagreed.
  */
  #define CXR_REFUSE(why)                                                     \
    do { if (m_cxrwhy.size() <= ng) m_cxrwhy.resize(ng+1);                    \
         if (m_cxrwhy[ng].empty()) {                                          \
           m_cxrwhy[ng] = (why);                                              \
           msg_Info()<<"CEEX: COMIX_REAL declines "<<ng<<" photons: "          \
                     <<(why)<<" (reported once)\n"; }                         \
         return false; } while (0)

  if (ng < 1 || m_pceex.size() < nl || m_PhoHel.size() < ng)
    CXR_REFUSE("no photons, or momenta/helicities not set");
  /*
    CEEX's own container is NOT the constraint here, and an earlier version of
    this check said it was. Amplitude holds 2^nl entries - one per FERMION
    helicity - because the photon helicities are drawn per event rather than
    stored (Amplitude::SetLegs is called with m_flavs.size()). What has to
    hold 2^(nl+ng) is Comix's Spin_Amplitudes, which is Comix's array and is
    checked against that size below. So the only limit worth testing here is
    that the index does not overflow, and everything else is decided by
    whether a process with nl+ng legs was built at all.
  */
  if (nl + ng >= 8*sizeof(size_t))
    CXR_REFUSE("helicity index would overflow");
  YFS::Real_Correction *prov(RealProvider(ng));
  if (prov == NULL)
    CXR_REFUSE("no real ME provider registered at this multiplicity");
  if (prov->p_proc == NULL)
    CXR_REFUSE("the provider has no process attached");

  // The real process must be this process plus ng photons, in that order -
  // YFS_Process builds it by pushing photons onto the final state, so it is,
  // but a mapped or reordered process would silently misindex every helicity.
  const Flavour_Vector &pf(prov->p_proc->Flavours());
  if (pf.size() != nl + ng)
    CXR_REFUSE("the provider's process has the wrong number of legs");
  for (size_t i(0); i < nl; ++i) if (pf[i] != m_flavs[i])
    CXR_REFUSE("the provider's process has different fermions");
  for (size_t i(nl); i < nl + ng; ++i)
    if (pf[i].Kfcode() != kf_photon)
      CXR_REFUSE("the provider's extra legs are not photons");

  Vec4D_Vector p(nl + ng);
  for (size_t i(0); i < nl; ++i) p[i] = m_pceex[i];
  for (size_t j(0); j < ng; ++j) p[nl + j] = m_allphotons[j];

  // m_pceex uses the LAB pair, which is the only one that balances against
  // the photons (see the comment on m_pceex). Check it rather than trust it:
  // Comix at a non-conserving point returns a number, not an error.
  Vec4D bal(p[0] + p[1]);
  for (size_t i(2); i < p.size(); ++i) bal -= p[i];
  const double scale(sqrt(dabs(m_svarQ)) + 1.);
  double worst(0.);
  for (int c(0); c < 4; ++c) worst = Max(worst, dabs(bal[c]));
  if (worst > 1e-6*scale) { ++m_cxrfail; CXR_REFUSE("momenta do not balance"); }

  const std::vector<METOOLS::Spin_Amplitudes> *amps
    (prov->ComixAmplitudes(p));
  if (amps == NULL || amps->empty())
    { ++m_cxrfail; CXR_REFUSE("Comix returned no amplitudes (KeepAmplitudes?)"); }
  if (amps->size() > 1) {
    static bool warned(false);
    if (!warned) {
      warned = true;
      msg_Error()<<METHOD<<"(): Comix returned "<<amps->size()
                 <<" colour configurations; CEEX uses the first. Reported "
                 <<"once."<<std::endl;
    }
  }
  const METOOLS::Spin_Amplitudes &sa((*amps)[0]);
  /*
    Count the legs that CARRY two helicity states, not the legs. METOOLS packs
    mixed-radix over the per-leg spin multiplicities, so a scalar contributes a
    factor one: e+e- -> H mu+ mu- gamma has six legs but 2x2x1x2x2x2 = 32
    entries, not 2^6 = 64. Amplitude::s_nlegs is already that count for the
    fermion side (SetLegs takes the flavours), and each photon adds one.
  */
  const size_t np(Amplitude::s_nlegs);
  const size_t nhel(((size_t)1) << (np + ng));
  if (sa.size() != nhel) {
    static bool warned(false);
    if (!warned) {
      warned = true;
      msg_Error()<<METHOD<<"(): expected "<<nhel<<" helicity entries for "
                 <<(np+ng)<<" packed legs, got "<<sa.size()
                 <<". CEEX: COMIX_REAL disabled for this run."<<std::endl;
    }
    m_comixreal = 0;
    return false;
  }

  /*
    CEEX draws a helicity per photon (MakePhotonHel) and both the numerator
    and the crude denominator use the drawn ones, so the average is right;
    Comix has computed every combination, and we pick the drawn one. The
    photons occupy the bits above the fermion legs, photon j in bit nl+j, so
    the drawn combination is one word. +1 -> eps+ -> 0.
  */
  size_t hg(0);
  for (size_t j(0); j < ng; ++j)
    if (m_PhoHel[j] <= 0) hg |= (((size_t)1) << j);
  /*
    The two codes do NOT label the helicity index of every leg the same way,
    and the mismatch is not a phase - it is a permutation of the array, so it
    survives squaring and cannot be argued away. It is visible directly in
    which entries are non-zero: the hand-coded Born gates on hel1 == -hel2,
    i.e. index h0 != h1, while Comix's non-zero entries have h0 == h1. One
    leg of each fermion pair is labelled oppositely.

    m_comixflip says which: bit i flips leg i (0,1 = beams; 2,3 = the
    outgoing pair, 4 = the photon). It is a SETTING, CEEX: COMIX_REAL_FLIP,
    and it was determined by MEASUREMENT in two steps, because no single
    measurement fixes all of it:

      - the four FERMION bits come from the soft probe, where the reference
        is the clean (eikonal x Born). Measured: mask 10 - the incoming
        POSITRON and the outgoing ANTI-fermion, i.e. both antiparticle legs
        - scores 3e-4 against 9e-2 for the complementary assignment. Three
        orders of magnitude; that is not a fit.
      - the PHOTON bit cannot come from there. The soft factor is a scalar in
        helicity space and |s_+| = |s_-|, so the wrong photon helicity is an
        overall constant in the soft limit and the probe is exactly
        degenerate in it (measured: identical metric to six digits). It comes
        instead from the hard-photon comparison against the hand-coded
        amplitude, which is decisive there: 0.046 with the bit set against
        0.84 without, on 563 of 570 events.

    Hence 26 = 11010b. Change it against the @@@ CEEXSOFTFLIP and @@@ CEEXFLIP
    lines, not by reasoning about conventions.
  */
  /*
    The fermion bits of the map are derived at the Born (CalibrateComixMap),
    where there is no photon to confuse them; the photon bit cannot come from
    there and is COMIX_REAL_PHOTON_FLIP. Photons are identical particles, so
    that one bit applies to every photon - there is nothing to tell them
    apart that could give them different conventions.
  */
  const size_t fmask((((size_t)1) << np) - 1);
  const size_t ffl((size_t)m_comixflip & fmask);
  const size_t gfl(m_comixphoflip ? ((((size_t)1) << ng) - 1) : 0);
  const size_t hgf(hg ^ gfl);
  const size_t nf(((size_t)1) << np);
  for (size_t f(0); f < nf; ++f)
    m_comixM1.m_A[f] = RealNorm() * sa[(f ^ ffl) | (hgf << np)];
  /*
    The raw table, kept so the flip scan can look at maps other than the one
    in force without a second Comix evaluation. It has one plane per photon
    helicity, so it only means anything for a single photon; above that the
    scan is not available and the diagnostics that read it say so.
  */
  if (ng == 1)
    for (int h = 0; h <= 1; ++h)
      for (size_t f(0); f < nf; ++f)
        m_comixraw[h].m_A[f] = RealNorm() * sa[f | ((size_t)h << np)];
  m_comixdrawnhel = (int)hg;

  m_havecomixreal = true;
  return true;
  #undef CXR_REFUSE
}


void Ceex_Base::ApplyComixReal()
{
  m_havecomixreal = false;
  if (!m_comixreal) return;
  /*
    Comix supplies the EXACT n-photon amplitude, which replaces the whole
    2^n partition sum rather than correcting it: the partition sum is CEEX's
    way of approximating that amplitude from soft factors times a reduced
    Born, so where the amplitude itself is available there is nothing left
    for the partitions to do. At n = 1 this is the beta_1 matching the file
    header describes; above it the same statement just carries more photons.
    Multiplicities with no real process built for them fall back to the
    hand-coded path, which is what FetchComixReal returning false does.
  */
  if (m_allphotons.empty()) return;
  /*
    The per-photon route replaces this one rather than joining it. Both write
    the real into m_AmpExpo1, so running them together would apply Comix twice
    by two different constructions.
  */
  if (m_perphoton) return;
  /*
    Above one photon this is a DIFFERENT MATCHING SCHEME, not more of the same
    one, and it is off by default until that is settled.

    The hand-coded routine it replaces is called InfraredSubtractedME_1_0, and
    the name is the argument: photon j's own soft factor is divided out
    (nrm = sProd/Sactu) and replaced by the exact emission amplitude, so the
    soft content already carried by exp(Y) and the crude S-factors is not
    counted twice. Comix's raw n-photon amplitude is not subtracted at all -
    it carries the full soft content of all n photons - so substituting it
    wholesale double counts the soft region. At n = 1 that is arranged for,
    because the Born partition terms and beta_1 are replaced together. At
    n > 1 nothing arranges it, and the exact n-photon amplitude additionally
    contains genuine double-hard beta_2 that a BETA: 1 run does not have.

    What the per-photon structure needs instead is Comix's ONE-photon
    amplitude evaluated once per photon at mapped momenta - the shape
    NLO_Base::MapMomenta already builds for YFS.NLO.
  */
  static const bool multi(ATOOLS::Settings::GetMainSettings()["CEEX"]
                          ["COMIX_REAL_MULTIPHOTON"].Get<int>() != 0);
  if (m_allphotons.size() > 1 && !multi) {
    CountComixReal(m_allphotons.size(), false);
    static bool warned(false);
    if (!warned) {
      warned = true;
      msg_Info()<<"CEEX: COMIX_REAL is applied at one photon only. Above that"
                <<" the exact n-photon amplitude is a different matching"
                <<" scheme, not this one - see ApplyComixReal. Set CEEX:"
                <<" COMIX_REAL_MULTIPHOTON to explore it anyway.\n";
    }
    return;
  }
  if (!m_comixcalibrated) DeriveComixMap();
  if (!m_comixreal) return;   // calibration may have switched it off
  const bool got(FetchComixReal());
  CountComixReal(m_allphotons.size(), got);
  if (!got) return;

  static const bool cxchk(ATOOLS::Settings::GetMainSettings()["CEEX"]["COMIX_CHECK"].Get<int>() != 0);
  const Amplitude hand1(m_AmpExpo1);   // Born + virtual + hand-coded real

  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f) {
          const Complex a0(m_AmpExpo0.m_A[f]);
          const Complex bv(m_AmpBornVirt.m_A[f]);
          // The virtual as the per-helicity multiplicative factor it is.
          // Both sides of this ratio are hand-coded, so every spinor-phase
          // convention cancels out of it; that is the whole point.
          const Complex V(std::abs(a0) > 0. ? bv/a0 : Complex(1., 0.));
          const Complex M1(m_comixM1.m_A[f]);
          m_AmpExpo1.m_A[f]    = V*M1;
          m_AmpBornReal.m_A[f] = M1;
          // The real INCREMENT over the Born partition sum, for the SHPART
          // diagnostics. Phase-sensitive by construction - it is the one
          // place the two conventions are subtracted from one another - and
          // it feeds nothing that is squared.
          m_snapReal.m_A[f]    = M1 - a0;
  }

  {
    // Norm-weighted, not worst-case: a worst case over 16 helicities is
    // always ~1, because the mass-suppressed configurations are ~0 in one
    // code and 1e-12 in the other however right everything is.
    double mnum(0.), mden(0.);
    double sh(0.), sc(0.), worst(0.);
    double sdiff(0.), sabsdiff(0.), s0(0.);
    for (int f = 0; f < nh; ++f) {
            const Complex H(hand1.m_A[f]);
            const Complex C(m_comixM1.m_A[f]);
            const Complex A0(m_AmpExpo0.m_A[f]);
            sh += std::norm(H); sc += std::norm(C);
            const double den(std::abs(H) + std::abs(C));
            if (den > 0.) worst = Max(worst, std::abs(H - C)/den);
            mnum     += sqr(std::abs(H) - std::abs(C));
            mden     += std::norm(H) + std::norm(C);
            s0       += std::norm(A0);
            sdiff    += std::norm(C - A0);
            sabsdiff += sqr(std::abs(C) - std::abs(A0));
    }
    const size_t ng(m_allphotons.size());
    if (m_cxrnormbyn.size() <= ng) { m_cxrnormbyn.resize(ng+1, 0.);
                                     m_cxrmetbyn.resize(ng+1, 0.); }
    if (sh > 0.) { m_cxrnormsum += sc/sh; ++m_cxrn; m_cxrnormbyn[ng] += sc/sh; }
    if (mden > 0.) { const double m(sqrt(mnum/mden));
                     m_cxrmetsum += m; m_cxrmetbyn[ng] += m; }
    if (!cxchk) return;
    SoftProbe();
    /*
      Which of the 2^(nlegs+1) index maps brings the two labellings into line.

      MAGNITUDES only, so the question "do the two codes mean the same
      helicity by this index" is answered without the separate question "do
      they use the same spinor phase" getting in the way. The metric is
      norm-weighted rather than a worst case: a worst case over 16 entries is
      always ~1, because an entry that is zero in one code and 1e-9 in the
      other has relative difference 1 however right the map is.
    */
    const int nlg(Amplitude::s_nlegs);   // helicity legs, not flavours
    const int fmaskx((1 << nlg) - 1);
    int bestmask(0); double bestmetric(1e30); double m0metric(-1.);
    for (int mask(0); mask < (1 << (nlg + 1)); ++mask) {
      double num(0.), den(0.);
      const int hgm(m_comixdrawnhel ^ ((mask >> nlg) & 1));
      for (int f = 0; f < nh; ++f) {
              const double H(std::abs(hand1.m_A[f]));
              const double C(std::abs(m_comixraw[hgm].m_A[f ^ (mask & fmaskx)]));
              num += sqr(H - C); den += H*H + C*C;
      }
      const double met(den > 0. ? sqrt(num/den) : 0.);
      if (mask == m_comixflip) m0metric = met;
      if (met < bestmetric) { bestmetric = met; bestmask = mask; }
    }
    std::cerr<<"@@@ CEEXFLIP inuse="<<m_comixflip<<" metric="<<m0metric
             <<" best="<<bestmask<<" bestmetric="<<bestmetric<<std::endl;
    const double ecm((m_momenta[0]+m_momenta[1]).Mass());
    const double x(ecm > 0. ? 2.*m_allphotons[0][0]/ecm : 0.);
    std::cerr<<"@@@ CEEXCMP nphot=1 x="<<x
             <<" sum_hand="<<sh<<" sum_comix="<<sc
             <<" ratio="<<(sc != 0. ? sh/sc : 0.)
             <<" worst_elem_reldiff="<<worst<<std::endl;
    std::cerr<<"@@@ CEEXSOFT x="<<x
             <<" rel="<<(s0 > 0. ? sqrt(sdiff/s0) : -1.)
             <<" relabs="<<(s0 > 0. ? sqrt(sabsdiff/s0) : -1.)
             <<" handrel="<<(s0 > 0. ? sqrt(std::abs(sh-s0)/s0) : -1.)
             <<std::endl;
    static int nrow(0);
    if (nrow++ < 3)
      for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
          for (int c = 0; c <= 1; ++c)
            for (int d = 0; d <= 1; ++d) {
              const Complex H(hand1.m_A[Idx(a,b,c,d)]);
              const Complex C(m_comixM1.m_A[Idx(a,b,c,d)]);
              std::cerr<<"@@@ CEEXHEL "<<a<<b<<c<<d
                       <<" absH="<<std::abs(H)<<" absC="<<std::abs(C)
                       <<" H=("<<H.real()<<","<<H.imag()<<")"
                       <<" C=("<<C.real()<<","<<C.imag()<<")"<<std::endl;
            }
  }
}


void Ceex_Base::CountComixReal(size_t ng, bool ok)
{
  std::vector<long> &v(ok ? m_cxrnbyn : m_cxrfailbyn);
  if (v.size() <= ng) v.resize(ng+1, 0);
  ++v[ng];
}


void Ceex_Base::ReportComixReal() const
{
  if (!m_comixreal) return;
  msg_Info()<<"CEEX: COMIX_REAL supplied the real amplitude on "<<m_cxrn
            <<" events";
  if (m_cxrn > 0)
    msg_Info()<<", mean |M1_comix|^2/|A1_hand|^2 = "<<(m_cxrnormsum/m_cxrn)
              <<", mean per-helicity magnitude mismatch "
              <<(m_cxrmetsum/m_cxrn);
  if (m_cxrfail)
    msg_Info()<<"; "<<m_cxrfail<<" events fell back to the hand-coded real";
  msg_Info()<<".\n";
  if (m_perphoton) {
    msg_Info()<<"CEEX: COMIX_REAL_PER_PHOTON replaced the one-photon"
              <<" amplitude on "<<m_cxppn<<" photons";
    if (m_cxppn > 0)
      msg_Info()<<", mean |C-H|/|C,H| = "<<(m_cxppdev/m_cxppn);
    if (m_cxppfail)
      msg_Info()<<"; "<<m_cxppfail<<" photons kept the hand-coded amplitude";
    msg_Info()<<".\n";
    return;
  }
  /*
    Broken down by photon multiplicity, because a run in which every event
    had one photon exercises none of the n-photon path and must not be read
    as evidence that it works. "taken" is where the Comix amplitude replaced
    the partition sum; "refused" is where FetchComixReal declined - no real
    process built at that multiplicity, momenta that do not balance, or more
    legs than the amplitude container holds - and the hand-coded path ran.
  */
  const size_t nmax(Max(m_cxrnbyn.size(), m_cxrfailbyn.size()));
  for (size_t n(1); n < nmax; ++n) {
    const long t(n < m_cxrnbyn.size() ? m_cxrnbyn[n] : 0);
    const long f(n < m_cxrfailbyn.size() ? m_cxrfailbyn[n] : 0);
    if (!t && !f) continue;
    msg_Info()<<"  "<<n<<" photon"<<(n>1?"s":"")<<": "<<t<<" taken, "
              <<f<<" refused";
    if (t > 0 && n < m_cxrnormbyn.size())
      msg_Info()<<", <|M_comix|^2/|A_hand|^2> = "<<(m_cxrnormbyn[n]/t)
                <<", <per-helicity mismatch> = "<<(m_cxrmetbyn[n]/t);
    msg_Info()<<"\n";
  }
}


/*!
  The soft limit, measured rather than argued.

  beta_1 must vanish as the photon softens: the exact one-photon amplitude
  goes to the eikonal times the Born,

      M_1(k) -> e^2 [ S_ini(k) + S_fin(k) ] M_0 ,

  with S the very soft factors CEEX already builds in CalculateSfactors. So
  the ratio

      R(lambda) = sum_hel |M1_comix(lambda k)|^2
                / sum_hel |e^2 (S_ini + S_fin) Born_hand|^2

  must go to 1 as lambda -> 0, and it does so regardless of process, spinor
  phase convention or gauge reference, because both sides are squared. What R
  tests is everything that is NOT a phase: the overall normalisation of
  Comix's amplitudes against CEEX's, the helicity INDEX map, and the choice
  of photon helicity - get any of those wrong and R tends to something else,
  or to nothing at all.

  It has to be a probe rather than a check on the event as generated, because
  the events do not supply the limit: a one-photon CEEX event is one whose
  single RESOLVED photon is hard, and every one of them measured here had
  x = 2E/sqrt(s) > 0.1. So the photon is scaled down by hand, the outgoing
  pair rebuilt to keep momentum conserved and both masses on shell, and the
  ladder in lambda printed.

  Costs one Comix evaluation per rung, and perturbs the random sequence
  through GeneratePoint(), so it runs only under SHERPA_CEEX_COMIX and only
  for the first few events.
*/
void Ceex_Base::SoftProbe()
{
  static int ndone(0);
  static const int nmax(ATOOLS::Settings::GetMainSettings()["CEEX"]
                        ["SOFT_PROBE_EVENTS"].Get<int>());
  if (ndone >= nmax) return;
  if (!m_havecomixreal || m_allphotons.size() != 1) return;
  YFS::Real_Correction *p_realprov(RealProvider(1));
  if (p_realprov == NULL) return;
  const Vec4D p1(m_pceex[0]), p2(m_pceex[1]), k0(m_allphotons[0]);
  const Vec4D Q0(m_pceex[m_if1] + m_pceex[m_if2]);
  const double m3(m_flavs[m_if1].Mass()), m4(m_flavs[m_if2].Mass());
  Vec4D d3(m_pceex[m_if1]);
  { Poincare c0(Q0); c0.Boost(d3); }
  const double n3(d3.PSpat());
  if (n3 <= 0.) return;
  ++ndone;

  const double sp_save(m_sp);
  const Complex qratio(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
  const int hg(m_PhoHel[0]);

  for (int it(0); it <= 6; ++it) {
    const double lam(pow(10., -double(it)));
    const Vec4D k(lam*k0);
    const Vec4D Q(p1 + p2 - k);
    const double sQ(Q.Abs2());
    if (sQ <= sqr(m3 + m4)) continue;
    const double rs(sqrt(sQ));
    const double E3((sQ + m3*m3 - m4*m4)/(2.*rs));
    const double pm(sqrt(Max(0., E3*E3 - m3*m3)));
    Vec4D q3(E3, pm*d3[1]/n3, pm*d3[2]/n3, pm*d3[3]/n3);
    Vec4D q4(rs - E3, -q3[1], -q3[2], -q3[3]);
    Poincare cms(Q);
    cms.BoostBack(q3); cms.BoostBack(q4);

    Vec4D_Vector pp(5);
    pp[0]=p1; pp[1]=p2; pp[2]=q3; pp[3]=q4; pp[4]=k;
    const std::vector<METOOLS::Spin_Amplitudes> *amps
      (p_realprov->ComixAmplitudes(pp));
    if (amps == NULL || amps->empty() || (*amps)[0].size() != 32) continue;
    const METOOLS::Spin_Amplitudes &sa((*amps)[0]);

    // hand side: the eikonal times the Born, at THIS configuration
    Vec4D_Vector bp(4);
    bp[0]=p1; bp[1]=p2; bp[2]=q3; bp[3]=q4;
    m_sp = sQ;
    MakeProp();
    MakePropT(bp);
    Amplitude B;
    BornAmplitude(bp, B, -1., -1., -1);
    const Complex Si(Sfactor(p1, p2, k, hg));
    const Complex Sf(qratio*Sfactor(q3, q4, k, hg));
    const Complex soft(m_e*m_e*(Si + Sf));

    const int fl(m_comixflip);
    const int nlg(Amplitude::s_nlegs);   // helicity legs, not flavours
    const int fmaskx((1 << nlg) - 1);
    const int nh(Amplitude::NHel());
    const int hgi(((hg > 0 ? 0 : 1) ^ ((fl >> nlg) & 1)));
    double sc(0.), shd(0.), scall(0.);
    for (int f = 0; f < nh; ++f) {
      const Complex C(sa[FlatIndex((size_t)(f ^ (fl & fmaskx)), hgi, nlg)]);
      sc  += std::norm(RealNorm()*C);
      shd += std::norm(soft*B.m_A[f]);
    }
    for (size_t i(0); i < sa.size(); ++i) scall += std::norm(sa[i]);

    /*
      The same soft statement tested WITH PHASE, and entirely inside Comix's
      convention.

      Everything above compares magnitudes, and says so. That is enough to fix
      the index map and the normalisation, and not enough to subtract: to form
      the infrared-subtracted remainder M_1 - S B - which is the object CEEX's
      per-photon routine actually holds - the eikonal term has to cancel
      against M_1 as a COMPLEX number, per helicity. Subtracting a hand-coded
      Born from a Comix amplitude was measured to make |M_1| grow rather than
      shrink, because the two carry a per-helicity phase that is not constant.

      So both terms are Comix's own here. The Born needs its own 2 -> 2 point:
      the fermion legs of the 2 -> 3 configuration do NOT conserve momentum on
      their own (P1 + P2 = q3 + q4 + k), and handing that to Comix produced
      NaN. The beams are therefore rebuilt at (q3+q4)^2.

      m_comixnorm cancels between numerator and denominator, so this test is
      independent of the normalisation it would be used with - it isolates the
      phase alone.
    */
    {
      const Vec4D Q2(q3 + q4);
      const double s2(Q2.Abs2());
      const double mm1(m_flavs[0].Mass()), mm2(m_flavs[1].Mass());
      const double lam2(sqr(s2) + sqr(mm1*mm1) + sqr(mm2*mm2)
                        - 2.*s2*mm1*mm1 - 2.*s2*mm2*mm2 - 2.*mm1*mm1*mm2*mm2);
      if (s2 > 0. && lam2 > 0.) {
        Vec4D b3(q3), b4(q4);
        Poincare bq(Q2); bq.Boost(b3); bq.Boost(b4);
        const double lcm(0.5*sqrt(lam2/s2));
        const double sgn(m_bornmomenta[0][3] < 0 ? -1. : 1.);
        Vec4D B1(lcm*sqrt(1.+mm1*mm1/sqr(lcm)), 0., 0.,  sgn*lcm);
        Vec4D B2(lcm*sqrt(1.+mm2*mm2/sqr(lcm)), 0., 0., -sgn*lcm);
        Poincare bq2(Q2); bq2.BoostBack(B1); bq2.BoostBack(B2);
        Vec4D_Vector b2{B1, B2, q3, q4};
        Amplitude Bc;
        if (ComixBornAmplitude(b2, Bc)) {
          double amin(-1.), amax(-1.), phmin(0.), phmax(0.);
          int nl2(0);
          double nc2(0.);
          for (int f = 0; f < nh; ++f)
            nc2 += std::norm(sa[FlatIndex((size_t)(f ^ (fl & fmaskx)),
                                          hgi, nlg)]);
          nc2 = sqrt(nc2);
          for (int f = 0; f < nh; ++f) {
            const Complex Cf(sa[FlatIndex((size_t)(f ^ (fl & fmaskx)),
                                          hgi, nlg)]);
            const Complex Df(soft*Bc.m_A[f ^ (fl & fmaskx)]);
            if (std::abs(Cf) < 1e-4*nc2 || std::abs(Df) == 0.) continue;
            const Complex r(Cf/Df);
            const double a(std::abs(r)), ph(std::arg(r));
            if (nl2++ == 0) { amin = amax = a; phmin = phmax = ph; }
            else { amin = Min(amin,a); amax = Max(amax,a);
                   phmin = Min(phmin,ph); phmax = Max(phmax,ph); }
          }
          /*
            Carried alongside so the events whose phase is NOT 0 or pi can be
            told apart from the ones whose is. A relative SIGN between the
            two processes would put every event at one or the other; anything
            else has to correlate with something physical.
          */
          const double ct(q3.PSpat() > 0. && k.PSpat() > 0.
                          ? (Vec3D(q3)*Vec3D(k))/(q3.PSpat()*k.PSpat()) : 0.);
          const double cb(k.PSpat() > 0. ? k[3]/k.PSpat() : 0.);
          std::cerr<<"@@@ CEEXSOFTPHASE lam="<<lam<<" nlive="<<nl2
                   <<" |r|=["<<amin<<","<<amax<<"]"
                   <<" arg(r)=["<<phmin<<","<<phmax<<"]"
                   <<" hel="<<hg<<" cos_kf="<<ct<<" cos_kbeam="<<cb
                   <<" phi_k="<<atan2(k[2], k[1])
                   <<" sQ="<<s2<<std::endl;
        } else {
          std::cerr<<"@@@ CEEXSOFTPHASE lam="<<lam
                   <<" Comix Born unavailable"<<std::endl;
        }
      }
    }
    /*
      The index map, determined where it can be determined: in the soft
      limit the hand side is (eikonal x Born), a clean reference with the
      same helicity structure as M_1, so the map that lines the two up is
      the right one rather than the least wrong one. Magnitudes only.
    */
    int bm(0); double bmet(1e30), inusemet(-1.), phometh(-1.);
    for (int mask(0); mask < (1 << (nlg + 1)); ++mask) {
      double num(0.), den(0.);
      const int hgm(((hg > 0 ? 0 : 1) ^ ((mask >> nlg) & 1)));
      for (int f = 0; f < nh; ++f) {
        const double C(std::abs(RealNorm()*
                                sa[FlatIndex((size_t)(f ^ (mask & fmaskx)),
                                             hgm, nlg)]));
        const double H(std::abs(soft*B.m_A[f]));
        num += sqr(H - C); den += H*H + C*C;
      }
      const double met(den > 0. ? sqrt(num/den) : 0.);
      if (mask == m_comixflip) inusemet = met;
      if (mask == (m_comixflip ^ 16)) phometh = met;
      if (met < bmet) { bmet = met; bm = mask; }
    }
    std::cerr<<"@@@ CEEXSOFTFLIP lam="<<lam<<" inuse="<<m_comixflip
             <<" metric="<<inusemet<<" photonflipped="<<phometh
             <<" hel="<<hg
             <<" best="<<bm<<" bestmetric="<<bmet<<std::endl;
    /*
      NLO_Base::CheckRealCollinearSub's test, applied to Comix's real: divide
      the raw squared ME by the CLOSED-FORM massive eikonal of the two charged
      initial legs. Pure kinematics - no spinor convention and no CEEX soft
      factor, so it cannot fail in the same way as the thing being tested. The
      eikonal theorem forces this to be FLAT along the ladder; flat is the
      verdict, not convergence to any particular value.

      It has to be a scan at FIXED kinematics. Taken across events instead,
      the ratio spans eight orders - but that is the Born and the Z propagator
      moving between events, not the eikonal.
    */
    { const double mel(m_flavs[0].Mass());
      const double pk1(p1*k), pk2(p2*k), p12(p1*p2);
      const double Scl(pk1 != 0. && pk2 != 0.
                       ? 2.*p12/(pk1*pk2) - mel*mel/(pk1*pk1)
                                          - mel*mel/(pk2*pk2) : 0.);
      double nb(0.);
      for (int f = 0; f < nh; ++f) nb += std::norm(B.m_A[f]);
      std::cerr<<"@@@ CEEXEIKFLAT lam="<<lam
               <<" x="<<(2.*k[0]/sqrt(m_s))
               <<" ratio="<<((Scl != 0. && nb > 0.) ? sc/(Scl*nb) : -1.)
               <<std::endl; }
    std::cerr<<"@@@ CEEXSOFTSCAN lam="<<lam
             <<" x="<<2.*k[0]/sqrt(dabs((p1+p2).Abs2()))
             <<" R="<<(shd > 0. ? sc/shd : -1.)
             <<" sum_comix="<<sc<<" sum_eikborn="<<shd
             <<" allhel/diff="<<(p_realprov->m_lastcomix != 0. ?
                                 scall/p_realprov->m_lastcomix : 0.)
             <<std::endl;
  }
  m_sp = sp_save;
  MakeProp();
  MakePropT(m_pceex);
}


bool Ceex_Base::CalibrateComixMap(ComixCalib &c)
{
  if (p_bornproc == NULL || m_bornmomenta.size() < 4) return false;
  // The Born configuration, in its own rest frame.
  Vec4D_Vector bp(m_bornmomenta);
  Poincare bcms(bp[0] + bp[1]);
  for (size_t i(0); i < bp.size(); ++i) bcms.Boost(bp[i]);

  Amplitude cx, hand;
  const double sp_save(m_sp);
  m_sp = (bp[2] + bp[3]).Abs2();
  MakeProp();
  BornAmplitude(bp, hand);
  double cxme2(0.);
  const bool ok(ComixBornAmplitude(bp, cx, &cxme2));
  m_sp = sp_save;
  MakeProp();
  if (!ok) return false;

  const int nh(Amplitude::NHel());
  double sh(0.), sc(0.);
  for (int f = 0; f < nh; ++f) {
    sh += std::norm(m_e*m_e*hand.m_A[f]);
    sc += std::norm(cx.m_A[f]);
  }
  if (!(sh > 0.) || !(sc > 0.)) return false;

  const double rh(1./sqrt(sh)), rc(1./sqrt(sc));
  int bestmask(-1);
  double bestmet(-1.), nextmet(-1.);
  for (int m = 0; m < nh; ++m) {
    double met(0.);
    for (int f = 0; f < nh; ++f)
      met += sqr(std::abs(m_e*m_e*hand.m_A[f])*rh - std::abs(cx.m_A[f ^ m])*rc);
    met = sqrt(met);
    if (bestmask < 0 || met < bestmet) { nextmet = bestmet;
                                         bestmet = met; bestmask = m; }
    else if (nextmet < 0. || met < nextmet) nextmet = met;
  }

  /*
    With the winning mask in hand, |C/H| has to be the SAME number in every
    helicity slot that carries weight, or the two objects are not the same
    amplitude in two conventions. The threshold below is deliberately coarse:
    with a massive electron Comix populates the helicity-FLIP entries at the
    m_e/E level, which the hand-coded Born sets to exactly zero by its
    hel1 == -hel2 gate, so a loose cut mixes a convention test with a
    mass-suppression test and the ratio wanders for a reason that has nothing
    to do with conventions.
  */
  double rmin(0.), rmax(0.);
  int nlive(0);
  for (int f = 0; f < nh; ++f) {
    const Complex H(m_e*m_e*hand.m_A[f]), C(cx.m_A[f ^ bestmask]);
    if (std::abs(H) > 1e-3*sqrt(sh) && std::abs(C) > 1e-3*sqrt(sc)) {
      const double a(std::abs(C/H));
      if (nlive++ == 0) rmin = rmax = a;
      else { rmin = Min(rmin, a); rmax = Max(rmax, a); }
    }
  }

  c.mask = bestmask; c.met = bestmet; c.next = nextmet;
  c.N = sqrt(sc/sh);  c.sh = sh; c.sc = sc; c.me2 = cxme2;
  c.rmin = rmin; c.rmax = rmax; c.nlive = nlive;
  return true;
}


/*!
  Adopt the calibration, for whichever of the two constants was left to be
  derived (COMIX_REAL_FLIP < 0, COMIX_REAL_NORM <= 0). One attempt per run:
  the map is a property of the two codes' conventions, not of the event.
*/
void Ceex_Base::DeriveComixMap()
{
  m_comixcalibrated = true;
  const bool needflip(m_comixflip < 0), neednorm(!(m_comixnorm > 0.));
  // a normalisation set by hand in the card is taken as given
  if (!neednorm) m_normexact = true;
  if (!needflip && !neednorm) return;

  ComixCalib c;
  if (!CalibrateComixMap(c)) {
    // Falling back to the fitted values would hide the failure behind numbers
    // that happen to be right for THIS process, so refuse instead.
    msg_Error()<<METHOD<<"(): Born calibration failed; the Comix real cannot"
               <<" be normalised. Set CEEX: COMIX_REAL_NORM and"
               <<" COMIX_REAL_FLIP by hand, or leave COMIX_REAL off.\n";
    m_comixreal = 0;
    return;
  }

  const size_t nl(m_flavs.size());
  /*
    The photon bit sits above the FERMION helicity slots, and there are
    Amplitude::s_nlegs of those - legs with two spin states. m_flavs.size()
    counts flavours, which is the same number only when every leg is a
    fermion: for H l+ l- it is 5 against 4, and the bit landed at 32 instead
    of 16, outside the helicity range entirely. Invisible at 2 -> 2.
  */
  if (needflip) m_comixflip = c.mask | (m_comixphoflip ? Amplitude::NHel() : 0);

  /*
    N = 2 is not a coincidence and not a fit. sum_hel |A_comix|^2 comes out at
    four times Comix's own matrix element, and Comix's matrix element is
    spin-AVERAGED; sum_hel |e^2 A_hand|^2 reproduces that averaged element
    directly. So the hand-coded CEEX amplitude carries sqrt(1/4) - the square
    root of the initial state's spin multiplicity - inside its normalisation,
    which is KKMC's convention.

    That makes the value EXACT, and the measurement its confirmation rather
    than its source: the constant used is 1/sqrt(2s+1 per incoming leg), and
    the measured N only has to agree. Taking the measured number instead would
    make the whole run depend on which event happened to calibrate it, at the
    1e-11 level where the two Born constructions differ over the electron
    mass - reproducibility thrown away for nothing.
  */
  double nspin(1.);
  for (size_t i(0); i < 2 && i < nl; ++i) nspin *= m_flavs[i].IntSpin() + 1.;
  const double expected(nspin > 0. ? 1./sqrt(nspin) : 0.);
  if (neednorm) {
    const bool agrees(expected > 0. && c.N > 0.
                      && std::abs(c.N*expected - 1.) < 1e-2);
    m_comixnorm = agrees ? expected : (c.N > 0. ? 1./c.N : 0.);
    m_normexact = agrees;
    if (!agrees)
      msg_Error()<<METHOD<<"(): measured N = "<<c.N<<" is not"
                 <<" sqrt(initial spin states) = "<<(expected>0.?1./expected:0.)
                 <<"; falling back to the measured value, but the amplitudes"
                 <<" are then not in the convention assumed here.\n";
  }

  msg_Info()<<METHOD<<"(): Comix -> CEEX map derived from the Born:\n"
            <<"  flip mask   = "<<m_comixflip<<"  (fermion bits "<<c.mask
            <<", metric "<<c.met<<" against "<<c.next<<" for the runner-up"
            <<(m_comixphoflip ? "; photon bit set" : "; photon bit clear")<<")\n"
            <<"  normalisation = "<<m_comixnorm<<"  (N = "<<c.N
            <<"; 1/sqrt(initial spin states) = "<<expected<<")\n"
            <<"  |C/H| over the "<<c.nlive<<" live helicities: ["
            <<c.rmin<<", "<<c.rmax<<"]\n";
}


/*!
  The hand-coded one-photon amplitude at a given configuration, ISR + FSR,
  returned rather than accumulated.

  This is the same spinor construction InfraredSubtractedME_1_0 and its FSR
  partner build, with the partition weights (sProd/Sactu) and the kinematic
  (1 - CKine) term left out: those belong to the partition, not to the photon,
  and what is wanted here is the photon's own amplitude so that it can be
  compared with Comix's. Summing the two stages is what makes it the FULL
  one-photon amplitude rather than one stage's share of it.
*/
bool Ceex_Base::HandOnePhotonAmplitude(const Vec4D_Vector &bp, const Vec4D &k,
                                       int hel, Amplitude &M1)
{
  if (bp.size() < 4) return false;
  const double m3(m_flavs[2].Mass()), m4(m_flavs[3].Mass());
  const double p1k(k*bp[0]), p2k(k*bp[1]), p3k(k*bp[2]), p4k(k*bp[3]);
  if (IsZero(p1k) || IsZero(p2k) || IsZero(p3k) || IsZero(p4k)) return false;
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f) M1.m_A[f] = Complex(0., 0.);

  // --- initial state: the beam spinors carry the photon index in turn
  {
    const Vec4D_Vector arm1{k, bp[1], bp[2], bp[3]};
    const Vec4D_Vector arm2{bp[0], k, bp[2], bp[3]};
    Amplitude AmpBornU, AmpBornV, AmpU, AmpV;
    BornAmplitude(arm1, AmpBornU, -1., -1., -1);
    BornAmplitude(arm2, AmpBornV, -1., -1., -1);
    const double mI(m_flavs[0].Mass());
    UGamma(k, bp[0], k, hel, AmpU, 0., mI);
    VGamma(bp[1], k, k, hel, AmpV, mI, 0.);
    const double gI(m_qe * m_e * m_e * m_e);
    AddU(M1, AmpBornU, AmpU,  gI / p1k / 2.);
    AddV(M1, AmpBornV, AmpV, -gI / p2k / 2.);
  }
  // --- final state: the same, on the outgoing pair
  {
    const Vec4D_Vector arm1{bp[0], bp[1], k, bp[3]};
    const Vec4D_Vector arm2{bp[0], bp[1], bp[2], k};
    Amplitude AmpBornU, AmpBornV, AmpU, AmpV;
    BornAmplitude(arm1, AmpBornU, 0., m4, -1);
    BornAmplitude(arm2, AmpBornV, m3, 0., -1);
    UGamma(bp[2], k, k, hel, AmpU, m3, 0.);
    VGamma(k, bp[3], k, hel, AmpV, 0., m4);
    const double gF(m_qf * m_e * m_e * m_e);
    AddUF(M1, AmpU, AmpBornU,  gF / p3k / 2.);
    AddVF(M1, AmpBornV, AmpV, -gF / p4k / 2.);
  }
  return true;
}


/*!
  Comix's one-photon amplitude for ONE photon of the event, with the others
  absorbed into the recoil of the outgoing pair.

  Comix needs a momentum-conserving on-shell point and the hand-coded spinor
  algebra does not: it is handed spinor ARGUMENTS, not a momentum set, which
  is why it can evaluate at m_pceex plus one photon while ignoring the others.
  So a configuration has to be built. The one used here is the soft probe's,
  which is already verified to reproduce the eikonal limit: keep the beams,
  take Q = p1 + p2 - k, and rebuild the pair on shell in Q's frame along the
  physical pair's direction. The other photons then live in the recoil rather
  than as legs, and momentum conservation is exact by construction.

  \param bp receives the reduced fermion configuration, so that the hand-coded
            amplitude can be built at the SAME point - a ratio between
            amplitudes at different momenta would mean nothing.
*/
/*!
  Comix's one-photon amplitude at a momentum set supplied verbatim, masked and
  normalised into CEEX's convention. No configuration is built here, so the
  caller can hold the fermions fixed and move only the photon - which is what
  the beta_1 subtraction needs, and which requires
  COMIX: MOMENTUM_PROJECTION: 0 because such a set does not conserve momentum.
*/
bool Ceex_Base::ComixRealAt(const Vec4D_Vector &pp, int hel, Amplitude &M1,
                            const double propscale)
{
  YFS::Real_Correction *prov(RealProvider(1));
  if (prov == NULL) return false;
  /*
    nlg indexes HELICITY, so it counts legs with two spin states, not
    flavours. The two agree only when every leg is a fermion; for H l+ l- they
    are 4 and 5, and the flavour count put the photon bit and the whole
    FlatIndex stride one power of two too high. Invisible at 2 -> 2.
  */
  const int nlg(Amplitude::s_nlegs);
  if (pp.size() != m_flavs.size() + 1) return false;
  const std::vector<METOOLS::Spin_Amplitudes> *amps
    (prov->ComixAmplitudes(pp, propscale));
  if (amps == NULL || amps->empty()) return false;
  const METOOLS::Spin_Amplitudes &sa((*amps)[0]);
  if ((int)sa.size() != (1 << (nlg + 1))) return false;
  const int fl(m_comixflip), fmaskx((1 << nlg) - 1), nh(Amplitude::NHel());
  const int hgi((hel > 0 ? 0 : 1) ^ ((fl >> nlg) & 1));
  for (int f = 0; f < nh; ++f)
    M1.m_A[f] = RealNorm()
      * sa[FlatIndex((size_t)(f ^ (fl & fmaskx)), hgi, nlg)];
  return true;
}

bool Ceex_Base::ComixOnePhotonAmplitude(const Vec4D &k, int hel,
                                        Vec4D_Vector &bp, Amplitude &M1,
                                        Vec4D &kmap)
{
  YFS::Real_Correction *prov(RealProvider(1));
  if (prov == NULL || prov->p_proc == NULL || m_pceex.size() < 4) return false;
  const int nlg((int)m_flavs.size());
  if (nlg != 4) return false;   // the reduction below is the 2 -> 2 one

  /*
    Take Q from the OUTGOING side. Building it as p1 + p2 - k instead leaves
    the other photons' momentum inside the pair, which then carries energy
    that belongs to them: measured, the pair came out at the full s on every
    soft-photon configuration, and the hand-coded amplitude collapsed by
    orders of magnitude against Comix's while the two agreed to ~10% on hard
    photons. It is the mirror of the Real_Map fault - that one subtracted the
    photon twice, this one did not subtract the others at all.

    So the reduced system is exactly this photon and the pair, and it is the
    BEAMS that are rebuilt to match it. That is NLO_Base::MapMomenta's
    construction, which is also what YFS.NLO feeds its real ME, so the two
    agree on the kinematics instead of each inventing one.
  */
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  const double m3(m_flavs[2].Mass()), m4(m_flavs[3].Mass());
  Vec4D q3(m_pceex[2]), q4(m_pceex[3]), kk(k);
  const Vec4D Q(q3 + q4 + kk);
  const double sQ(Q.Abs2());
  if (sQ <= sqr(m3 + m4)) return false;
  // Kaellen function, written out rather than pulled in from Dipole.H
  const double lam(sqr(sQ) + sqr(m1*m1) + sqr(m2*m2)
                   - 2.*sQ*m1*m1 - 2.*sQ*m2*m2 - 2.*m1*m1*m2*m2);
  if (lam <= 0.) return false;

  Poincare boostLab(m_bornmomenta[0] + m_bornmomenta[1]);
  Poincare pRot(m_bornmomenta[0], Vec4D(0., 0., 0., 1.));
  Poincare boostQ(Q);
  pRot.RotateBack(q3); pRot.RotateBack(q4); pRot.RotateBack(kk);
  boostQ.Boost(q3);    boostQ.Boost(q4);    boostQ.Boost(kk);

  const double sign_z(m_bornmomenta[0][3] < 0 ? -1. : 1.);
  const double lamCM(0.5*sqrt(lam/sQ));
  const double E1(lamCM*sqrt(1. + m1*m1/sqr(lamCM)));
  const double E2(lamCM*sqrt(1. + m2*m2/sqr(lamCM)));
  Vec4D P1(E1, 0., 0.,  sign_z*lamCM);
  Vec4D P2(E2, 0., 0., -sign_z*lamCM);

  Poincare pRot2(m_bornmomenta[0], Vec4D(0., 0., 0., 1.));
  Vec4D *legs[5] = {&P1, &P2, &q3, &q4, &kk};
  for (int i(0); i < 5; ++i) { pRot2.Rotate(*legs[i]);
                               boostLab.BoostBack(*legs[i]); }

  // Comix at a non-conserving point returns a number, not an error.
  { Vec4D bal(P1 + P2 - q3 - q4 - kk);
    double worst(0.);
    for (int c(0); c < 4; ++c) worst = Max(worst, dabs(bal[c]));
    if (worst > 1e-6*(sqrt(dabs(sQ)) + 1.)) return false; }

  bp.assign({P1, P2, q3, q4});
  Vec4D_Vector pp{P1, P2, q3, q4, kk};
  const std::vector<METOOLS::Spin_Amplitudes> *amps(prov->ComixAmplitudes(pp));
  if (amps == NULL || amps->empty()) return false;
  const METOOLS::Spin_Amplitudes &sa((*amps)[0]);
  if ((int)sa.size() != (1 << (nlg + 1))) return false;

  const int fl(m_comixflip), fmaskx((1 << nlg) - 1), nh(Amplitude::NHel());
  const int hgi((hel > 0 ? 0 : 1) ^ ((fl >> nlg) & 1));
  for (int f = 0; f < nh; ++f)
    M1.m_A[f] = RealNorm()
      * sa[FlatIndex((size_t)(f ^ (fl & fmaskx)), hgi, nlg)];
  kmap = kk;
  return true;
}


/*!
  One ratio per photon, built before the partition loop.

  Applied inside the loop it turns every partition's share of that photon into
  Comix's amplitude, because the ratio is constant in the partition index and
  the shares add: sum_s w_s B_s r = r sum_s w_s B_s, and sum_s w_s B_s is the
  hand-coded amplitude the ratio divides by.

  The approximation this makes, stated plainly: r is one number per helicity,
  so the SAME correction is applied to the photon's initial-state and
  final-state shares. It would be exact if Comix could be asked for the
  amplitude with only one stage's emitters radiating, which it cannot. The
  two shares agree in the soft limit, where both reduce to the eikonal, so
  the approximation is in the hard region - the same region where the
  hand-coded and Comix amplitudes differ by the 9% that is not yet explained.
*/
void Ceex_Base::BuildComixPhotonRatios()
{
  const size_t ng(m_allphotons.size());
  m_cxratio.assign(ng, Amplitude());
  m_cxcsub.assign(ng, Amplitude());
  m_cxratiook.assign(ng, 0);
  if (!m_comixreal || !m_perphoton || ng == 0) return;
  if (!m_comixcalibrated) DeriveComixMap();
  if (!m_comixreal) return;

  const int nh(Amplitude::NHel());
  const double sp_save(m_sp);
  /*
    The scale factor is chosen per photon, not fixed. What matters is that the
    scaled photon sits deep in the EIKONAL regime - x = 2E/sqrt(s) well below
    one - and pushing it far past that only costs precision: the spinor
    products of a 1e-11 GeV photon against 125 GeV beams retain about three
    digits, and the cancellation needed to expose beta_1 wants seven. Measured
    with a fixed lambda = 1e-5, soft photons subtracted only 13% instead of
    the six orders required.

    A photon already below X_EIK is its own eikonal: beta_1 there is smaller
    than the hand-coded value by the same ratio it would be corrected by, so
    the hand-coded amplitude is kept rather than a difference of two nearly
    equal numbers being trusted.
  */
  static const double X_EIK(1e-4);

  for (size_t j(0); j < ng; ++j) {
    const Vec4D k0(m_allphotons[j]);
    const int hel(m_PhoHel[j]);

    /*
      Comix's amplitude at this photon, and at the SAME photon scaled down.
      Both are Comix evaluations, and Amplitude::SetGauge uses a fixed global
      reference vector rather than one built from the momenta, so the two are
      phase comparable to each other. That is what makes the subtraction below
      well formed where subtracting a Born was not.
    */
    Vec4D_Vector bp, bps;
    Vec4D kmap, kmaps;
    Amplitude C, Cs;
    /*
      CEEX's OWN arguments, unmapped. Its beta_1 is built from spinor
      arguments - arm1 = {k, p1, p2, p3}, the physical momenta with a beam
      spinor replaced by the photon - which do not conserve momentum once
      other photons are present. That is the same situation as the partition
      Born, where handing Comix exactly those arguments reproduced CEEX's
      B(X) to 1e-8 once MOMENTUM_PROJECTION was off.

      Mapping to a conserving configuration was therefore never necessary,
      and it was harmful: rebuilding the pair moves the Born by O(E_gamma),
      the same order as beta_1 itself.
    */
    static const bool ceexargs(ATOOLS::Settings::GetMainSettings()["CEEX"]
                               ["COMIX_REAL_CEEX_ARGS"].Get<int>() != 0);
    if (ceexargs) {
      bp.assign(m_pceex.begin(), m_pceex.begin()+4);
      kmap = k0;
      Vec4D_Vector ppc{bp[0], bp[1], bp[2], bp[3], kmap};
      if (!ComixRealAt(ppc, hel, C, m_svarQ)) { ++m_cxppfail; continue; }
    } else
    if (!ComixOnePhotonAmplitude(k0, hel, bp, C, kmap)) { ++m_cxppfail; continue; }
    /*
      The SAME fermion configuration, with only the photon scaled. Rebuilding
      the pair for the scaled photon - which is what this did - moves the Born
      by O(E_gamma), and that is the order of beta_1 itself, so the
      subtraction carried an O(1) relative error on the very thing it was
      extracting. The Born must stay put while the eikonal grows.

      The set no longer conserves momentum, by (1-lambda) k. Comix accepts
      that with MOMENTUM_PROJECTION off, which is the whole reason this is
      now expressible.
    */
    const double xg(2.*kmap[0]/sqrt(m_s));
    if (!(xg > 10.*X_EIK)) { ++m_cxppfail; continue; }  // already eikonal
    const double lam(X_EIK/xg);
    kmaps = lam*kmap;
    bps   = bp;
    Vec4D_Vector pps{bp[0], bp[1], bp[2], bp[3], kmaps};
    if (!ComixRealAt(pps, hel, Cs, ceexargs ? m_svarQ : -1.))
      { ++m_cxppfail; continue; }

    /*
      The eikonal RATIO, not the eikonal. CEEX's soft factor carries the
      photon polarisation phase, which is fixed by the photon's DIRECTION -
      measured: |r| = 1/e^2 in every angular band while arg(r) tracks the
      angle, and arg flips as r_- = -conj(r_+) under a helicity pin, which is
      the eps_- = -eps_+* relation. Scaling k does not rotate it, so that
      phase is identical in numerator and denominator and cancels exactly.
      Nothing here needs it to be known.
    */
    const Complex qratio(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
    const Complex Sk (Sfactor(bp[0],  bp[1],  kmap,  hel)
                      + qratio*Sfactor(bp[2],  bp[3],  kmap,  hel));
    const Complex Sks(Sfactor(bps[0], bps[1], kmaps, hel)
                      + qratio*Sfactor(bps[2], bps[3], kmaps, hel));
    if (std::abs(Sks) == 0.) { ++m_cxppfail; continue; }
    const Complex srat(Sk/Sks);

    Amplitude Csub;   // filled below, once Comix's own Born is in hand

    /*
      Bring it into CEEX's convention. The Born fixes this completely: at the
      same 2 -> 2 point, a_f = e^2 B_hand(f) / B_comix(f) is the per-helicity
      generalisation of COMIX_REAL_NORM - its modulus IS that constant, 0.5,
      and its phase is the fermion-spinor convention difference, which the
      Born calibration showed is NOT constant across helicities and so cannot
      be a single number.

      The Born needs its own 2 -> 2 point: the fermion legs of the 2 -> 3
      configuration do not conserve momentum by themselves, and handing that
      to Comix returns NaN.
    */
    Amplitude Bc, Bh, H;
    bool okb(false), okh(false);
    if (ceexargs) {
      /*
        The Born at the SAME arguments and the same scale, so both terms of
        the subtraction carry identical Born content. No 2 -> 2 point is
        constructed: CEEX's arguments are used verbatim, non-conserving, which
        is the whole point - the previous version built a conserving
        configuration and thereby moved the Born by O(E_gamma).
      */
      m_sp = m_svarQ;
      MakeProp();
      MakePropT(bp);
      okb = ComixBornAmplitude(bp, Bc, NULL, m_svarQ);
      BornAmplitude(bp, Bh, -1., -1., -1);
      okh = HandOnePhotonAmplitude(bp, kmap, hel, H);
      m_sp = sp_save;
      MakeProp();
      MakePropT(m_pceex);
    } else {
      /*
        The mapped route: the Born needs its own 2 -> 2 point, because the
        fermion legs of the 2 -> 3 configuration do not conserve momentum by
        themselves and handing that to Comix returns NaN.
      */
      const Vec4D Q2(bp[2] + bp[3]);
      const double s2(Q2.Abs2());
      const double mm1(m_flavs[0].Mass()), mm2(m_flavs[1].Mass());
      const double lam2(sqr(s2) + sqr(mm1*mm1) + sqr(mm2*mm2)
                        - 2.*s2*mm1*mm1 - 2.*s2*mm2*mm2 - 2.*mm1*mm1*mm2*mm2);
      if (s2 <= 0. || lam2 <= 0.) { ++m_cxppfail; continue; }
      const double lcm(0.5*sqrt(lam2/s2));
      const double sgn(m_bornmomenta[0][3] < 0 ? -1. : 1.);
      Vec4D B1(lcm*sqrt(1.+mm1*mm1/sqr(lcm)), 0., 0.,  sgn*lcm);
      Vec4D B2(lcm*sqrt(1.+mm2*mm2/sqr(lcm)), 0., 0., -sgn*lcm);
      { Poincare bq2(Q2); bq2.BoostBack(B1); bq2.BoostBack(B2); }
      Vec4D_Vector b2{B1, B2, bp[2], bp[3]};
      m_sp = s2;
      MakeProp();
      MakePropT(b2);
      okb = ComixBornAmplitude(b2, Bc);
      BornAmplitude(b2, Bh, -1., -1., -1);
      okh = HandOnePhotonAmplitude(bp, kmap, hel, H);
      m_sp = sp_save;
      MakeProp();
      MakePropT(m_pceex);
    }
    if (!okb || !okh) { ++m_cxppfail; continue; }

    const int fmaskx(Amplitude::NHel() - 1);
    /*
      beta_1 = M_1 - S B, with BOTH terms Comix's own: M_1 from the 2 -> 3
      evaluation and B from the 2 -> 2 one, so no hand-coded object is ever
      added to a Comix one. No second real evaluation and no scaled photon -
      the Born is the one already computed, which is what keeps it fixed while
      the eikonal grows.
    */
    /*
      beta_1 = M_1(k) - s(k) M_0, with the eikonal term taken from Comix's OWN
      amplitude rather than written out.

      Under k -> lambda k the eikonal scales EXACTLY: the polarisation vector
      depends only on the photon's direction, which scaling does not change,
      and every p.k in the denominator picks up one factor of lambda. So
      s(lambda k) = s(k)/lambda, and

          s(k) M_0  =  lambda M_1(lambda k)   as lambda -> 0.

      That removes the gauge reference from the problem entirely. Writing the
      eikonal out instead means choosing polarisation vectors, and CEEX's
      reference (m_zeta, m_eta, m_b) is not Comix's (the fixed gauge vector in
      Amplitude::SetGauge) - the two differ by a little-group phase, which is
      exactly the per-photon phase measured earlier (|r| = 1/e^2,
      r_- = -conj(r_+)). Both terms here are Comix's own amplitude at CEEX's
      own arguments, so no such phase exists to compensate.

      The fermions are held fixed between the two evaluations and the
      propagator pole is pinned to the same scale, so the Born content is
      common and cancels rather than merely nearly cancelling.
    */
    for (int f = 0; f < nh; ++f)
      Csub.m_A[f] = C.m_A[f] - lam*Cs.m_A[f];
    m_cxcsub[j] = Csub;
    double nrmH(0.), num(0.), den(0.);
    for (int f = 0; f < nh; ++f) nrmH += std::norm(H.m_A[f]);
    nrmH = sqrt(nrmH);
    int nal(0);
    for (int f = 0; f < nh; ++f) {
      const Complex Bcf(Bc.m_A[f ^ (m_comixflip & fmaskx)]);
      const Complex Bhf(m_e*m_e*Bh.m_A[f]);
      if (std::abs(Bcf) == 0. || std::abs(H.m_A[f]) <= 1e-8*nrmH) {
        m_cxratio[j].m_A[f] = Complex(1., 0.);
        continue;
      }
      const Complex a(Bhf/Bcf);                 // Comix -> CEEX, per helicity
      const Complex cs(a*Csub.m_A[f]);
      m_cxratio[j].m_A[f] = cs/H.m_A[f];
      num += std::norm(cs - H.m_A[f]);
      den += std::norm(cs) + std::norm(H.m_A[f]);
      ++nal;
    }
    /*
      Does the Comix remainder behave like the hand-coded one? Both are the
      IR-SUBTRACTED beta_1: they must vanish together as the photon softens,
      not diverge like the full amplitude. Printed against the photon energy,
      which is what separates "wrong normalisation" from "wrong object".
    */
    { static const bool d1(ATOOLS::Settings::GetMainSettings()["CEEX"]
                           ["BETA1_CHECK"].Get<int>() != 0);
      static int nd(0);
      if (d1 && nd < 14) { ++nd;
        double nc(0.), nh2(0.), ncs(0.);
        for (int f = 0; f < nh; ++f) {
          nc  += std::norm(C.m_A[f]);
          ncs += std::norm(Csub.m_A[f]);
          nh2 += std::norm(H.m_A[f]);
        }
        /*
          The YFS soft theorem at THIS photon, squared so no convention phase
          enters: |M_1|^2 -> S~ |B|^2 as the photon softens. Both sides are
          Comix's own - C from the 2 -> 3 evaluation, Bc from the 2 -> 2 one -
          and the eikonal is CEEX's, whose polarisation phase drops out of the
          modulus. If this does not go to 1, the eikonal reference is wrong and
          no subtraction built from it can work.
        */
        double snum(0.), sden(0.);
        for (int f = 0; f < nh; ++f) {
          snum += std::norm(C.m_A[f]);
          sden += std::norm(Sk*Bc.m_A[f ^ (m_comixflip & fmaskx)]);
        }
        /*
          The same test NLO_Base::CheckRealCollinearSub uses for the squared
          reals: divide by the CLOSED-FORM massive eikonal of the two charged
          initial legs, which is pure kinematics - no spinor convention, no
          CEEX soft factor, nothing that could be wrong in the same way as the
          thing being tested. The eikonal theorem then forces the ratio to be
          FLAT, and flat is the verdict; it does not have to tend to zero or
          to one. A departure from flat is the generator, not the physics.
        */
        const Vec4D q1(m_pceex[0]), q2(m_pceex[1]), kk(m_allphotons[j]);
        const double mel(m_flavs[0].Mass());
        const double pk1(q1*kk), pk2(q2*kk), p12(q1*q2);
        const double Scl(pk1 != 0. && pk2 != 0.
                         ? 2.*p12/(pk1*pk2) - mel*mel/(pk1*pk1)
                                            - mel*mel/(pk2*pk2) : 0.);
        double nb(0.);
        for (int f = 0; f < nh; ++f)
          nb += std::norm(Bc.m_A[f ^ (m_comixflip & fmaskx)]);
        std::cerr<<"@@@ SOFTRAT Ek="<<m_allphotons[j][0]
                 <<" x="<<(2.*m_allphotons[j][0]/sqrt(m_s))
                 <<" |M1|^2/(S~|B|^2)="<<(sden>0.? snum/sden : -1.)
                 <<" eik_flat="<<((Scl != 0. && nb > 0.)
                                  ? snum/(Scl*nb) : -1.)<<std::endl;
        /*
          The same subtraction without a second amplitude evaluation: the Born
          is the one already in hand, Comix's own at the 2 -> 2 point, and the
          eikonal is CEEX's soft factor. |C| = |Sk||Bc| in the soft limit
          (measured: the ratio is 1/e^2 with soft carrying an explicit e^2,
          and m_e IS e here, so the factors cancel). Printed as moduli AND as
          the subtracted modulus, because the two differ by CEEX's photon
          polarisation phase, which is a convention and cannot cancel here.
        */
        double nsb(0.), ndir(0.);
        for (int f = 0; f < nh; ++f) {
          // Bc is RAW Comix; C already carries m_comixnorm. Without it here
          // the subtraction term is half the size and |M1|/|S B| sits at 0.5.
          const Complex SB(RealNorm()*Sk*Bc.m_A[f ^ (m_comixflip & fmaskx)]);
          nsb  += std::norm(SB);
          ndir += std::norm(C.m_A[f] - SB);
        }
        std::cerr<<"@@@ BETA1DIR Ek="<<m_allphotons[j][0]
                 <<" |M1|="<<sqrt(nc)<<" |S*B|="<<sqrt(nsb)
                 <<" |M1|/|S*B|="<<(nsb>0.? sqrt(nc/nsb) : -1.)
                 <<" |M1-S*B|/|H|="<<(nh2>0.? sqrt(ndir/nh2) : -1.)<<std::endl;
        std::cerr<<"@@@ BETA1 Ek="<<m_allphotons[j][0]
                 <<" |M1_comix|="<<sqrt(nc)
                 <<" |Csub|="<<sqrt(ncs)
                 <<" |H_hand|="<<sqrt(nh2)
                 <<" Csub/H="<<(nh2>0.? sqrt(ncs/nh2) : -1.)<<std::endl; } }
    if (nal == 0) { ++m_cxppfail; continue; }
    m_cxratiook[j] = 1;
    if (den > 0.) m_cxppdev += sqrt(num/den);
    ++m_cxppn;
  }
  m_sp = sp_save;
  MakeProp();
  MakePropT(m_pceex);
}


/*!
  The real-subtraction validation figure, for CEEX's beta_1.

  Same construction as NLO_Base::CheckRealSub, which produced the equivalent
  plots for the squared reals: scan the photon energy downward at fixed
  direction and write |beta_1|/beta_0 against E_gamma. The physics content is
  that the subtracted real must VANISH as the photon softens - the curve going
  to zero is the statement that the infrared subtraction is right, and a curve
  that flattens instead is a subtraction that has missed a term.

  Two curves, on one normalisation so they overlay: CEEX's hand-coded beta_1
  and the one built from Comix's amplitude as
  beta_1 = M_1(k) - lambda M_1(lambda k).
*/
void Ceex_Base::Beta1Scan()
{
  static bool done(false);
  if (done || m_allphotons.empty() || m_pceex.size() < 4) return;
  if (!m_comixcalibrated) DeriveComixMap();
  if (m_comixflip < 0) return;
  YFS::Real_Correction *prov(RealProvider(1));
  if (prov == NULL || p_bornproc == NULL) return;
  const Vec4D k0(m_allphotons[0]);
  if (k0[0] <= 0.) return;
  done = true;

  std::string tag;
  for (size_t i(0); i < m_flavs.size(); ++i) { tag += "_";
                                               tag += m_flavs[i].IDName(); }
  std::ofstream fc(("CeexBeta1_comix" + tag + ".txt").c_str());
  std::ofstream fh(("CeexBeta1_hand"  + tag + ".txt").c_str());
  fc << std::setprecision(10);
  fh << std::setprecision(10);

  Vec4D_Vector bp(m_pceex.begin(), m_pceex.begin()+4);
  const int nh(Amplitude::NHel());
  const int fmaskx(Amplitude::NHel() - 1);
  const int hel(m_PhoHel.empty() ? 1 : m_PhoHel[0]);
  const double sp_save(m_sp);

  // the Born both curves are normalised to, once
  m_sp = m_svarQ; MakeProp(); MakePropT(bp);
  Amplitude Bc, Bh;
  const bool okb(ComixBornAmplitude(bp, Bc, NULL, m_svarQ));
  BornAmplitude(bp, Bh, -1., -1., -1);
  m_sp = sp_save; MakeProp(); MakePropT(m_pceex);
  if (!okb) return;
  double nbc(0.), nbh(0.);
  for (int f = 0; f < nh; ++f) {
    nbc += std::norm(Bc.m_A[f]);
    nbh += std::norm(m_e*m_e*Bh.m_A[f]);
  }
  if (!(nbc > 0.) || !(nbh > 0.)) return;
  /*
    What the two Borns differ by, before anything is divided by them. The
    derived map says |B_comix| = 2 e^2 |B_hand| exactly (the 2 being
    1/sqrt(initial spin states)), so this should print 2 - and if it does,
    normalising the two beta_1 curves to their OWN Born puts a spurious factor
    of 2 between them, which is a property of the plot and not of the
    amplitudes.
  */
  msg_Info()<<"CEEX: beta_1 scan Born check: |B_comix|/(2 e^2 |B_hand|) = "
            <<(sqrt(nbc)/(2.*sqrt(nbh)))<<"  (e^2 |B_hand| = "<<sqrt(nbh)
            <<", |B_comix| = "<<sqrt(nbc)<<")\n";
  /*
    Both curves are normalised to the SAME Born from here, Comix's, so their
    ratio is |beta_1_comix|/|beta_1_hand| and nothing else. Dividing each by
    its own Born is what put an unexplained offset between them.
  */
  const double nbref(nbc);

  /*
    A clean logarithmic scan rather than CheckRealSub's k = k/i, which
    compounds: after 25 steps it is at 1e-6 GeV and after 190 at 1e-176, so
    the range cannot be chosen and the tail is numerical noise. Here the
    number of decades is a setting, which is what makes "how far down does it
    stay linear" an answerable question.
  */
  static const int ndec(ATOOLS::Settings::GetMainSettings()["CEEX"]
                        ["BETA1_SCAN_DECADES"].SetDefault(8).Get<int>());
  static const int nppd(ATOOLS::Settings::GetMainSettings()["CEEX"]
                        ["BETA1_SCAN_PPD"].SetDefault(4).Get<int>());
  for (int it(0); it <= ndec*nppd; ++it) {
    const Vec4D k(k0*pow(10., -double(it)/double(nppd)));
    const double xg(2.*k[0]/sqrt(m_s));
    // CheckRealSub stops at the ISR cutoff; without a floor k = k/i compounds
    // and the scan runs to 1e-176 GeV and NaN.
    if (xg <= 0.) break;

    /*
      beta_1 alone tends to a CONSTANT as the photon softens - the
      k-substituted Born carries one power of E against the 1/(2 p.k), so the
      two cancel. That is not what YFS requires to vanish. What enters the
      master formula is beta_1(k)/s(k), the code's nrm = sProd/Sactu, and with
      s ~ 1/E that ratio vanishes linearly. So the soft factor is divided out
      here, which is also what makes this curve the same object as the one
      NLO_Base::CheckRealSub plots for the squared reals.

      Moduli only, so CEEX's polarisation convention in s(k) does not enter.
    */
    const Complex qrat(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
    const Complex Sk(Sfactor(bp[0], bp[1], k, hel)
                     + qrat*Sfactor(bp[2], bp[3], k, hel));
    const double aS(std::abs(Sk));
    if (!(aS > 0.)) continue;

    // --- Comix: beta_1 = M_1(k) - lambda M_1(lambda k), same fermions, same
    //     propagator scale, so the Born content is common and cancels
    const double lam(Min(0.1, 1e-4/Max(xg, 1e-30)));
    Amplitude C, Cs;
    Vec4D_Vector p1v{bp[0], bp[1], bp[2], bp[3], k};
    Vec4D_Vector p2v{bp[0], bp[1], bp[2], bp[3], lam*k};
    if (ComixRealAt(p1v, hel, C, m_svarQ) &&
        ComixRealAt(p2v, hel, Cs, m_svarQ)) {
      double nc(0.);
      for (int f = 0; f < nh; ++f)
        nc += std::norm(C.m_A[f] - lam*Cs.m_A[f]);
      fc << k[0] << "," << sqrt(nc/nbref)/aS << std::endl;
    }

    // --- the hand-coded beta_1 at the same point
    m_sp = m_svarQ; MakeProp(); MakePropT(bp);
    Amplitude H;
    const bool okh(HandOnePhotonAmplitude(bp, k, hel, H));
    m_sp = sp_save; MakeProp(); MakePropT(m_pceex);
    if (okh) {
      double nhh(0.);
      for (int f = 0; f < nh; ++f) nhh += std::norm(H.m_A[f]);
      // e^2 restores the coupling the hand-coded Born omits, so that both
      // numerators are in the same normalisation before the common Born
      fh << k[0] << "," << sqrt(nhh/nbref)/aS << std::endl;
    }
  }
  fc.close(); fh.close();
  msg_Info()<<"CEEX: beta_1 scan written to CeexBeta1_{comix,hand}"<<tag
            <<".txt\n";
}
