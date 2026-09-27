/*!
  \file Ceex_Comix.C

  The O(alpha) REAL correction taken from Comix instead of from hand-coded
  spinor products.

*/

#include "YFS/CEEX/Ceex_Base.H"
#include "YFS/NLO/NLO_Base.H"   // MapMomenta, for the reduced beta_1 kinematics
#include "YFS/NLO/Real_Correction.H"

#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"
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
           msg_Debugging()<<"CEEX: COMIX_REAL declines "<<ng<<" photons: "          \
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
  /*
    YFS: ME_PROBE - the momenta CEEX hands Comix for its one-photon M_1: the
    photon-lepton angle in units of that lepton's m/E, the lepton energy,
    the balance, and sum |A|^2 over both photon planes (map-free).
  */
  { static const bool mp(ATOOLS::Settings::GetMainSettings()["YFS"]
                         ["ME_PROBE"].SetDefault(0).Get<int>() != 0);
    if (mp && ng == 1 && amps && !amps->empty()) {
      const Vec4D &g(p.back());
      double best(1e99), el(0.);
      for (size_t i(2); i < nl && i < m_flavs.size(); ++i) {
        if (!m_flavs[i].IsChargedLepton()) continue;
        const double ct(Vec3D(p[i])*Vec3D(g)/(Vec3D(p[i]).Abs()*Vec3D(g).Abs()));
        const double th(acos(Max(-1., Min(1., ct)))/(m_flavs[i].Mass()/p[i][0]));
        if (th < best) { best = th; el = p[i][0]; }
      }
      double all(0.);
      for (size_t j(0); j < (*amps)[0].size(); ++j) all += std::norm((*amps)[0][j]);
      std::ostringstream o;
      o<<std::setprecision(6)<<"@@@ CEEXPP x="<<2.*g[0]/sqrt(m_s)
       <<std::setprecision(10)<<" th_pp="<<best<<" El_pp="<<el
       <<" bal="<<worst<<" sumA2="<<all<<" g="<<g<<"\n";
      std::cerr<<o.str();
    } }
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

    Hence 26 = 11010b. Change it against the CEEXSOFTFLIP and CEEXFLIP
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
    m_comixM1.m_A[f] = RealNorm() * ComixPhotonCoupling(ng) * sa[(f ^ ffl) | (hgf << np)];
  /*
    The raw table, kept so the flip scan can look at maps other than the one
    in force without a second Comix evaluation. It has one plane per photon
    helicity, so it only means anything for a single photon; above that the
    scan is not available and the diagnostics that read it say so.
  */
  if (ng == 1)
    for (int h = 0; h <= 1; ++h)
      for (size_t f(0); f < nf; ++f)
        m_comixraw[h].m_A[f] = RealNorm() * ComixPhotonCoupling(ng) * sa[f | ((size_t)h << np)];
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
    // ApplyComixReal no longer substitutes anything (the one-photon
    // substitution was removed with the flux fix); it is a CXREAL diagnostic
    // at n = 1 only, so above one photon there is simply nothing to do.
    return;
  }
  if (!m_comixcalibrated) DeriveComixMap();
  if (!m_comixreal) return;   // calibration may have switched it off
  const bool got(FetchComixReal());
  CountComixReal(m_allphotons.size(), got);
  if (!got) return;

  static const bool cxchk(ATOOLS::Settings::GetMainSettings()["CEEX"]["COMIX_CHECK"].Get<int>() != 0);
  const Amplitude hand1(m_AmpExpo1);   // the partition loop's A1
  if (m_b1trace && m_allphotons.size() == 1) {
    // the two routes to the same one-photon amplitude, side by side
    Amplitude Mex;
    Vec4D_Vector pp(m_pceex); pp.push_back(m_allphotons[0]);
    const double rn(RealNorm());
    double nf(0.), na(0.), nx(0.), nd(0.);
    const bool okx(rn > 0. && ComixRealAt(pp, m_PhoHel[0], Mex, -1.));
    for (int f(0); f < Amplitude::NHel(); ++f) {
      nf += std::norm(m_comixM1.m_A[f]); na += std::norm(m_AmpExpo1.m_A[f]);
      if (okx) { const Complex ex(m_cxbalign.m_A[f]*Mex.m_A[f]/rn);
                 nx += std::norm(ex); nd += std::norm(ex - m_comixM1.m_A[f]); }
    }
    msg_Error()<<std::setprecision(8)<<"CXREAL |M1_fetch|="<<sqrt(nf)
             <<" |A1_loop|="<<sqrt(na)<<" |align*Mex/rn|="<<sqrt(nx)
             <<" |diff|/|fetch|="<<(nf>0.? sqrt(nd/nf) : -1.)
             <<" rn="<<rn<<" m_comixnorm="<<m_comixnorm<<" m_cxrnorm="<<m_cxrnorm
             <<" normexact="<<m_normexact<<" flip="<<m_comixflip
             <<" hel="<<m_PhoHel[0]<<" drawnhel="<<m_comixdrawnhel<<std::endl;
  }
  /*
    This used to REPLACE m_AmpExpo1 by V*M1, V = AmpBornVirt/AmpExpo0, and
    ran right before MakeRho. That was the one-photon matching before
    ComixInfraredSubtracted_1_0 existed; now the partition loop's beta_1 IS
    Comix's one-photon amplitude and A1 = M1 holds at one photon
    algebraically, so the replacement had nothing to do - until
    NO_PSEUDOFLUX: 2 put the flux into rho_0 but not rho_1, when V became
    1/flux and the substitution silently divided every one-photon A1 by it
    (seed-3 point: real/Born -0.586 -> -0.918). What remains is the
    comparison below, now a closure check: the mean |M1_comix|^2/|A1|^2 it
    reports must be 1.
  */
  const int nh(Amplitude::NHel());

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
    msg_Error()<<"CEEXFLIP inuse="<<m_comixflip<<" metric="<<m0metric
             <<" best="<<bestmask<<" bestmetric="<<bestmetric<<std::endl;
    const double ecm((m_momenta[0]+m_momenta[1]).Mass());
    const double x(ecm > 0. ? 2.*m_allphotons[0][0]/ecm : 0.);
    msg_Error()<<"CEEXCMP nphot=1 x="<<x
             <<" sum_hand="<<sh<<" sum_comix="<<sc
             <<" ratio="<<(sc != 0. ? sh/sc : 0.)
             <<" worst_elem_reldiff="<<worst<<std::endl;
    msg_Error()<<"CEEXSOFT x="<<x
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
              msg_Error()<<"CEEXHEL "<<a<<b<<c<<d
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
  msg_Debugging()<<"CEEX: COMIX_REAL supplied the real amplitude on "<<m_cxrn
            <<" events";
  if (m_cxrn > 0)
    msg_Debugging()<<", mean |M1_comix|^2/|A1_hand|^2 = "<<(m_cxrnormsum/m_cxrn)
              <<", mean per-helicity magnitude mismatch "
              <<(m_cxrmetsum/m_cxrn);
  if (m_cxrfail)
    msg_Debugging()<<"; "<<m_cxrfail<<" events fell back to the hand-coded real";
  msg_Debugging()<<".\n";
  if (m_perphoton) {
    msg_Debugging()<<"CEEX: COMIX_REAL_PER_PHOTON replaced the one-photon"
              <<" amplitude on "<<m_cxppn<<" photons";
    if (m_cxppn > 0)
      msg_Debugging()<<", mean |C-H|/|C,H| = "<<(m_cxppdev/m_cxppn);
    if (m_cxppfail)
      msg_Debugging()<<"; "<<m_cxppfail<<" photons kept the hand-coded amplitude";
    msg_Debugging()<<".\n";
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
    msg_Debugging()<<"  "<<n<<" photon"<<(n>1?"s":"")<<": "<<t<<" taken, "
              <<f<<" refused";
    if (t > 0 && n < m_cxrnormbyn.size())
      msg_Debugging()<<", <|M_comix|^2/|A_hand|^2> = "<<(m_cxrnormbyn[n]/t)
                <<", <per-helicity mismatch> = "<<(m_cxrmetbyn[n]/t);
    msg_Debugging()<<"\n";
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
      sc  += std::norm(RealNorm()*ComixPhotonCoupling()*C);
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
        Poincare pRot(m_bornmomenta[0], Vec4D(0., 0., 0., 1.));
        const double lcm(0.5*sqrt(lam2/s2));
        const double sgn(m_bornmomenta[0][3] < 0 ? -1. : 1.);
        Vec4D B1(lcm*sqrt(1.+mm1*mm1/sqr(lcm)), 0., 0.,  sgn*lcm);
        Vec4D B2(lcm*sqrt(1.+mm2*mm2/sqr(lcm)), 0., 0., -sgn*lcm);
        Poincare bq2(Q2); 
        pRot.RotateBack(B1); pRot.RotateBack(B2); 
        bq2.BoostBack(B1); bq2.BoostBack(B2);
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
          msg_Error()<<"CEEXSOFTPHASE lam="<<lam<<" nlive="<<nl2
                   <<" |r|=["<<amin<<","<<amax<<"]"
                   <<" arg(r)=["<<phmin<<","<<phmax<<"]"
                   <<" hel="<<hg<<" cos_kf="<<ct<<" cos_kbeam="<<cb
                   <<" phi_k="<<atan2(k[2], k[1])
                   <<" sQ="<<s2<<std::endl;
        } else {
          msg_Error()<<"CEEXSOFTPHASE lam="<<lam
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
        const double C(std::abs(RealNorm()*ComixPhotonCoupling()*
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
    msg_Error()<<"CEEXSOFTFLIP lam="<<lam<<" inuse="<<m_comixflip
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
      msg_Error()<<"CEEXEIKFLAT lam="<<lam
               <<" x="<<(2.*k[0]/sqrt(m_s))
               <<" ratio="<<((Scl != 0. && nb > 0.) ? sc/(Scl*nb) : -1.)
               <<std::endl; }
    msg_Error()<<"CEEXSOFTSCAN lam="<<lam
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
  /*
    The t-channel propagators too, at the SAME point. Only MakeProp() was
    reset here, so for Bhabha the hand-coded Born carried whatever t-channel
    exchange was last set (or none, on the first call) while Comix's Born was
    evaluated at bp with its full t-channel: the calibration then measured
    N = 4.3 instead of 2, fell back to that number, and chose the helicity
    flip mask by matching two different amplitudes. MakePropT returns early
    for anything but Bhabha.
  */
  MakePropT(bp);
  BornAmplitude(bp, hand);
  double cxme2(0.);
  const bool ok(ComixBornAmplitude(bp, cx, &cxme2));
  m_sp = sp_save;
  MakeProp();
  if (m_pceex.size() >= 4) MakePropT(m_pceex);
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
  /*
    With YFS: USE_MODEL_ALPHA 0 the hand-coded Born carries alpha(0) and
    Comix's the model's, so N is sqrt(nspin) * alpha_model/alpha(0): 2.0774
    in the G_mu scheme. m_rescale_alpha = alpha(0)/alpha_model (1 otherwise).
  */
  const double expected(nspin > 0. ? m_rescale_alpha/sqrt(nspin) : 0.);
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

  msg_Debugging()<<METHOD<<"(): Comix -> CEEX map derived from the Born:\n"
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
    M1.m_A[f] = RealNorm() * ComixPhotonCoupling()
      * sa[FlatIndex((size_t)(f ^ (fl & fmaskx)), hgi, nlg)];
  return true;
}

bool Ceex_Base::ComixRealShifted(const Vec4D_Vector &pp, int hel,
                                 Amplitude &M1, const PropShifts &shifts)
{
  YFS::Real_Correction *prov(RealProvider(1));
  if (prov == NULL) return false;
  const int nlg(Amplitude::s_nlegs);
  if (pp.size() != m_flavs.size() + 1) return false;
  const std::vector<METOOLS::Spin_Amplitudes> *amps
    (prov->ComixAmplitudesShifts(pp, shifts));
  if (amps == NULL || amps->empty()) return false;
  const METOOLS::Spin_Amplitudes &sa((*amps)[0]);
  if ((int)sa.size() != (1 << (nlg + 1))) return false;
  const int fl(m_comixflip), fmaskx((1 << nlg) - 1), nh(Amplitude::NHel());
  const int hgi((hel > 0 ? 0 : 1) ^ ((fl >> nlg) & 1));
  for (int f = 0; f < nh; ++f)
    M1.m_A[f] = RealNorm() * ComixPhotonCoupling()
      * sa[FlatIndex((size_t)(f ^ (fl & fmaskx)), hgi, nlg)];
  return true;
}

Ceex_Base::PropShifts Ceex_Base::StageShifts(int iphot) const
{
  PropShifts sh;
  if (m_stage.size() != m_allphotons.size()) return sh;
  // one entry per stage: (the stage's Comix mask, its photons of this
  // partition). A W decay stage's mask is the W's daughters, so its photons
  // reach that W's line and nothing else; the production stage keeps the
  // initial-leg mask (StageShiftMask).
  for (size_t g(0); g < m_stagelegs.size(); ++g) {
    const size_t mask(StageShiftMask((int)g));
    if (mask == 0) continue;
    sh.push_back(std::make_pair(mask, StagePhotonSum((int)g, iphot)));
  }
  return sh;
}

/*
  The reduced legs of this partition, for the space-like exchange lines.

  Lines with one initial leg and part of the final state - Bhabha's
  t-channel boson, the t/u-channel electrons of e+e- -> gamma gamma, a
  t-channel neutrino - cannot be placed by a stage, and the stage rule left
  them at the unreduced invariant in every partition while M_1 shifted them
  (COMIX::Amplitude::SetPropShifts, class (b), has the measurements). They
  are put at the partition's REDUCED invariant: the beams rebuilt back to
  back at X_wp along the event's axis in X's frame, the pair at Y_wp along
  its own direction (LegsAt) - the construction the crude Born momenta come
  from, so that the crude and beta_0 carry the same t for the partition the
  generator actually labelled, and exact in both collinear limits. The
  spinors stay physical; only the pole moves. This passes, per Born leg,
  theta_i (p~_i - p_i) with theta = -1 incoming, +1 outgoing (Comix's
  all-outgoing current momenta); SetPropShifts detects the exchange lines
  from their leg content and sums the entries over their legs, so no
  process-specific mask is named here. iphot < 0: the Born
  (BornLegsAt(m_PXvec)); iphot >= 0: M_1 of that photon (PartitionLegs). At
  one photon PartitionLegs returns the physical legs and every entry is
  zero, so the n = 1 amplitude is untouched. CEEX: TCHANNEL_SHIFT: 0
  passes nothing.
*/
/*
  CEEX: TCHANNEL_SHIFT. The reduced legs come from LegsAt, a 2 -> 2
  construction: beams at X, "the radiating pair" m_if1/m_if2 at Y. Beyond
  2 -> 2 that pair is just the first two final fermions - the two neutrinos
  of e+e- -> nu_mu nubar_e mu+ e- - and the shift tilts the reduced beams and
  drags the single-W t-channel gamma* towards t = 0: YFS.CEEX 0.443 pb +-75%
  with it, 0.0404 +-1.9% without, YFS.NLO 0.0435 +-4.7% (161 GeV,
  2026-09-27, NOTES-w-stages-2026-09-27.md 1.3). So the default (-1) keeps it
  for 2 -> 2, where it is validated (Bhabha, gamma gamma), and drops it above.
*/
bool Ceex_Base::ExchangeLineShiftsOn() const
{
  static const int mode(ATOOLS::Settings::GetMainSettings()["CEEX"]
                        ["TCHANNEL_SHIFT"].Get<int>());
  if (mode >= 0) return mode != 0;
  return m_flavs.size() == 4;
}

void Ceex_Base::AddExchangeLineShifts(int iphot, PropShifts &sh) const
{
  if (!ExchangeLineShiftsOn() || m_pceex.size() < 4
      || m_flavs.size() != m_pceex.size()) return;
  Vec4D_Vector pb;
  const bool ok(iphot < 0 ? BornLegsAt(m_PXvec, pb) : PartitionLegs(iphot, pb));
  if (!ok || pb.size() != m_pceex.size()) return;
  for (size_t i(0); i < m_pceex.size(); ++i) {
    const Vec4D d((i < 2 ? -1. : 1.)*(pb[i] - m_pceex[i]));
    sh.push_back(std::make_pair((((size_t)1) << i)
                                | PHASIC::Process_Base::s_propshiftleg, d));
  }
}

bool Ceex_Base::BornHasExchangeLine()
{
  std::map<const PHASIC::Process_Base*, int>::const_iterator
    it(m_exchline.find(p_bornproc));
  if (it != m_exchline.end()) return it->second != 0;
  if (p_bornproc == NULL || m_pceex.size() < 4
      || m_flavs.size() != m_pceex.size()) return false;
  /*
    Leg entries only reach currents with one initial leg and a proper subset
    of the final legs; without such a current SetPropShifts sets every shift
    to zero and the amplitude is bit for bit the unshifted one. The probe
    shift is O(1e-3) of the beam energy, far above rounding where it acts.
  */
  Amplitude A0, A1;
  PropShifts none, leg;
  const double e(1e-3*m_pceex[0][0]);
  for (size_t i(0); i < m_pceex.size(); ++i)
    leg.push_back(std::make_pair((((size_t)1) << i)
                                 | PHASIC::Process_Base::s_propshiftleg,
                                 Vec4D(e, 0.3*e*(i+1), -0.2*e, 0.5*e)));
  if (!ComixBornShifted(m_pceex, A0, none) || !ComixBornShifted(m_pceex, A1, leg))
    return false;                          // not cached: try again next event
  int has(0);
  for (int f(0); f < Amplitude::NHel(); ++f)
    if (A0.m_A[f] != A1.m_A[f]) { has = 1; break; }
  m_exchline[p_bornproc] = has;
  msg_Info()<<"CEEX: Born "<<(p_bornproc ? p_bornproc->Name() : std::string("?"))
            <<(has ? " has" : " has no")<<" space-like exchange line"
            <<(has ? "s" : "")<<"."<<std::endl;
  return has != 0;
}

Vec4D Ceex_Base::PartitionShift(int iphot, bool reducing) const
{
  Vec4D d;
  if (m_stage.size() != m_allphotons.size()) return d;
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if ((int)i != iphot
        && (m_stagereduces[m_stage[i]] != 0) == reducing) d += m_allphotons[i];
  return d;
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
    M1.m_A[f] = RealNorm() * ComixPhotonCoupling()
      * sa[FlatIndex((size_t)(f ^ (fl & fmaskx)), hgi, nlg)];
  kmap = kk;
  return true;
}



/*
  Closure test for the lambda subtraction, with NO hand-coded reference.
  M1 = s beta_0 + beta_1 defines s, so solving
        s_implied(hel) = (M1(hel) - beta_1(hel)) / beta_0(hel)
  gives the soft factor the Comix amplitudes actually obey. If the lambda
  subtraction is right, s_implied is the SAME for every live helicity (the
  eikonal is a scalar, it cannot depend on hel) and its ratio to CEEX's own
  Sfactor is the convention constant.
*/
void Ceex_Base::Beta1Closure()
{
  if (m_allphotons.empty()) return;
  /*
    The SOFTEST photon of the event, not "events with one photon". One-photon
    events are inherently HARD - a single emission that survives took most of
    the energy, measured x_gamma 0.76-0.94 - and the identity being tested
    needs beta_0 at the reduced kinematics, which only coincides with the
    unreduced one when the photon is soft. ComixRealAt evaluates the
    one-photon process for whichever photon it is handed, so any photon of any
    event is a valid configuration for this test.
  */
  size_t js(0);
  for (size_t i(1); i < m_allphotons.size(); ++i)
    if (m_allphotons[i][0] < m_allphotons[js][0]) js = i;
  const int nh(Amplitude::NHel());
  const int msk(m_comixflip & (Amplitude::NHel()-1));
  const double rn(RealNorm());
  if (!(rn > 0.)) return;
  Amplitude B0, M1, B1;
  if (!ComixBornAmplitude(m_pceex, B0, NULL, m_svarQ)) return;
  Vec4D_Vector pp(m_pceex); pp.push_back(m_allphotons[js]);
  if (!ComixRealAt(pp, m_PhoHel[js], M1, m_svarQ)) return;
  if (!ComixBeta1At(m_allphotons[js], m_PhoHel[js], m_svarQ, B1)) return;
  const Complex qratio(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
  const Complex sceex(Sfactor(m_pceex[0], m_pceex[1],
                              m_allphotons[js], m_PhoHel[js])
                      + qratio*Sfactor(m_pceex[m_if1], m_pceex[m_if2],
                                       m_allphotons[js], m_PhoHel[js]));
  double nb(0.);
  for (int f(0); f < nh; ++f) nb += std::norm(B0.m_A[f]);
  nb = sqrt(nb);
  if (!(nb > 0.) || std::abs(sceex) == 0.) return;
  int nlive(0); Complex first(0.,0.); double worst(0.);
  for (int f(0); f < nh; ++f) {
    if (std::abs(B0.m_A[f]) < 1e-3*nb) continue;
    const Complex m1((M1.m_A[f ^ msk])/rn), b1((B1.m_A[f ^ msk])/rn);
    const Complex si((m1 - b1)/B0.m_A[f]);
    if (nlive++ == 0) first = si;
    else if (std::abs(first) > 0.)
      worst = Max(worst, std::abs(si/first - 1.));
  }
  if (nlive < 2 || std::abs(first) == 0.) return;
  const double xg(m_s > 0. ? 2.*m_allphotons[js][0]/sqrt(m_s) : -1.);
  msg_Error()<<"B1CLOSE xg="<<xg<<" nlive="<<nlive
           <<" hel_spread="<<worst
           <<" s_implied_over_Sceex="<<std::abs(first/sceex)
           <<" arg="<<std::arg(first/sceex)<<std::endl;
}


/*
  The YFS theorem fixes the amplitude-level soft factor with no reference to
  anyone's conventions: summing s(k,hel) s*(k,hel) over the two photon
  helicities must give the ordinary (squared-level) eikonal

     Stilde(k) = - sum_ij Q_i Q_j th_i th_j (p_i.p_j)/((p_i.k)(p_j.k))
                 + sum_i  Q_i^2 m_i^2 / (p_i.k)^2

  up to the overall e^2. So the ratio of the two IS the normalisation of
  Sfactor, measured rather than assumed - and once s is pinned,
  beta_1 = M1 - s M0 is fully determined.
*/
void Ceex_Base::SoftNormCheck(const Vec4D &k)
{
  if (m_pceex.size() < 4) return;
  // the charged legs with their YFS signs: incoming +, outgoing -
  std::vector<Vec4D>  pl;
  std::vector<double> Q, th;
  pl.push_back(m_pceex[0]);     Q.push_back(m_qe); th.push_back(+1.);
  pl.push_back(m_pceex[1]);     Q.push_back(-m_qe); th.push_back(+1.);
  pl.push_back(m_pceex[m_if1]); Q.push_back(m_qf); th.push_back(-1.);
  pl.push_back(m_pceex[m_if2]); Q.push_back(-m_qf); th.push_back(-1.);
  double st(0.);
  for (size_t i(0); i < pl.size(); ++i) {
    const double pik(pl[i]*k);
    if (pik == 0.) return;
    for (size_t j(0); j < pl.size(); ++j) {
      const double pjk(pl[j]*k);
      if (pjk == 0.) return;
      st -= Q[i]*Q[j]*th[i]*th[j]*(pl[i]*pl[j])/(pik*pjk);
    }
    st += Q[i]*Q[i]*pl[i].Abs2()/(pik*pik);
  }
  // sum over the two photon helicities of |s|^2, with CEEX's own Sfactor
  const Complex qratio(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
  double ss(0.);
  for (int h(-1); h <= 1; h += 2) {
    const Complex sa(Sfactor(m_pceex[0], m_pceex[1], k, h)
                     + qratio*Sfactor(m_pceex[m_if1], m_pceex[m_if2], k, h));
    ss += std::norm(sa);
  }
  if (st == 0.) return;
  msg_Error()<<"SOFTNORM xg="<<(m_s>0.?2.*k[0]/sqrt(m_s):-1.)
           <<" sum_hel|s|^2="<<ss<<" Stilde="<<st
           <<" ratio="<<(ss/st)<<std::endl;
}

/*
  The paper's own soft-limit test (hep-ph/0006359, the EEX discussion around
  eq.(single-initial)): hold the spectators fixed, scale k_j -> 0, and

      beta_1(k_j; X_wp) / ( s(k_j) beta_0(X_wp) )  ->  0

  beta_1 is formed DIRECTLY from its definition, M_1 - s beta_0, not by the
  lambda subtraction - the point is to test whether our M_1, s and beta_0 are
  mutually consistent, i.e. whether they share one prescription. A ratio that
  tends to a nonzero constant IS the mismatch, and its value measures it.

  X_wp is recomputed at every lambda: for an ISR photon X = P - lambda k, which
  is the "same extrapolation for beta_0 and beta_1" the paper requires. Holding
  X fixed while scaling k would itself manufacture a mismatch.

  beta_0 carries the pseudo-flux X^2/(p_c+p_d)^2, as eq.(305) defines it.
*/
void Ceex_Base::SoftLimitTest()
{
  /*
    The soft theorem on the BEAM-REDUCED configuration - the one beta_1 is
    actually evaluated on for a multi-photon event.

    An earlier version of this test held the spectators fixed and boosted the
    pair, and the theorem held: |M_1 - s B|/|s B| fell linearly in lambda. But
    that is not the configuration ComixInfraredSubtracted_1_0 uses. There the
    beams are reduced to S = (final legs) + k, so scaling k moves the beams
    too, and the measured subtraction floors at 10-30% instead of vanishing.

    So scan lambda here on exactly that construction, and report alongside it
    the angle between the reduced beam and the physical one. The beams are
    rebuilt back-to-back in S's rest frame about the axis found by boosting
    p_a into that frame; once the spectator photons carry transverse momentum
    that axis is rotated, and a rotation - unlike the rescaling that motivated
    this construction - does NOT leave the eikonal invariant. If r tracks the
    angle, that is the fault.
  */
  if (m_allphotons.empty() || !m_cxbalignok) return;
  const int nh(Amplitude::NHel()), flip(m_comixflip & (nh-1));
  const double rn(RealNorm());
  if (!(rn > 0.)) return;
  const Vec4D k0(m_allphotons[0]);
  const int hel(m_PhoHel[0]);
  std::string out;
  for (int e(0); e <= 6; ++e) {
    const double lam(pow(10., -(double)e));
    const Vec4D k(lam*k0);
    Vec4D_Vector pb;
    if (!ReducedBeams(k, pb)) continue;
    const Complex stot(TotalEikonal(pb, k, hel));
    if (std::abs(stot) == 0.) continue;
    Vec4D_Vector pp(pb);
    pp.push_back(k);
    Amplitude C, M1;
    if (!ComixBornAmplitude(pb, C, NULL, -1., -1.)) continue;
    if (!ComixRealAt(pp, hel, M1, -1.)) continue;
    double nb(0.), ne(0.);
    for (int f(0); f < nh; ++f) {
      const Complex sb(stot*C.m_A[f ^ flip]);
      nb += std::norm(M1.m_A[f]/rn - sb); ne += std::norm(sb);
    }
    if (!(ne > 0.)) continue;
    const Vec3D a(pb[0]), b(m_pceex[0]);
    const double ab(a.Abs()*b.Abs());
    const double ang(ab > 0. ? acos(Min(1., Max(-1., (a*b)/ab))) : -1.);
    out += " lam=" + ToString(lam) + " r=" + ToString(sqrt(nb/ne))
         + " ang=" + ToString(ang);
  }
  if (!out.empty())
    msg_Error()<<"SOFTLIM nphot="<<m_allphotons.size()<<out<<std::endl;
}

/*
  A real phase-space point whose invariant is s' = X^2.

  CEEX's beta_0(X_wp) has been the Born on the PHYSICAL spinors with only the
  propagator pole moved to X - a configuration, as the comment in
  InfraredSubtractedME_0_0 says, that no momenta realise. That is workable for
  hand-written spinor algebra. It is not workable next to a Comix amplitude,
  which necessarily lives at a real point: at one photon the two coincide, and
  beyond one photon no placement of the beta_1 subtraction can match both. The
  measurement either way was exact closure at n=1 with the multi-photon
  cancellation broken, or the reverse.

  So build the point instead. Beams back-to-back on shell with total X, the
  radiating pair back-to-back with total X minus the spectator final legs, and
  the spectators left alone. Directions are taken from the event, in X's frame,
  so nothing is invented; only the scale moves. The propagator then follows
  from the momenta rather than being pinned, which is the point.
*/
bool Ceex_Base::LegsAt(const Vec4D &X, const Vec4D &Y, Vec4D_Vector &pb) const
{
  pb = m_pceex;
  const double X2(X.Abs2());
  const double ma(m_flavs[0].Mass()), mb(m_flavs[1].Mass());
  if (X2 <= sqr(ma + mb)) return false;
  Poincare cmsX(X);
  /*
    The beam axis in X's frame. The generator builds its s' point with the
    reduced beams back to back along the LAB z axis as the pure boost
    Poincare(X) carries it (measured: the Comix Born at the generator's own
    point, m_plabmom, is m_born to six digits at every x, and its transverse
    momentum is parallel to X's). The direction of the boosted e- differs
    from that by aberration once the photons carry transverse momentum - for
    a 44 GeV photon at 147 degrees the two axes gave Borns 1400x apart on
    e+e- -> gamma gamma, whose Born is 1/(1 - cos^2). An s-channel Born
    hardly sees the axis; a space-like exchange line does, and the crude
    the CEEX weight divides by is the generator's, so its axis is the one
    to use. CEEX: REDUCED_AXIS: 0 keeps the boosted-beam axis.
  */
  static const int axis(ATOOLS::Settings::GetMainSettings()["CEEX"]
                        ["REDUCED_AXIS"].SetDefault(1).Get<int>());
  Vec3D da(0., 0., 1.);
  if (axis == 0) {
    Vec4D ra(m_pceex[0]);
    cmsX.Boost(ra);
    const double rap(Vec3D(ra).Abs());
    if (!(rap > 0.)) return false;
    da = Vec3D(ra)/rap;
  }
  const double kala(sqr(X2 - ma*ma - mb*mb) - 4.*ma*ma*mb*mb);
  if (kala <= 0.) return false;
  const double pa(sqrt(kala)/(2.*sqrt(X2)));
  Vec4D b1(sqrt(pa*pa + ma*ma), pa*da), b2(sqrt(pa*pa + mb*mb), -pa*da);
  cmsX.BoostBack(b1); cmsX.BoostBack(b2);
  pb[0] = b1; pb[1] = b2;
  // the radiating pair balances Y, along the direction it has in Y's frame
  const double Y2(Y.Abs2());
  const double m3(m_flavs[m_if1].Mass()), m4(m_flavs[m_if2].Mass());
  if (Y2 <= sqr(m3 + m4)) return false;
  Poincare cmsY(Y);
  Vec4D r3(m_pceex[m_if1]);
  cmsY.Boost(r3);
  const double r3p(Vec3D(r3).Abs());
  if (!(r3p > 0.)) return false;
  const Vec3D d3(Vec3D(r3)/r3p);
  const double kalf(sqr(Y2 - m3*m3 - m4*m4) - 4.*m3*m3*m4*m4);
  if (kalf <= 0.) return false;
  const double pf(sqrt(kalf)/(2.*sqrt(Y2)));
  Vec4D q3(sqrt(pf*pf + m3*m3), pf*d3), q4(sqrt(pf*pf + m4*m4), -pf*d3);
  cmsY.BoostBack(q3); cmsY.BoostBack(q4);
  pb[m_if1] = q3; pb[m_if2] = q4;
  return true;
}

bool Ceex_Base::GeneratorBornAt(const Vec4D &R, Vec4D_Vector &pb) const
{
  const size_t nl(m_flavs.size());
  if (nl < 4 || m_prefsr.size() != nl) return false;
  Vec4D Q;
  for (size_t i(2); i < nl; ++i) Q += m_prefsr[i];
  const double s2(R.Abs2());
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (!(s2 > sqr(m1 + m2)) || !(Q.Abs2() > 0.) || !(Q[0] > 0.)) return false;
  Poincare toQ(Q);
  std::vector<Vec3D> q; std::vector<double> mm2; double msum(0.);
  for (size_t i(2); i < nl; ++i) {
    Vec4D qi(m_prefsr[i]); toQ.Boost(qi); q.push_back(Vec3D(qi));
    mm2.push_back(sqr(m_flavs[i].Mass())); msum += m_flavs[i].Mass();
  }
  const double M(sqrt(s2));
  if (msum >= M) return false;
  auto etot = [&](double xi) { double e(0.);
    for (size_t j(0); j < q.size(); ++j) e += sqrt(mm2[j] + xi*xi*q[j].Sqr());
    return e; };
  double lo(0.), hi(1.);
  while (etot(hi) < M && hi < 1e6) hi *= 2.;
  for (int it(0); it < 200; ++it) { const double mid(0.5*(lo+hi)); (etot(mid) < M ? lo : hi) = mid; }
  const double xi(0.5*(lo+hi));
  pb.assign(nl, Vec4D());
  const double sgn(m_bornmomenta.size() > 0 && m_bornmomenta[0][3] < 0. ? -1. : 1.);
  const double lam(0.5*sqrt(Max(0., sqr(s2 - m1*m1 - m2*m2) - 4.*m1*m1*m2*m2)/s2));
  pb[0] = Vec4D(sqrt(lam*lam + m1*m1), 0., 0.,  sgn*lam);
  pb[1] = Vec4D(sqrt(lam*lam + m2*m2), 0., 0., -sgn*lam);
  for (size_t i(2); i < nl; ++i) {
    const Vec3D v(xi*q[i-2]);
    pb[i] = Vec4D(sqrt(mm2[i-2] + v.Sqr()), v);
  }
  return true;
}

bool Ceex_Base::BornLegsAt(const Vec4D &X, Vec4D_Vector &pb) const
{
  Vec4D Y(X);
  for (size_t i(2); i < m_pceex.size(); ++i)
    if (i != m_if1 && i != m_if2) Y -= m_pceex[i];
  return LegsAt(X, Y, pb);
}

/*
  The legs to hand Comix for photon iphot in the CURRENT partition.

  Where the other photons' momentum goes is fixed by the partition, not by a
  kinematic preference: a photon the partition counts as initial-state came
  out of the beams, one it counts as final-state went into the radiating
  pair. Placing them so, the beams carry X_wp(others) and the final state
  carries X_wp(others) - k, and the one-photon amplitude at that point has
  its ISR-attachment propagator at (P' - k)^2 = X_wp with this photon counted
  initial and its FSR-attachment propagator at P'^2 = X_wp with it counted
  final - which is what the CEEX expansion asks of M_1^I and M_1^F in this
  partition. Comix cannot separate the two attachments and does not need to.

  At one photon there are no others and this is the physical point up to
  round-off. The spectators (a Higgs, say) are never touched.
*/
bool Ceex_Base::PartitionLegs(int iphot, Vec4D_Vector &pb) const
{
  if (iphot < 0 || iphot >= (int)m_allphotons.size()
      || m_stage.size() != m_allphotons.size()) return false;
  Vec4D X(m_pceex[0] + m_pceex[1]);
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if ((int)i != iphot && m_stagereduces[m_stage[i]]) X -= m_allphotons[i];
  Vec4D Y(X - m_allphotons[iphot]);
  for (size_t i(2); i < m_pceex.size(); ++i)
    if (i != m_if1 && i != m_if2) Y -= m_pceex[i];
  return LegsAt(X, Y, pb);
}

/*
  The legs to hand Comix for photon k, when the event has other photons too.

  Comix is a Berends-Giele recursion: it contracts currents on the momenta it
  is given and never checks that they balance. Hand it the event's legs plus
  ONE of n photons and it returns a number for a configuration missing the
  other n-1 photons' momenta - not gauge invariant, not an amplitude. At one
  photon the set happens to conserve, which is why a one-photon test looks
  healthy while the same code is nonsense at twelve.

  Where the missing momentum belongs is not a free choice. Taking it out of the
  final pair (boost the pair to balance P-k) moves the pair's invariant up to
  near sqrt(s), and at the Z pole that changes the Born by the width - the
  measured result was 5e4 pb against 888. It belongs in the BEAMS, because
  that is where it physically went: the spectator ISR photons reduced the
  energy entering the hard process, which is what X_wp means.

  So the final legs and the photon stay exactly as they are, and the beams are
  rebuilt back-to-back on shell with invariant (sum of final legs + k)^2. Two
  things then hold that nothing else gives:

    - the pair is untouched, so beta_0 here is the same Born, at the same
      spinors, that InfraredSubtractedME_0_0 adds to the partition sum;
    - the eikonal is almost unchanged, because the leg current p/(p.k) is
      invariant under p -> x p, and reducing a beam is a rescaling plus the
      small rotation that the spectators' transverse momentum forces.

  The earlier warning against "rebuilding the beams" was about rebuilding them
  for the BORN, where an O(E_gamma) shift is the size of beta_1 itself. Here
  the Born is left alone and only the legs radiating into it are reduced.

  At one photon the final legs plus k already are P, so this returns the
  physical beams and the n=1 stream is unchanged.
*/
bool Ceex_Base::ReducedBeams(const Vec4D &k, Vec4D_Vector &pb) const
{
  pb = m_pceex;
  Vec4D S(k);
  for (size_t i(2); i < m_pceex.size(); ++i) S += m_pceex[i];
  const double S2(S.Abs2());
  const double ma(m_flavs[0].Mass()), mb(m_flavs[1].Mass());
  if (S2 <= sqr(ma + mb)) return false;
  Poincare cms(S);
  Vec4D pa(m_pceex[0]);
  cms.Boost(pa);
  const double pap(Vec3D(pa).Abs());
  if (!(pap > 0.)) return false;
  const Vec3D dir(Vec3D(pa)/pap);
  const double kal(sqr(S2 - ma*ma - mb*mb) - 4.*ma*ma*mb*mb);
  if (kal <= 0.) return false;
  const double rs(sqrt(S2)), ps(sqrt(kal)/(2.*rs));
  Vec4D pa2(sqrt(ps*ps + ma*ma), ps*dir), pb2(sqrt(ps*ps + mb*mb), -ps*dir);
  cms.BoostBack(pa2); cms.BoostBack(pb2);
  pb[0] = pa2; pb[1] = pb2;
  return true;
}

/*
  The total eikonal of a given leg configuration: sum over every external leg
  of w * SfactorLeg, w = Q * theta. This is the same sum CalculateSfactors
  checks its stage currents against (eq. 6.4), so on the physical legs it
  returns sum_g m_Sfac[g][j] identically. It is recomputed on the REDUCED legs
  here so that the subtraction below cancels against the amplitude Comix
  actually evaluated, rather than against the one the partition loop tabulated.
*/
Complex Ceex_Base::TotalEikonal(const Vec4D_Vector &p, const Vec4D &k,
                                int hel)
{
  Complex tot(0., 0.);
  for (size_t l(0); l < m_flavs.size() && l < p.size(); ++l) {
    const double w(m_flavs[l].Charge() * (l < 2 ? -1. : +1.));
    if (w != 0.) tot += w * SfactorLeg(p[l], k, hel);
  }
  return tot;
}

Complex Ceex_Base::StageEikonal(int stage, const Vec4D_Vector &p,
                                const Vec4D &k, int hel, int iphot)
{
  // every leg of the stage, a reconstructed resonance included: its
  // momentum is built from its daughters in p (StageLegMomentum)
  return StageCurrent(stage, iphot, k, hel, p);
}

/*
  beta_1 for one photon, from Comix, with the eikonal subtracted - and with no
  partition decomposition of its own.

  The YFS theorem gives the subtraction directly:

      beta_1(k) = M_1(k) - s_tot(k) B,

  where s_tot is the eikonal of EVERY leg and B the Born on the same legs. That
  object vanishes as k -> 0 identically, needs no stage assignment, no
  per-stage propagator scale, and no pseudo-flux.

  The partition sum is not what beta_1 needs; it is what beta_0 needs. Summing
  over wp resums the propagator shift - which X the resonance sees once some
  photons are counted as initial-state - and that is a beta_0 statement. KKMC
  splits beta_1 across the stages as well, which forces a pseudo-flux on beta_0
  and a compensating (1-CKine) on the final-state beta_1 so the two cancel.
  Both are artefacts of that split. Taking beta_1 from the amplitude, neither
  is needed:

      M_n  ~  sum_wp prod s^wp beta_0(X_wp)  +  sum_j prod_{i!=j} s_tot(k_i)
                                                * beta_1(k_j).

  Earlier attempts to keep the stage split are what produced the scan over
  where to put the flux; every placement was wrong because the object being
  patched should not exist.

  beta_1 is therefore partition independent and is evaluated at Comix's own
  propagator scales, not at m_sp. Only its weight, the eikonal product of the
  OTHER photons, varies over partitions. Adding it on one designated stage of
  photon j visits each assignment of the others exactly once:
  sum_{wp: wp_j = 0} prod_{i!=j} s^{wp_i} = prod_{i!=j} s_tot(k_i).
*/
bool Ceex_Base::ComixInfraredSubtracted_1_0(const Vec4D &k, int hel,
                                            int iphot)
{
  if (!m_cxbalignok || iphot < 0) return false;
  if (iphot >= (int)m_stage.size() || (int)m_Sfac.size() != m_nstages
      || iphot >= (int)m_Sfac[m_stage[iphot]].size()) return false;
  /*
    Below the eikonal threshold beta_1 vanishes and what is left of the
    subtraction is the difference of two nearly equal large numbers. A
    twelve-photon event is mostly such photons. Zero is both the honest answer
    and the stable one.
  */
  /*
    Only ONE partition may carry M_1, so that summing over partitions visits
    each assignment of the other photons exactly once; the beta_0 subtraction,
    by contrast, has to happen in EVERY partition, because that is where _0_0
    put it. Hence a flag rather than an early return.

    Which stage carries M_1 is arbitrary - it cancels - but it must be one that
    exists for every photon. Index 0 is not that: the ordering comes from
    BuildStages and need not put the initial state first, and for a neutral
    final state (e+e- -> nu nu) index 0 is the non-radiating stage, so every
    photon returned with no contribution and beta_1 silently vanished. Name the
    stage by what it is.
  */
  /*
    A collapsed (fixed-stage, x < SOFT_PARTITION_CUT) photon is never
    enumerated, so it never visits the reducing stage: with the rule below
    alone it received the beta_0 subtraction in every partition and M_1 in
    none. That is infrared-unsafe whenever the photon is above BETA1_XCUT in
    y = 2k.S/S^2 while below the collapse cut in x = 2E/sqrt(s) - which is
    the rule after a hard ISR photon, S^2 << s. Measured at 250 GeV, mu mu:
    a 0.07 GeV photon (x = 5.6e-4, y = 1.6e-3) next to a 111 GeV
    radiative-return photon gave |beta_1| = 56 |A_0| in every helicity and a
    CEEX weight 3163 times the crude; eight of the nine events above 20 were
    of this kind and CEEX came out 3x Born+real with a 48% error.

    For a fixed photon every partition visits each assignment of the OTHERS
    exactly once, so M_1 is added in every partition, against the TOTAL
    eikonal times the partition Born - the very product m_sactu put into
    beta_0 for it. M_1's attachment lines sit at X_{j in I} and X_{j in F}
    by the shifts; the collapsed Born is at X_fixed, an O(k) difference the
    collapse already accepts.
  */
  const bool fixed(m_fixedstage.size() == m_stage.size()
                   && m_fixedstage[iphot] != 0);
  const bool addm1(fixed || m_stagereduces[m_stage[iphot]]);
  const double rn(RealNorm());
  if (!(rn > 0.)) return false;
  /*
    What M_1 is evaluated on. Writing out the O(alpha^1) CEEX expansion and
    summing the two partitions that differ only in where THIS photon sits
    gives, per assignment of the others,

      prod_{i!=j} s_i^{wp_i} [ M_1^I(k_j; X_{j in I}) + M_1^F(k_j; X_{j in F})
                              - s_I(k_j) beta_0(X_{j in I})
                              - s_F(k_j) beta_0(X_{j in F}) ],

    with M_1^{I,F} at the PHYSICAL spinors and the propagator at the
    partition's X (KKMC, hep-ph/0006359). Comix reproduces exactly that from
    the physical legs plus k_j when the initial-side lines carry the other
    initial-state photons and the final-side lines the other final-state
    ones (BETA1_LEGS: 0, ComixRealShifted): the line with the photon on the
    initial state then sits at X_{j in I}, the one with it on the final state
    at X_{j in F}, and the soft limit is s(physical legs) times the very
    Born the partition loop added.

    BETA1_LEGS: 1 keeps the earlier construction for comparison: a balanced
    real point with the other initial-state photons taken out of the beams
    and the other final-state ones absorbed by the pair (PartitionLegs). It
    has the right propagators but the wrong spinors: for a radiative-return
    photon at x = 0.85 the amplitude reduces to the eikonal times the Born
    with UNREDUCED beam spinors to 4%, and to 2.6x the reduced-beam Born.
    With that Born as beta_0 the nu nu cross section came out 7x the NLO.

    Before either, the beams were reduced by EVERY other photon in every
    partition - M_1 as if the others were all initial-state - while the
    subtraction ran over both assignments. Measured on an n = 2 point with
    two final-state photons at the Z pole: real/Born +2.26 against KKMC's
    -0.17.
  */
  static const int legsmode(ATOOLS::Settings::GetMainSettings()["CEEX"]
                            ["BETA1_LEGS"].Get<int>());
  const bool rebuild(legsmode != 0 || m_redborn);
  Vec4D_Vector pb;
  Vec4D dI;
  PropShifts shifts;
  if (rebuild) { if (!PartitionLegs(iphot, pb)) return false; }
  else {
    pb = m_pceex;
    shifts = StageShifts(iphot);
    AddExchangeLineShifts(iphot, shifts); // space-like exchange lines
    dI = PartitionShift(iphot, true);   // the initial stage's photons, for S
  }
  /*
    Softness is measured against the system this photon is radiated from -
    the beams less the other initial-state photons - not against the beams:
    2k.S/S^2 is its energy in S's rest frame over half the invariant mass.
    Below the cut beta_1 vanishes and the whole photon is skipped - no M_1,
    no subtraction - which is what beta_1 = 0 means. S does not depend on
    where THIS photon sits, so the two partitions that together make one
    beta_1 term are cut together.
  */
  Vec4D S(rebuild ? pb[0] + pb[1] : m_pceex[0] + m_pceex[1] - dI);
  /*
    A photon on a W DECAY stage is radiated from that W, not from the
    beams: its softness is 2k.P_W/P_W^2 with P_W the W as produced
    (daughters + the partition's other decay photons + this one). Only
    with W stages; the ISR/FSR case keeps S, where X = Q + K_F anyway.
  */
  if (WStagesActive() && !rebuild && m_stage[iphot] != m_initstage
      && !fixed) {
    const int g(m_stage[iphot]);
    for (size_t l(0); l < m_stagelegs[g].size(); ++l)
      if (IsResonanceLeg(m_stagelegs[g][l].leg))
        S = StageSystemMomentum(g) + StagePhotonSum(g, -1);
  }
  const double S2(S.Abs2());
  const double y(S2 > 0. ? 2.*(k*S)/S2 : 0.);
  static const double xcut(ATOOLS::Settings::GetMainSettings()["CEEX"]
                           ["BETA1_XCUT"].Get<double>());
  if (!(y > xcut)) return true;     // m_b2on[iphot] stays 0: beta_1 = 0 here
  /*
    M_1 is needed on one partition per assignment of the others - the one
    with this photon on the reducing stage - and is evaluated only there.

    The subtraction is this partition's Born times the SAME soft factors that
    M_1 reduces to: the other photons' from the partition, and THIS photon's
    stage eikonal on the legs Comix was handed (the physical ones in the
    default mode). The Born is the one InfraredSubtractedME_0_0 added, taken
    from where it left it: the (1 - n) cancellation between the n
    subtractions and the one addition is exact only if they are the same
    numbers.
  */
  Amplitude M1;
  if (addm1) {
    Vec4D_Vector pp(pb);
    pp.push_back(k);
    const bool ok(rebuild ? ComixRealAt(pp, hel, M1, -1.)
                          : ComixRealShifted(pp, hel, M1, shifts));
    if (!ok) return false;
  }
  // the other photons' eikonal product, as a product - never as sProd/s_j,
  // which manufactures a divergence at the zeros of s_j
  Complex w(1., 0.);
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if ((int)i != iphot)
      w *= (m_sactu.size() == m_allphotons.size() ? m_sactu[i] : m_Sfac[m_stage[i]][i]);
  const Complex sj(fixed ? TotalEikonal(pb, k, hel)
                         : StageEikonal(m_stage[iphot], pb, k, hel, iphot));
  const int nh(Amplitude::NHel());
  /*
    CEEX: ORDER 2 - record exactly what this partition's beta_1 of this
    photon is, for beta_2's subtraction (Ceex_Beta2.C): the aligned M_1 if
    carried here, and the eikonal subtracted with. Stores only; the O(alpha)
    arithmetic below is untouched.
  */
  if (m_order == 2 && iphot < (int)m_b2on.size()) {
    m_b2on[iphot] = 1;
    m_b2sj[iphot] = sj;
    m_b2hasM1[iphot] = addm1 ? 1 : 0;
    if (addm1)
      for (int f(0); f < nh; ++f)
        m_b2M1[iphot].m_A[f] = m_cxbalign.m_A[f] * M1.m_A[f]/rn;
  }
  double nsub(0.), nm1(0.), nv(0.);
  for (int f(0); f < nh; ++f) {
    const Complex sub(w * sj * m_partborn0.m_A[f]);
    Complex v(-sub);
    if (addm1) {
      const Complex m1(w * m_cxbalign.m_A[f] * M1.m_A[f]/rn);
      v += m1;
      nm1 += std::norm(m1);
    }
    nsub += std::norm(sub); nv += std::norm(v);
    m_AmpExpo1.m_A[f]    += v;
    m_AmpBornReal.m_A[f] += v;
    if (iphot < (int)m_realphot.size()) m_realphot[iphot].m_A[f] += v;
    m_snapReal.m_A[f]    += v;
    m_beta10 += v;
  }
  if (m_b1trace) {
    std::string st;
    for (size_t i(0); i < m_stage.size(); ++i) st += ToString(m_stage[i]);
    const Complex sp(m_Sfac[m_stage[iphot]][iphot]);
    /*
      Does Comix's M_1 reduce to the total eikonal times Comix's own Born on
      the legs it was handed, with the same propagator treatment, and how does
      that Born compare with the partition Born being subtracted?
    */
    auto born_here = [&](const Vec4D_Vector &legs, Amplitude &C) {
      return rebuild ? ComixBornAmplitude(legs, C, NULL, -1., -1.)
                     : ComixBornShifted(legs, C, shifts); };
    auto real_here = [&](const Vec4D_Vector &legs, const Vec4D &kk,
                         Amplitude &M) {
      Vec4D_Vector pp(legs); pp.push_back(kk);
      return rebuild ? ComixRealAt(pp, hel, M, -1.)
                     : ComixRealShifted(pp, hel, M, shifts); };
    const int flip(m_comixflip & (nh-1));
    double rsoft(-1.), bratio(-1.);
    { Amplitude Cc;
      if (born_here(pb, Cc)) {
        const Complex stot(TotalEikonal(pb, k, hel));
        double nd(0.), ne(0.), nc(0.), nb(0.);
        for (int f(0); f < nh; ++f) {
          const Complex sb(stot*Cc.m_A[f ^ flip]);
          if (addm1) { nd += std::norm(M1.m_A[f]/rn - sb); ne += std::norm(sb); }
          nc += std::norm(m_cxbalign.m_A[f]*Cc.m_A[f ^ flip]);
          nb += std::norm(m_partborn0.m_A[f]);
        }
        if (addm1 && ne > 0.) rsoft = sqrt(nd/ne);
        if (nc > 0.) bratio = sqrt(nb/nc);
      } }
    /*
      The lambda scan only means something on legs that balance: scaling k
      alone leaves (1-lambda)k unaccounted for, which the shifts do not know
      about, so in the default mode only lambda = 1 is reported.
    */
    if (addm1) {
      std::string sc;
      for (int e(0); e <= (rebuild ? 3 : 0); ++e) {
        const double lam(pow(10., -(double)e));
        const Vec4D kl(lam*k);
        Vec4D_Vector pl(pb);
        if (rebuild) {
          Vec4D X(pb[0] + pb[1]), Y(X - kl);
          for (size_t i(2); i < m_pceex.size(); ++i)
            if (i != m_if1 && i != m_if2) Y -= m_pceex[i];
          if (!LegsAt(X, Y, pl)) continue;
        }
        Amplitude Cl, Ml;
        if (!born_here(pl, Cl) || !real_here(pl, kl, Ml)) continue;
        const Complex stl(TotalEikonal(pl, kl, hel));
        double nd(0.), ne(0.), nm(0.);
        Complex ov(0., 0.);
        for (int f(0); f < nh; ++f) {
          const Complex sb(stl*Cl.m_A[f ^ flip]);
          nd += std::norm(Ml.m_A[f]/rn - sb); ne += std::norm(sb);
          nm += std::norm(Ml.m_A[f]/rn);
          ov += (Ml.m_A[f]/rn) * std::conj(sb);
        }
        sc += " lam=" + ToString(lam) + " r=" + ToString(ne>0.? sqrt(nd/ne) : -1.)
            + " phi=" + ToString(std::arg(ov))
            + " |M1|=" + ToString(sqrt(nm)) + " |sC|=" + ToString(sqrt(ne));
      }
      std::cerr<<"B1LAM wp="<<st<<" j="<<iphot<<sc<<std::endl;
    }
    std::cerr<<std::setprecision(6)
             <<"B1TRACE wp="<<st<<" j="<<iphot<<" addm1="<<addm1
             <<" fixed="<<fixed<<" y="<<y<<" |w|="<<std::abs(w)
             <<" |s_cfg|="<<std::abs(sj)<<" |s_phys|="<<std::abs(sp)
             <<" s_cfg/s_phys="<<(std::abs(sp)>0.? sj/sp : Complex(0.,0.))
             <<" |sub|="<<sqrt(nsub)<<" |M1|="<<sqrt(nm1)
             <<" |beta1|="<<sqrt(nv)
             <<" |beta1|/|sub|="<<(nsub>0.? sqrt(nv/nsub) : -1.)
             <<" S2="<<S2<<" rsoft="<<rsoft<<" |B0wp|/|Ccfg|="<<bratio
             <<std::endl;
  }
  return true;
}

bool Ceex_Base::ComixBeta1At(const Vec4D &k, int hel,
                             const double propscale, Amplitude &B1)
{
  const double xg(m_s > 0. ? 2.*k[0]/sqrt(m_s) : 0.);
  static const double X_EIK(1e-4);
  // Already eikonal: beta_1 is below the precision of the difference, and it
  // vanishes there anyway. Zero is the honest answer, not a hand-coded one.
  if (!(xg > 10.*X_EIK)) {
    for (int f(0); f < Amplitude::NHel(); ++f) B1.m_A[f] = Complex(0., 0.);
    return true;
  }
  /*
    The definition, directly:  beta_1(k) = M_1(k) - s(k) beta_0.

    This used to be done as a lambda-subtraction, M_1(k) - lambda M_1(lambda k),
    which uses s(lambda k) = s(k)/lambda to reach the same place without ever
    writing s down. That detour existed because M_1 and s beta_0 did not agree
    in the soft limit, so the relative normalisation between them could not be
    trusted. The cause was a sign: SfactorLeg/Sfactor built the hel=-1 eikonal
    from Sminus rather than conj(Splus), the opposite convention for epsilon_-
    to the one Comix uses, and no |s|^2 or helicity-summed test can see it.

    With that fixed the soft limit holds - scaling k -> lambda k on conserving
    kinematics gives |M_1 - s beta_0|/|s beta_0| falling linearly in lambda
    (0.62, 0.064, 0.0065, 0.0013) until the 1/lambda cancellation runs out of
    double precision - so the subtraction can be done where it is defined.

    That also removes the lambda-subtraction's own defect: it evaluated
    M_1(lambda k) at m_pceex plus lambda k, a momentum set short by
    (1-lambda)k. Comix contracts currents on whatever it is handed, so that
    second term was neither gauge invariant nor an amplitude.

    Indexing: M_1 is Comix-ordered, beta_0 is CEEX-ordered, and B1 is returned
    Comix-ordered for the caller to map - hence the mask on beta_0 and the
    RealNorm on the eikonal term rather than on M_1.
  */
  Amplitude B0, M1;
  Vec4D_Vector pp(m_pceex);
  pp.push_back(k);
  /*
    The Born beta_1 subtracts must be the SAME object the partition sum adds
    with the soft factor, or rho_1 is not |M_1|^2 at one photon. With
    CEEX: BORN_AT_SPRIME: 1 the partition Born is the reduced-leg Born
    (BornLegsAt), so it is subtracted here too; the pinned physical-spinor
    Born was subtracted regardless, which broke the closure and gave
    e+e- -> gamma gamma a CEEX column of 7e7 pb with that switch on.
  */
  { static const bool sprime(ATOOLS::Settings::GetMainSettings()["CEEX"]
                             ["BORN_AT_SPRIME"].Get<int>() != 0);
    bool ok(false);
    if (sprime || m_redborn) {
      Vec4D_Vector pb;
      ok = BornLegsAt(m_PXvec, pb) && ComixBornAmplitude(pb, B0, NULL, -1., -1.);
    }
    if (!ok && !ComixBornAmplitude(m_pceex, B0, NULL, propscale)) return false; }
  if (!ComixRealAt(pp, hel, M1, propscale)) return false;
  const double rn(RealNorm());
  if (!(rn > 0.)) return false;
  const Complex qratio(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
  const Complex sf(Sfactor(m_pceex[0], m_pceex[1], k, hel)
                   + qratio*Sfactor(m_pceex[m_if1], m_pceex[m_if2], k, hel));
  // the paper's pseudo-flux X^2/(p_c+p_d)^2; exactly 1 for a single ISR photon
  const double pflux(propscale > 0. && m_svarQ > 0. ? propscale/m_svarQ : 1.);
  const int nh(Amplitude::NHel()), msk(m_comixflip & (nh-1));
  for (int f(0); f < nh; ++f)
    B1.m_A[f] = M1.m_A[f] - rn*sf*pflux*B0.m_A[f ^ msk];
  return true;
}

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
        msg_Error()<<"SOFTRAT Ek="<<m_allphotons[j][0]
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
        msg_Error()<<"BETA1DIR Ek="<<m_allphotons[j][0]
                 <<" |M1|="<<sqrt(nc)<<" |S*B|="<<sqrt(nsb)
                 <<" |M1|/|S*B|="<<(nsb>0.? sqrt(nc/nsb) : -1.)
                 <<" |M1-S*B|/|H|="<<(nh2>0.? sqrt(ndir/nh2) : -1.)<<std::endl;
        msg_Error()<<"BETA1 Ek="<<m_allphotons[j][0]
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

  Same construction as NLO_Base::CheckRealSub
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
  msg_Debugging()<<"CEEX: beta_1 scan Born check: |B_comix|/(2 e^2 |B_hand|) = "
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
  msg_Debugging()<<"CEEX: beta_1 scan written to CeexBeta1_{comix,hand}"<<tag
            <<".txt\n";
}
