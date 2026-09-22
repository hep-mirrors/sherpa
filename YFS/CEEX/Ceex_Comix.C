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
  if (p_realprov == NULL || p_realprov->p_proc == NULL) return false;
  const size_t nl(m_flavs.size());
  if (m_allphotons.size() != 1 || m_pceex.size() < nl || m_PhoHel.empty())
    return false;

  // The real process must be this process plus one photon, in that order -
  // YFS_Process builds it by pushing a photon onto the final state, so it is,
  // but a mapped or reordered process would silently misindex every helicity.
  const Flavour_Vector &pf(p_realprov->p_proc->Flavours());
  if (pf.size() != nl + 1) return false;
  for (size_t i(0); i < nl; ++i) if (pf[i] != m_flavs[i]) return false;
  if (pf[nl].Kfcode() != kf_photon) return false;

  Vec4D_Vector p(nl + 1);
  for (size_t i(0); i < nl; ++i) p[i] = m_pceex[i];
  p[nl] = m_allphotons[0];

  // m_pceex uses the LAB pair, which is the only one that balances against
  // the photons (see the comment on m_pceex). Check it rather than trust it:
  // Comix at a non-conserving point returns a number, not an error.
  Vec4D bal(p[0] + p[1]);
  for (size_t i(2); i < p.size(); ++i) bal -= p[i];
  const double scale(sqrt(dabs(m_svarQ)) + 1.);
  double worst(0.);
  for (int c(0); c < 4; ++c) worst = Max(worst, dabs(bal[c]));
  if (worst > 1e-6*scale) { ++m_cxrfail; return false; }

  const std::vector<METOOLS::Spin_Amplitudes> *amps
    (p_realprov->ComixAmplitudes(p));
  if (amps == NULL || amps->empty()) { ++m_cxrfail; return false; }
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
  const size_t nhel(((size_t)1) << (nl + 1));
  if (sa.size() != nhel) {
    static bool warned(false);
    if (!warned) {
      warned = true;
      msg_Error()<<METHOD<<"(): expected "<<nhel<<" helicity entries for "
                 <<(nl+1)<<" legs, got "<<sa.size()
                 <<". CEEX: COMIX_REAL disabled for this run."<<std::endl;
    }
    m_comixreal = 0;
    return false;
  }

  // CEEX draws ONE photon helicity per event (MakePhotonHel) and both the
  // numerator and the crude denominator use it, so the average is right;
  // Comix has computed both, and we pick the drawn one. +1 -> eps+ -> 0.
  const int hg(m_PhoHel[0] > 0 ? 0 : 1);
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
  const int fl(m_comixflip);
  // NOTE the mask is per-process: the photon bit sits above the fermion legs,
  // so a value fitted at 4 legs does NOT carry over to 6.
  const size_t fmask((((size_t)1) << nl) - 1);
  const int    hgf(hg ^ ((fl >> nl) & 1));
  const size_t nf(((size_t)1) << nl);
  for (size_t f(0); f < nf; ++f)
    m_comixM1.m_A[f] = m_comixnorm
      * sa[FlatIndex(f ^ (fl & fmask), hgf, (int)nl)];
  // The raw table, kept so the flip scan can look at maps other than the one
  // in force without a second Comix evaluation.
  for (int h = 0; h <= 1; ++h)
    for (size_t f(0); f < nf; ++f)
      m_comixraw[h].m_A[f] = m_comixnorm * sa[FlatIndex(f, h, (int)nl)];
  m_comixdrawnhel = hg;

  m_havecomixreal = true;
  return true;
}


void Ceex_Base::ApplyComixReal()
{
  m_havecomixreal = false;
  if (!m_comixreal) return;
  // Only the O(alpha) real, i.e. exactly one photon - see the file header.
  if (m_allphotons.size() != 1) return;
  if (!m_comixcalibrated) DeriveComixMap();
  if (!m_comixreal) return;   // calibration may have switched it off
  if (!FetchComixReal()) return;

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
    if (sh > 0.) { m_cxrnormsum += sc/sh; ++m_cxrn; }
    if (mden > 0.) m_cxrmetsum += sqrt(mnum/mden);
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
    const int nlg((int)m_flavs.size());
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


void Ceex_Base::ReportComixReal() const
{
  if (!m_comixreal) return;
  msg_Info()<<"CEEX: COMIX_REAL supplied the O(alpha) real on "<<m_cxrn
            <<" one-photon events";
  if (m_cxrn > 0)
    msg_Info()<<", mean |M1_comix|^2/|A1_hand|^2 = "<<(m_cxrnormsum/m_cxrn)
              <<", mean per-helicity magnitude mismatch "
              <<(m_cxrmetsum/m_cxrn);
  if (m_cxrfail)
    msg_Info()<<"; "<<m_cxrfail<<" events fell back to the hand-coded real";
  msg_Info()<<".\n";
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
  if (ndone >= 10) return;
  if (!m_havecomixreal || m_allphotons.size() != 1) return;
  if (p_realprov == NULL) return;
  const Vec4D p1(m_pceex[0]), p2(m_pceex[1]), k0(m_allphotons[0]);
  const Vec4D Q0(m_pceex[2] + m_pceex[3]);
  const double m3(m_flavs[2].Mass()), m4(m_flavs[3].Mass());
  Vec4D d3(m_pceex[2]);
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
    const int nlg((int)m_flavs.size());
    const int fmaskx((1 << nlg) - 1);
    const int nh(Amplitude::NHel());
    const int hgi(((hg > 0 ? 0 : 1) ^ ((fl >> nlg) & 1)));
    double sc(0.), shd(0.), scall(0.);
    for (int f = 0; f < nh; ++f) {
      const Complex C(sa[FlatIndex((size_t)(f ^ (fl & fmaskx)), hgi, nlg)]);
      sc  += std::norm(m_comixnorm*C);
      shd += std::norm(soft*B.m_A[f]);
    }
    for (size_t i(0); i < sa.size(); ++i) scall += std::norm(sa[i]);
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
        const double C(std::abs(m_comixnorm*
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
  if (needflip) m_comixflip = c.mask | (m_comixphoflip ? (1 << (int)nl) : 0);

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
