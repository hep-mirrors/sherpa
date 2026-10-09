#include "YFS/ACRAIC/Acraic_Base.H"
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
    Flat Spin_Amplitudes index: bit i = leg i (0 -> +1, 1 -> -1), photons
    above the fermion legs. Same packing as Amplitude::Idx, so relabelling a
    set of legs is one FlipLegs.
  */
  inline size_t FlatIndex(size_t ferm, int hg, int nlegs)
  { return JoinLegs(ferm, (size_t)hg, nlegs); }

  //! The largest |component| of v: how far a momentum sum is from balance.
  double LargestComponent(const Vec4D &v)
  {
    double largest(0.);
    for (int c(0); c < 4; ++c) largest = Max(largest, dabs(v[c]));
    return largest;
  }

}


bool Acraic_Base::FetchComixReal()
{
  m_havecomixreal = false;
  const size_t nl(m_flavs.size());
  const size_t ng(m_allphotons.size());
  // Every refusal reports its reason, once per multiplicity.
  #define CXR_REFUSE(why)                                                     \
    do { if (m_cxrwhy.size() <= ng) m_cxrwhy.resize(ng+1);                    \
         if (m_cxrwhy[ng].empty()) {                                          \
           m_cxrwhy[ng] = (why);                                              \
           msg_Debugging()<<"ACRAIC: COMIX_REAL declines "<<ng<<" photons: "          \
                     <<(why)<<" (reported once)\n"; }                         \
         return false; } while (0)

  if (ng < 1 || m_pacraic.size() < nl || m_PhoHel.size() < ng)
    CXR_REFUSE("no photons, or momenta/helicities not set");
  // Amplitude holds only the 2^nl fermion helicities (photon helicities are
  // drawn); Comix's 2^(nl+ng) array is checked below, so only guard overflow.
  if (nl + ng >= 8*sizeof(size_t))
    CXR_REFUSE("helicity index would overflow");
  YFS::Real_Correction *prov(RealProvider(ng));
  if (prov == NULL)
    CXR_REFUSE("no real ME provider registered at this multiplicity");
  if (prov->p_proc == NULL)
    CXR_REFUSE("the provider has no process attached");

  // The real process must be this process plus ng photons, in that order;
  // a mapped or reordered process would misindex every helicity.
  const Flavour_Vector &pf(prov->p_proc->Flavours());
  if (pf.size() != nl + ng)
    CXR_REFUSE("the provider's process has the wrong number of legs");
  for (size_t i(0); i < nl; ++i) if (pf[i] != m_flavs[i])
    CXR_REFUSE("the provider's process has different fermions");
  for (size_t i(nl); i < nl + ng; ++i)
    if (pf[i].Kfcode() != kf_photon)
      CXR_REFUSE("the provider's extra legs are not photons");

  Vec4D_Vector p(nl + ng);
  for (size_t i(0); i < nl; ++i) p[i] = m_pacraic[i];
  for (size_t j(0); j < ng; ++j) p[nl + j] = m_allphotons[j];

  // Comix at a non-conserving point returns a number, not an error.
  Vec4D bal(p[0] + p[1]);
  for (size_t i(2); i < p.size(); ++i) bal -= p[i];
  const double scale(sqrt(dabs(m_svarQ)) + 1.);
  const double worst(LargestComponent(bal));
  if (worst > 1e-6*scale) { ++m_cxrfail; CXR_REFUSE("momenta do not balance"); }

  const std::vector<METOOLS::Spin_Amplitudes> *amps
    (prov->ComixAmplitudes(p));
  ProbeComixMomenta(p, amps, ng, nl, worst);
  if (amps == NULL || amps->empty())
    { ++m_cxrfail; CXR_REFUSE("Comix returned no amplitudes (KeepAmplitudes?)"); }
  if (amps->size() > 1) {
    static bool warned(false);
    if (!warned) {
      warned = true;
      msg_Error()<<METHOD<<"(): Comix returned "<<amps->size()
                 <<" colour configurations; ACRAIC uses the first. Reported "
                 <<"once."<<std::endl;
    }
  }
  const METOOLS::Spin_Amplitudes &sa((*amps)[0]);
  // Count legs with two helicity states: a scalar contributes a factor one
  // (e+e- -> H mu+ mu- gamma has 32 entries, not 64).
  const size_t np(Amplitude::s_nlegs);
  const size_t nhel(((size_t)1) << (np + ng));
  if (sa.size() != nhel) {
    static bool warned(false);
    if (!warned) {
      warned = true;
      msg_Error()<<METHOD<<"(): expected "<<nhel<<" helicity entries for "
                 <<(np+ng)<<" packed legs, got "<<sa.size()
                 <<". ACRAIC: COMIX_REAL disabled for this run."<<std::endl;
    }
    m_comixreal = false;
    return false;
  }

  // Pick the drawn photon helicities (MakePhotonHel) out of Comix's table;
  // photon j sits in bit nl+j, +1 -> 0.
  size_t hg(0);
  for (size_t j(0); j < ng; ++j)
    if (m_PhoHel[j] <= 0) hg = AddLeg(hg, j);
  /*
    The two codes label one leg of each fermion pair oppositely: a
    permutation of the array, not a phase. m_comixflip (bit i flips leg i)
    holds the fermion bits, derived at the Born (CalibrateComixMap); the
    photon bit is ACRAIC: COMIX_REAL_PHOTON_FLIP, one bit for all photons.
  */
  const size_t ffl(LowLegs((size_t)m_comixflip, np));
  const size_t gfl(m_comixphoflip ? AllLegs(ng) : 0);
  const size_t hgf(FlipLegs(hg, gfl));
  const size_t nf(NumLegSubsets(np));
  for (size_t f(0); f < nf; ++f)
    m_comixM1.m_A[f] = RealNorm() * ComixPhotonCoupling(ng) * sa[JoinLegs(FlipLegs(f, ffl), hgf, np)];
  // Raw table for the flip scan; one plane per photon helicity, so one photon only.
  if (ng == 1)
    for (int h = 0; h <= 1; ++h)
      for (size_t f(0); f < nf; ++f)
        m_comixraw[h].m_A[f] = RealNorm() * ComixPhotonCoupling(ng) * sa[JoinLegs(f, (size_t)h, np)];
  m_comixdrawnhel = (int)hg;

  m_havecomixreal = true;
  return true;
  #undef CXR_REFUSE
}


/*!
  ACRAIC: COMIX_REAL. Compares Comix's exact n-photon amplitude with the
  partition-loop A1 (CXREAL diagnostics); multiplicities without a real
  process fall back to the hand-coded path.
*/
void Acraic_Base::ApplyComixReal()
{
  m_havecomixreal = false;
  if (!m_comixreal) return;
  if (m_allphotons.empty()) return;
  // The per-photon route replaces this one: both write m_AmpExpo1.
  if (m_perphoton) return;

  static const bool multi(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["COMIX_REAL_MULTIPHOTON"].Get<bool>());
  if (m_allphotons.size() > 1 && !multi) {
    CountComixReal(m_allphotons.size(), false);
    return;
  }
  if (!m_comixcalibrated) DeriveComixMap();
  if (!m_comixreal) return;   // calibration may have switched it off
  const bool got(FetchComixReal());
  CountComixReal(m_allphotons.size(), got);
  if (!got) return;

  static const bool cxchk(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["COMIX_CHECK"].Get<int>() != 0);
  const Amplitude hand1(m_AmpExpo1);   // the partition loop's A1
  if (m_b1trace && m_allphotons.size() == 1) {
    // the two routes to the same one-photon amplitude, side by side
    Amplitude Mex;
    Vec4D_Vector pp(m_pacraic); pp.push_back(m_allphotons[0]);
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
  
  const int nh(Amplitude::NHel());

  // Norm-weighted, not worst-case: mass-suppressed entries make the worst case ~1.
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
  // Scan the 2^(nlegs+1) index maps on magnitudes only, so spinor phases
  // do not enter; norm-weighted metric.
  const int nlg(Amplitude::s_nlegs);   // helicity legs, not flavours
  int bestmask(0); double bestmetric(1e30); double m0metric(-1.);
  for (int mask(0); mask < NumLegSubsets<int>(nlg + 1); ++mask) {
    double num(0.), den(0.);
    const int hgm(FlipLegs(m_comixdrawnhel, LegBitOf(mask, nlg)));
    for (int f = 0; f < nh; ++f) {
      const double H(std::abs(hand1.m_A[f]));
      const double C(std::abs(m_comixraw[hgm].m_A[FlipLegs(f, LowLegs(mask, nlg))]));
      num += sqr(H - C); den += H*H + C*C;
    }
    const double met(den > 0. ? sqrt(num/den) : 0.);
    if (mask == m_comixflip) m0metric = met;
    if (met < bestmetric) { bestmetric = met; bestmask = mask; }
  }
  msg_Error()<<"ACRAICFLIP inuse="<<m_comixflip<<" metric="<<m0metric
           <<" best="<<bestmask<<" bestmetric="<<bestmetric<<std::endl;
  const double ecm((m_momenta[0]+m_momenta[1]).Mass());
  const double x(ecm > 0. ? 2.*m_allphotons[0][0]/ecm : 0.);
  msg_Error()<<"ACRAICCMP nphot=1 x="<<x
           <<" sum_hand="<<sh<<" sum_comix="<<sc
           <<" ratio="<<(sc != 0. ? sh/sc : 0.)
           <<" worst_elem_reldiff="<<worst<<std::endl;
  msg_Error()<<"ACRAICSOFT x="<<x
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
            msg_Error()<<"ACRAICHEL "<<a<<b<<c<<d
                     <<" absH="<<std::abs(H)<<" absC="<<std::abs(C)
                     <<" H=("<<H.real()<<","<<H.imag()<<")"
                     <<" C=("<<C.real()<<","<<C.imag()<<")"<<std::endl;
          }
}


void Acraic_Base::CountComixReal(size_t ng, bool ok)
{
  std::vector<long> &v(ok ? m_cxrnbyn : m_cxrfailbyn);
  if (v.size() <= ng) v.resize(ng+1, 0);
  ++v[ng];
}


void Acraic_Base::ReportComixReal() const
{
  if (!m_comixreal) return;
  msg_Debugging()<<"ACRAIC: COMIX_REAL supplied the real amplitude on "<<m_cxrn
            <<" events";
  if (m_cxrn > 0)
    msg_Debugging()<<", mean |M1_comix|^2/|A1_hand|^2 = "<<(m_cxrnormsum/m_cxrn)
              <<", mean per-helicity magnitude mismatch "
              <<(m_cxrmetsum/m_cxrn);
  if (m_cxrfail)
    msg_Debugging()<<"; "<<m_cxrfail<<" events fell back to the hand-coded real";
  msg_Debugging()<<".\n";
  if (m_perphoton) {
    msg_Debugging()<<"ACRAIC: COMIX_REAL_PER_PHOTON replaced the one-photon"
              <<" amplitude on "<<m_cxppn<<" photons";
    if (m_cxppn > 0)
      msg_Debugging()<<", mean |C-H|/|C,H| = "<<(m_cxppdev/m_cxppn);
    if (m_cxppfail)
      msg_Debugging()<<"; "<<m_cxppfail<<" photons kept the hand-coded amplitude";
    msg_Debugging()<<".\n";
    return;
  }
  // Per photon multiplicity: "refused" is where FetchComixReal declined and
  // the hand-coded path ran.
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
  Soft-limit check of the Comix real: M_1(k) -> e^2 [S_ini + S_fin] M_0, so

      R(lambda) = sum_hel |M1_comix(lambda k)|^2
                / sum_hel |e^2 (S_ini + S_fin) Born_hand|^2  ->  1 .

  Tests normalisation, helicity index map and photon helicity, not phases.
*/
void Acraic_Base::SoftProbe()
{
  static int ndone(0);
  static const int nmax(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                        ["SOFT_PROBE_EVENTS"].Get<int>());
  if (ndone >= nmax) return;
  if (!m_havecomixreal || m_allphotons.size() != 1) return;
  YFS::Real_Correction *p_realprov(RealProvider(1));
  if (p_realprov == NULL) return;
  const Vec4D p1(m_pacraic[0]), p2(m_pacraic[1]), k0(m_allphotons[0]);
  const Vec4D Q0(m_pacraic[m_if1] + m_pacraic[m_if2]);
  const double m3(m_flavs[m_if1].Mass()), m4(m_flavs[m_if2].Mass());
  Vec4D d3(m_pacraic[m_if1]);
  Poincare c0(Q0);
  c0.Boost(d3);
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

    // hand side: eikonal times the Born at this configuration
    Vec4D_Vector bp(4);
    bp[0]=p1; bp[1]=p2; bp[2]=q3; bp[3]=q4;
    m_sp = sQ;
    MakeProp(false);
    MakePropT(bp);
    Amplitude B;
    BornAmplitude(bp, B, -1., -1., -1);
    const Complex Si(Sfactor(p1, p2, k, hg));
    const Complex Sf(qratio*Sfactor(q3, q4, k, hg));
    const Complex soft(m_e*m_e*(Si + Sf));

    const int nlg(Amplitude::s_nlegs);   // helicity legs, not flavours
    const int nh(Amplitude::NHel());
    const int hgi(FlipLegs(HelicityBit<int>(hg), ComixPhotonFlip(nlg)));
    double sc(0.), shd(0.), scall(0.);
    for (int f = 0; f < nh; ++f) {
      const Complex C(sa[FlatIndex((size_t)FlipLegs(f, ComixLegFlips(nlg)), hgi, nlg)]);
      sc  += std::norm(RealNorm()*ComixPhotonCoupling()*C);
      shd += std::norm(soft*B.m_A[f]);
    }
    for (size_t i(0); i < sa.size(); ++i) scall += std::norm(sa[i]);

    // The same soft limit with phase, entirely in Comix's convention.
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
          nc2 += std::norm(sa[FlatIndex((size_t)FlipLegs(f, ComixLegFlips(nlg)),
                                        hgi, nlg)]);
        nc2 = sqrt(nc2);
        for (int f = 0; f < nh; ++f) {
          const Complex Cf(sa[FlatIndex((size_t)FlipLegs(f, ComixLegFlips(nlg)),
                                        hgi, nlg)]);
          const Complex Df(soft*Bc.m_A[FlipLegs(f, ComixLegFlips(nlg))]);
          if (std::abs(Cf) < 1e-4*nc2 || std::abs(Df) == 0.) continue;
          const Complex r(Cf/Df);
          const double a(std::abs(r)), ph(std::arg(r));
          if (nl2++ == 0) { amin = amax = a; phmin = phmax = ph; }
          else { amin = Min(amin,a); amax = Max(amax,a);
                 phmin = Min(phmin,ph); phmax = Max(phmax,ph); }
        }
        // Angles printed to correlate phases other than 0 or pi with kinematics.
        const double ct(q3.PSpat() > 0. && k.PSpat() > 0.
                        ? (Vec3D(q3)*Vec3D(k))/(q3.PSpat()*k.PSpat()) : 0.);
        const double cb(k.PSpat() > 0. ? k[3]/k.PSpat() : 0.);
        msg_Error()<<"ACRAICSOFTPHASE lam="<<lam<<" nlive="<<nl2
                 <<" |r|=["<<amin<<","<<amax<<"]"
                 <<" arg(r)=["<<phmin<<","<<phmax<<"]"
                 <<" hel="<<hg<<" cos_kf="<<ct<<" cos_kbeam="<<cb
                 <<" phi_k="<<atan2(k[2], k[1])
                 <<" sQ="<<s2<<std::endl;
      } else {
        msg_Error()<<"ACRAICSOFTPHASE lam="<<lam
                 <<" Comix Born unavailable"<<std::endl;
      }
    }
    // Index-map scan against the clean eikonal x Born reference; magnitudes only.
    int bm(0); double bmet(1e30), inusemet(-1.), phometh(-1.);
    for (int mask(0); mask < NumLegSubsets<int>(nlg + 1); ++mask) {
      double num(0.), den(0.);
      const int hgm(FlipLegs(HelicityBit<int>(hg), LegBitOf(mask, nlg)));
      for (int f = 0; f < nh; ++f) {
        const double C(std::abs(RealNorm()*ComixPhotonCoupling()*
                                sa[FlatIndex((size_t)FlipLegs(f, LowLegs(mask, nlg)),
                                             hgm, nlg)]));
        const double H(std::abs(soft*B.m_A[f]));
        num += sqr(H - C); den += H*H + C*C;
      }
      const double met(den > 0. ? sqrt(num/den) : 0.);
      if (mask == m_comixflip) inusemet = met;
      if (mask == FlipLegs(m_comixflip, LegBit<int>(4))) phometh = met;
      if (met < bmet) { bmet = met; bm = mask; }
    }
    msg_Error()<<"ACRAICSOFTFLIP lam="<<lam<<" inuse="<<m_comixflip
             <<" metric="<<inusemet<<" photonflipped="<<phometh
             <<" hel="<<hg
             <<" best="<<bm<<" bestmetric="<<bmet<<std::endl;
    // Raw Comix |M|^2 over the closed-form massive initial-state eikonal times
    // the Born (as NLO_Base::CheckRealCollinearSub): must be flat in lambda.
    const double mel(m_flavs[0].Mass());
    const double pk1(p1*k), pk2(p2*k), p12(p1*p2);
    const double Scl(pk1 != 0. && pk2 != 0.
                     ? 2.*p12/(pk1*pk2) - mel*mel/(pk1*pk1)
                                        - mel*mel/(pk2*pk2) : 0.);
    double nb(0.);
    for (int f = 0; f < nh; ++f) nb += std::norm(B.m_A[f]);
    msg_Error()<<"ACRAICEIKFLAT lam="<<lam
             <<" x="<<(2.*k[0]/sqrt(m_s))
             <<" ratio="<<((Scl != 0. && nb > 0.) ? sc/(Scl*nb) : -1.)
             <<std::endl;
    msg_Error()<<"ACRAICSOFTSCAN lam="<<lam
             <<" x="<<2.*k[0]/sqrt(dabs((p1+p2).Abs2()))
             <<" R="<<(shd > 0. ? sc/shd : -1.)
             <<" sum_comix="<<sc<<" sum_eikborn="<<shd
             <<" allhel/diff="<<(p_realprov->m_lastcomix != 0. ?
                                 scall/p_realprov->m_lastcomix : 0.)
             <<std::endl;
  }
  m_sp = sp_save;
  MakeProp(false);
  MakePropT(m_pacraic);
}


bool Acraic_Base::CalibrateComixMap(ComixCalib &c)
{
  if (p_bornproc == NULL || m_bornmomenta.size() < 4) return false;
  // The Born configuration, in its own rest frame.
  Vec4D_Vector bp(m_bornmomenta);
  Poincare bcms(bp[0] + bp[1]);
  for (size_t i(0); i < bp.size(); ++i) bcms.Boost(bp[i]);

  Amplitude cx, hand;
  const double sp_save(m_sp);
  m_sp = (bp[2] + bp[3]).Abs2();
  MakeProp(false);
  // MakePropT too: Comix's Born carries the t-channel at bp, so the hand Born
  // must as well (no-op except for Bhabha).
  MakePropT(bp);
  BornAmplitude(bp, hand);
  double cxme2(0.);
  const bool ok(ComixBornAmplitude(bp, cx, &cxme2));
  m_sp = sp_save;
  MakeProp(false);
  if (m_pacraic.size() >= 4) MakePropT(m_pacraic);
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
      met += sqr(std::abs(m_e*m_e*hand.m_A[f])*rh - std::abs(cx.m_A[FlipLegs(f, m)])*rc);
    met = sqrt(met);
    if (bestmask < 0 || met < bestmet) { nextmet = bestmet;
                                         bestmet = met; bestmask = m; }
    else if (nextmet < 0. || met < nextmet) nextmet = met;
  }

  // |C/H| must be the same in every live slot. The threshold is coarse on
  // purpose: Comix fills the helicity-flip entries at O(m_e/E), the hand Born
  // sets them to zero.
  double rmin(0.), rmax(0.);
  int nlive(0);
  for (int f = 0; f < nh; ++f) {
    const Complex H(m_e*m_e*hand.m_A[f]), C(cx.m_A[FlipLegs(f, bestmask)]);
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
  Adopt the Born calibration for COMIX_REAL_FLIP < 0 and/or
  COMIX_REAL_NORM <= 0. Once per run: the map is a convention, not per event.
*/
void Acraic_Base::DeriveComixMap()
{
  m_comixcalibrated = true;
  const bool needflip(m_comixflip < 0), neednorm(!(m_comixnorm > 0.));
  // a normalisation set by hand in the card is taken as given
  if (!neednorm) m_normexact = true;
  if (!needflip && !neednorm) return;

  ComixCalib c;
  if (!CalibrateComixMap(c)) {
    // Refuse rather than fall back to values fitted for one process.
    msg_Error()<<METHOD<<"(): Born calibration failed; the Comix real cannot"
               <<" be normalised. Set ACRAIC: COMIX_REAL_NORM and"
               <<" COMIX_REAL_FLIP by hand, or leave COMIX_REAL off.\n";
    m_comixreal = false;
    return;
  }

  const size_t nl(m_flavs.size());
  // The photon bit sits above the s_nlegs helicity slots, not m_flavs.size()
  // (they differ when a leg is a scalar, e.g. H l+ l-).
  if (needflip) m_comixflip = c.mask | (m_comixphoflip ? Amplitude::NHel() : 0);

  // N = sqrt(initial spin states) exactly: the hand-coded amplitude is
  // spin-averaged (KKMC convention), Comix's sum_hel is not. The measured N
  // only has to agree.
  double nspin(1.);
  for (size_t i(0); i < 2 && i < nl; ++i) nspin *= m_flavs[i].IntSpin() + 1.;
  // YFS: USE_MODEL_ALPHA 0: hand Born at alpha(0), Comix at alpha_model, so
  // m_rescale_alpha = alpha(0)/alpha_model enters (1 otherwise).
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

  msg_Debugging()<<METHOD<<"(): Comix -> ACRAIC map derived from the Born:\n"
            <<"  flip mask   = "<<m_comixflip<<"  (fermion bits "<<c.mask
            <<", metric "<<c.met<<" against "<<c.next<<" for the runner-up"
            <<(m_comixphoflip ? "; photon bit set" : "; photon bit clear")<<")\n"
            <<"  normalisation = "<<m_comixnorm<<"  (N = "<<c.N
            <<"; 1/sqrt(initial spin states) = "<<expected<<")\n"
            <<"  |C/H| over the "<<c.nlive<<" live helicities: ["
            <<c.rmin<<", "<<c.rmax<<"]\n";
}


//! The initial-state half of HandOnePhotonAmplitude, added to M1: the beam
//! spinors carry the photon index in turn.
void Acraic_Base::AddInitialStateArms(const Vec4D_Vector &bp, const Vec4D &k,
                                      int hel, Amplitude &M1)
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
  AddU(M1, AmpBornU, AmpU,  gI / (k*bp[0]) / 2.);
  AddV(M1, AmpBornV, AmpV, -gI / (k*bp[1]) / 2.);
}

//! The final-state half of HandOnePhotonAmplitude, added to M1: the same
//! construction on the outgoing pair.
void Acraic_Base::AddFinalStateArms(const Vec4D_Vector &bp, const Vec4D &k,
                                    int hel, Amplitude &M1)
{
  const double m3(m_flavs[2].Mass()), m4(m_flavs[3].Mass());
  const Vec4D_Vector arm1{bp[0], bp[1], k, bp[3]};
  const Vec4D_Vector arm2{bp[0], bp[1], bp[2], k};
  Amplitude AmpBornU, AmpBornV, AmpU, AmpV;
  BornAmplitude(arm1, AmpBornU, 0., m4, -1);
  BornAmplitude(arm2, AmpBornV, m3, 0., -1);
  UGamma(bp[2], k, k, hel, AmpU, m3, 0.);
  VGamma(k, bp[3], k, hel, AmpV, 0., m4);
  const double gF(m_qf * m_e * m_e * m_e);
  AddUF(M1, AmpU, AmpBornU,  gF / (k*bp[2]) / 2.);
  AddVF(M1, AmpBornV, AmpV, -gF / (k*bp[3]) / 2.);
}

/*!
  The hand-coded full one-photon amplitude (ISR + FSR) at a configuration:
  InfraredSubtractedME_1_0's spinor construction without the partition
  weights (sProd/Sactu) and the (1 - CKine) term, for comparison with Comix.
*/
bool Acraic_Base::HandOnePhotonAmplitude(const Vec4D_Vector &bp, const Vec4D &k,
                                       int hel, Amplitude &M1)
{
  if (bp.size() < 4) return false;
  const double p1k(k*bp[0]), p2k(k*bp[1]), p3k(k*bp[2]), p4k(k*bp[3]);
  if (IsZero(p1k) || IsZero(p2k) || IsZero(p3k) || IsZero(p4k)) return false;
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f) M1.m_A[f] = Complex(0., 0.);

  AddInitialStateArms(bp, k, hel, M1);
  AddFinalStateArms(bp, k, hel, M1);
  return true;
}


/*!
  Comix's one-photon amplitude at a momentum set taken verbatim, mapped into
  ACRAIC's convention. The caller may hold the fermions fixed and move only
  the photon; such a set needs COMIX: MOMENTUM_PROJECTION: 0.
*/
bool Acraic_Base::ComixRealAt(const Vec4D_Vector &pp, int hel, Amplitude &M1,
                            const double propscale)
{
  YFS::Real_Correction *prov(RealProvider(1));
  if (prov == NULL) return false;
  // Helicity legs (s_nlegs), not flavours: they differ for scalars (H l+ l-).
  const int nlg(Amplitude::s_nlegs);
  if (pp.size() != m_flavs.size() + 1) return false;
  const std::vector<METOOLS::Spin_Amplitudes> *amps
    (prov->ComixAmplitudes(pp, propscale));
  if (amps == NULL || amps->empty()) return false;
  const METOOLS::Spin_Amplitudes &sa((*amps)[0]);
  if ((int)sa.size() != NumLegSubsets<int>(nlg + 1)) return false;
  const int nh(Amplitude::NHel());
  const int hgi(FlipLegs(HelicityBit<int>(hel), ComixPhotonFlip(nlg)));
  for (int f = 0; f < nh; ++f)
    M1.m_A[f] = RealNorm() * ComixPhotonCoupling()
      * sa[FlatIndex((size_t)FlipLegs(f, ComixLegFlips(nlg)), hgi, nlg)];
  return true;
}

bool Acraic_Base::ComixRealShifted(const Vec4D_Vector &pp, int hel,
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
  if ((int)sa.size() != NumLegSubsets<int>(nlg + 1)) return false;
  const int nh(Amplitude::NHel());
  const int hgi(FlipLegs(HelicityBit<int>(hel), ComixPhotonFlip(nlg)));
  for (int f = 0; f < nh; ++f)
    M1.m_A[f] = RealNorm() * ComixPhotonCoupling()
      * sa[FlatIndex((size_t)FlipLegs(f, ComixLegFlips(nlg)), hgi, nlg)];
  return true;
}

Acraic_Base::PropShifts Acraic_Base::StageShifts(int iphot) const
{
  PropShifts sh;
  if (m_stage.size() != m_allphotons.size()) return sh;
  // One entry per stage: (Comix mask, its photons in this partition). A W
  // decay stage's mask is the W's daughters (StageShiftMask).
  for (size_t g(0); g < m_stagelegs.size(); ++g) {
    const size_t mask(StageShiftMask((int)g));
    if (mask == 0) continue;
    sh.push_back(std::make_pair(mask, StagePhotonSum((int)g, iphot)));
  }
  AddResonanceEntries(sh);   // ACRAIC: RESONANCE_EXPONENTIATION (else none)
  return sh;
}

/*
  ACRAIC: TCHANNEL_SHIFT. Automatic: on for 2 -> 2 only, since LegsAt's
  "radiating pair" beyond 2 -> 2 is just the first two final fermions and
  the shift then drags single-W t-channel photons towards t = 0.
*/
bool Acraic_Base::ExchangeLineShiftsOn() const
{
  static const tristate::code mode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                   ["TCHANNEL_SHIFT"].Get<tristate::code>());
  if (mode == tristate::automatic) return m_flavs.size() == 4;
  return mode == tristate::on;
}

/*
  The partition's reduced legs as leg entries theta_i (p~_i - p_i) (theta = -1
  incoming, +1 outgoing), from which COMIX::Amplitude::SetPropShifts places
  every space-like exchange line at the partition's reduced invariant; only
  the pole moves, the spinors stay physical. iphot < 0: the Born
  (BornLegsAt(m_PXvec)); iphot >= 0: M_1 of that photon (PartitionLegs, all
  zero at n = 1).
*/
void Acraic_Base::AddExchangeLineShifts(int iphot, PropShifts &sh) const
{
  if (!ExchangeLineShiftsOn() || m_pacraic.size() < 4
      || m_flavs.size() != m_pacraic.size()) return;
  Vec4D_Vector pb;
  const bool ok(iphot < 0 ? BornLegsAt(m_PXvec, pb) : PartitionLegs(iphot, pb));
  if (!ok || pb.size() != m_pacraic.size()) return;
  for (size_t i(0); i < m_pacraic.size(); ++i) {
    const Vec4D d((i < 2 ? -1. : 1.)*(pb[i] - m_pacraic[i]));
    sh.push_back(std::make_pair(LegUnion(LegBit(i),
                                PHASIC::Process_Base::s_propshiftleg), d));
  }
}

bool Acraic_Base::BornHasExchangeLine()
{
  std::map<const PHASIC::Process_Base*, int>::const_iterator
    it(m_exchline.find(p_bornproc));
  if (it != m_exchline.end()) return it->second != 0;
  if (p_bornproc == NULL || m_pacraic.size() < 4
      || m_flavs.size() != m_pacraic.size()) return false;
  // Probe: without an exchange line a leg shift leaves the amplitude bit for
  // bit unchanged.
  Amplitude A0, A1;
  PropShifts none, leg;
  const double e(1e-3*m_pacraic[0][0]);
  for (size_t i(0); i < m_pacraic.size(); ++i)
    leg.push_back(std::make_pair(LegUnion(LegBit(i),
                                 PHASIC::Process_Base::s_propshiftleg),
                                 Vec4D(e, 0.3*e*(i+1), -0.2*e, 0.5*e)));
  if (!ComixBornShifted(m_pacraic, A0, none) || !ComixBornShifted(m_pacraic, A1, leg))
    return false;                          // not cached: try again next event
  int has(0);
  for (int f(0); f < Amplitude::NHel(); ++f)
    if (A0.m_A[f] != A1.m_A[f]) { has = 1; break; }
  m_exchline[p_bornproc] = has;
  msg_Info()<<"ACRAIC: Born "<<(p_bornproc ? p_bornproc->Name() : std::string("?"))
            <<(has ? " has" : " has no")<<" space-like exchange line"
            <<(has ? "s" : "")<<"."<<std::endl;
  return has != 0;
}

Vec4D Acraic_Base::PartitionShift(int iphot, bool reducing) const
{
  Vec4D d;
  if (m_stage.size() != m_allphotons.size()) return d;
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if ((int)i != iphot
        && (m_stagereduces[m_stage[i]] != 0) == reducing) d += m_allphotons[i];
  return d;
}

/*!
  Comix's amplitude for one photon of the event, the others absorbed in the
  recoil. \param bp receives the reduced fermion configuration, so the
  hand-coded amplitude can be built at the same point.
*/
bool Acraic_Base::ComixOnePhotonAmplitude(const Vec4D &k, int hel,
                                        Vec4D_Vector &bp, Amplitude &M1,
                                        Vec4D &kmap)
{
  YFS::Real_Correction *prov(RealProvider(1));
  if (prov == NULL || prov->p_proc == NULL || m_pacraic.size() < 4) return false;
  const int nlg((int)m_flavs.size());
  if (nlg != 4) return false;   // the reduction below is the 2 -> 2 one

  // Q = q3 + q4 + k from the outgoing side (p1 + p2 - k would leave the other
  // photons inside the pair); the beams are rebuilt to match, as in
  // NLO_Base::MapMomenta.
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  const double m3(m_flavs[2].Mass()), m4(m_flavs[3].Mass());
  Vec4D q3(m_pacraic[2]), q4(m_pacraic[3]), kk(k);
  const Vec4D Q(q3 + q4 + kk);
  const double sQ(Q.Abs2());
  if (sQ <= sqr(m3 + m4)) return false;
  // Kaellen function
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
  if (LargestComponent(P1 + P2 - q3 - q4 - kk) > 1e-6*(sqrt(dabs(sQ)) + 1.))
    return false;

  bp.assign({P1, P2, q3, q4});
  Vec4D_Vector pp{P1, P2, q3, q4, kk};
  const std::vector<METOOLS::Spin_Amplitudes> *amps(prov->ComixAmplitudes(pp));
  if (amps == NULL || amps->empty()) return false;
  const METOOLS::Spin_Amplitudes &sa((*amps)[0]);
  if ((int)sa.size() != NumLegSubsets<int>(nlg + 1)) return false;

  const int nh(Amplitude::NHel());
  const int hgi(FlipLegs(HelicityBit<int>(hel), ComixPhotonFlip(nlg)));
  for (int f = 0; f < nh; ++f)
    M1.m_A[f] = RealNorm() * ComixPhotonCoupling()
      * sa[FlatIndex((size_t)FlipLegs(f, ComixLegFlips(nlg)), hgi, nlg)];
  kmap = kk;
  return true;
}



/*
  Closure test of the lambda subtraction, Comix only:
        s_implied(hel) = (M1(hel) - beta_1(hel)) / beta_0(hel)
  must be the same for every live helicity; its ratio to Sfactor is the
  convention constant.
*/
void Acraic_Base::Beta1Closure()
{
  if (m_allphotons.empty()) return;
  // Use the softest photon: beta_0 at the reduced point equals the unreduced
  // one only in the soft limit, and one-photon events are mostly hard.
  size_t js(0);
  for (size_t i(1); i < m_allphotons.size(); ++i)
    if (m_allphotons[i][0] < m_allphotons[js][0]) js = i;
  const int nh(Amplitude::NHel());
  const int msk(ComixLegFlips(Amplitude::s_nlegs));
  const double rn(RealNorm());
  if (!(rn > 0.)) return;
  Amplitude B0, M1, B1;
  if (!ComixBornAmplitude(m_pacraic, B0, NULL, m_svarQ)) return;
  Vec4D_Vector pp(m_pacraic); pp.push_back(m_allphotons[js]);
  if (!ComixRealAt(pp, m_PhoHel[js], M1, m_svarQ)) return;
  if (!ComixBeta1At(m_allphotons[js], m_PhoHel[js], m_svarQ, B1)) return;
  const Complex qratio(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
  const Complex sacraic(Sfactor(m_pacraic[0], m_pacraic[1],
                              m_allphotons[js], m_PhoHel[js])
                      + qratio*Sfactor(m_pacraic[m_if1], m_pacraic[m_if2],
                                       m_allphotons[js], m_PhoHel[js]));
  double nb(0.);
  for (int f(0); f < nh; ++f) nb += std::norm(B0.m_A[f]);
  nb = sqrt(nb);
  if (!(nb > 0.) || std::abs(sacraic) == 0.) return;
  int nlive(0); Complex first(0.,0.); double worst(0.);
  for (int f(0); f < nh; ++f) {
    if (std::abs(B0.m_A[f]) < 1e-3*nb) continue;
    const Complex m1((M1.m_A[FlipLegs(f, msk)])/rn), b1((B1.m_A[FlipLegs(f, msk)])/rn);
    const Complex si((m1 - b1)/B0.m_A[f]);
    if (nlive++ == 0) first = si;
    else if (std::abs(first) > 0.)
      worst = Max(worst, std::abs(si/first - 1.));
  }
  if (nlive < 2 || std::abs(first) == 0.) return;
  const double xg(m_s > 0. ? 2.*m_allphotons[js][0]/sqrt(m_s) : -1.);
  msg_Error()<<"B1CLOSE xg="<<xg<<" nlive="<<nlive
           <<" hel_spread="<<worst
           <<" s_implied_over_Sacraic="<<std::abs(first/sacraic)
           <<" arg="<<std::arg(first/sacraic)<<std::endl;
}


/*
  Normalisation of Sfactor: sum_hel |s(k,hel)|^2 must equal, up to e^2,

     Stilde(k) = - sum_ij Q_i Q_j th_i th_j (p_i.p_j)/((p_i.k)(p_j.k))
                 + sum_i  Q_i^2 m_i^2 / (p_i.k)^2
*/
void Acraic_Base::SoftNormCheck(const Vec4D &k)
{
  if (m_pacraic.size() < 4) return;
  // the charged legs with their YFS signs: incoming +, outgoing -
  std::vector<Vec4D>  pl;
  std::vector<double> Q, th;
  pl.push_back(m_pacraic[0]);     Q.push_back(m_qe); th.push_back(+1.);
  pl.push_back(m_pacraic[1]);     Q.push_back(-m_qe); th.push_back(+1.);
  pl.push_back(m_pacraic[m_if1]); Q.push_back(m_qf); th.push_back(-1.);
  pl.push_back(m_pacraic[m_if2]); Q.push_back(-m_qf); th.push_back(-1.);
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
  // sum over the two photon helicities of |s|^2, with ACRAIC's own Sfactor
  const Complex qratio(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
  double ss(0.);
  for (int h(-1); h <= 1; h += 2) {
    const Complex sa(Sfactor(m_pacraic[0], m_pacraic[1], k, h)
                     + qratio*Sfactor(m_pacraic[m_if1], m_pacraic[m_if2], k, h));
    ss += std::norm(sa);
  }
  if (st == 0.) return;
  msg_Error()<<"SOFTNORM xg="<<(m_s>0.?2.*k[0]/sqrt(m_s):-1.)
           <<" sum_hel|s|^2="<<ss<<" Stilde="<<st
           <<" ratio="<<(ss/st)<<std::endl;
}

/*
  Soft-limit test (hep-ph/0006359, eq.(single-initial)):

      | M_1(lambda k) - s_tot B | / | s_tot B |  ->  0   as lambda -> 0 ,

  on the beam-reduced legs (ReducedBeams) that beta_1 uses for multi-photon
  events. The angle between reduced and physical beam is printed, since a
  rotation of the beam axis does not leave the eikonal invariant.
*/
void Acraic_Base::SoftLimitTest()
{
  if (m_allphotons.empty() || !m_cxbalignok) return;
  const int nh(Amplitude::NHel()), flip(ComixLegFlips(Amplitude::s_nlegs));
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
      const Complex sb(stot*C.m_A[FlipLegs(f, flip)]);
      nb += std::norm(M1.m_A[f]/rn - sb); ne += std::norm(sb);
    }
    if (!(ne > 0.)) continue;
    const Vec3D a(pb[0]), b(m_pacraic[0]);
    const double ab(a.Abs()*b.Abs());
    const double ang(ab > 0. ? acos(Min(1., Max(-1., (a*b)/ab))) : -1.);
    out += " lam=" + ToString(lam) + " r=" + ToString(sqrt(nb/ne))
         + " ang=" + ToString(ang);
  }
  if (!out.empty())
    msg_Error()<<"SOFTLIM nphot="<<m_allphotons.size()<<out<<std::endl;
}

/*
  A real phase-space point with beams at total X and the radiating pair at
  total Y, each back to back on shell; directions from the event, spectators
  untouched. Comix needs a real point, not just a moved propagator pole.
*/
bool Acraic_Base::LegsAt(const Vec4D &X, const Vec4D &Y, Vec4D_Vector &pb) const
{
  pb = m_pacraic;
  const double X2(X.Abs2());
  const double ma(m_flavs[0].Mass()), mb(m_flavs[1].Mass());
  if (X2 <= sqr(ma + mb)) return false;
  Poincare cmsX(X);
  // ACRAIC: REDUCED_AXIS: beams along z in X's frame, as the generator
  // builds its s' point (matters for space-like exchange lines); false uses
  // the boosted-beam axis.
  static const bool reducedaxis(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                ["REDUCED_AXIS"].SetDefault(true).Get<bool>());
  Vec3D da(0., 0., 1.);
  if (!reducedaxis) {
    Vec4D ra(m_pacraic[0]);
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
  Vec4D r3(m_pacraic[m_if1]);
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

bool Acraic_Base::GeneratorBornAt(const Vec4D &R, Vec4D_Vector &pb) const
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

bool Acraic_Base::BornLegsAt(const Vec4D &X, Vec4D_Vector &pb) const
{
  Vec4D Y(X);
  for (size_t i(2); i < m_pacraic.size(); ++i)
    if (i != m_if1 && i != m_if2) Y -= m_pacraic[i];
  return LegsAt(X, Y, pb);
}

/*
  The legs for photon iphot in the current partition: the other photons the
  partition counts as initial-state come out of the beams (X_wp), the rest
  go into the radiating pair (X_wp - k). Physical point at one photon.
*/
bool Acraic_Base::PartitionLegs(int iphot, Vec4D_Vector &pb) const
{
  if (iphot < 0 || iphot >= (int)m_allphotons.size()
      || m_stage.size() != m_allphotons.size()) return false;
  Vec4D X(m_pacraic[0] + m_pacraic[1]);
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if ((int)i != iphot && m_stagereduces[m_stage[i]]) X -= m_allphotons[i];
  Vec4D Y(X - m_allphotons[iphot]);
  for (size_t i(2); i < m_pacraic.size(); ++i)
    if (i != m_if1 && i != m_if2) Y -= m_pacraic[i];
  return LegsAt(X, Y, pb);
}

tchmultiphoton::code Acraic_Base::TchannelMultiphotonMode() const
{
  static const tchmultiphoton::code mode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                         ["TCHANNEL_MULTIPHOTON"].Get<tchmultiphoton::code>());
  return mode;
}

/*
  ACRAIC: TCHANNEL_MULTIPHOTON >= 1. The one-photon point for photon iphot,
  as YFS.NLO's REAL_MAP 2 (NLO_Base::MapMomentaScaled): the event with only
  this photon, scaled by x = sqrt(Q^2/(P - k_j)^2) so the final state keeps
  Q^2, built in the P rest frame and carried back by the generator's boost
  Poincare(Q), so that its soft limit is exactly BornLegsAt(Q). Used instead
  of PartitionLegs, whose beams are tilted by a hard companion's pT.
  legs gets the Born legs, kpt the photon.
*/
bool Acraic_Base::OnePhotonScaledLegs(int iphot, Vec4D_Vector &legs,
                                    Vec4D &kpt) const
{
  if (iphot < 0 || iphot >= (int)m_allphotons.size()
      || m_pacraic.size() != m_flavs.size() || m_pacraic.size() < 4) return false;
  const Vec4D P(m_pacraic[0] + m_pacraic[1]);
  Vec4D Q;
  for (size_t i(2); i < m_pacraic.size(); ++i) Q += m_pacraic[i];
  const Vec4D R(P - m_allphotons[iphot]);
  const double P2(P.Abs2()), Q2(Q.Abs2()), R2(R.Abs2());
  if (!(P2 > 0.) || !(Q2 > 0.) || !(R2 > 0.) || !(Q[0] > 0.)) return false;
  const double x(Min(1., sqrt(Q2/R2)));
  const double ma(m_flavs[0].Mass()), mb(m_flavs[1].Mass());
  const double sx(x*x*P2);
  const double kal(sqr(sx - ma*ma - mb*mb) - 4.*ma*ma*mb*mb);
  if (!(sx > sqr(ma + mb)) || !(kal > 0.)) return false;
  // in the P rest frame: beams back to back along z at x sqrt(P^2), the
  // photon scaled by x, the final state carried from Q to Q' = x(P - k)
  Poincare toP(P), toQ(Q);
  Vec4D k(m_allphotons[iphot]);
  toP.Boost(k);
  k *= x;
  // leg 0 along +z, as LegsAt's reduced axis (REDUCED_AXIS 1) puts it
  const double pz(sqrt(kal)/(2.*sqrt(sx)));
  const Vec4D pa(sqrt(pz*pz + ma*ma), 0., 0.,  pz);
  const Vec4D pb(sqrt(pz*pz + mb*mb), 0., 0., -pz);
  const Vec4D Qp(pa + pb - k);
  if (!(Qp[0] > 0.) || !(Qp.Abs2() > 0.)) return false;
  Poincare toQp(Qp);
  legs = m_pacraic;
  legs[0] = pa; legs[1] = pb;
  for (size_t i(2); i < m_pacraic.size(); ++i) {
    Vec4D q(m_pacraic[i]);
    toQ.Boost(q);           // the generator's Born frame
    toQp.BoostBack(q);      // to Q' in the P rest frame
    legs[i] = q;
  }
  // back by the generator's boost: the soft limit is BornLegsAt(Q)
  for (Vec4D &l : legs) toQ.BoostBack(l);
  toQ.BoostBack(k);
  kpt = k;
  Vec4D bal(legs[0] + legs[1] - kpt);
  for (size_t i(2); i < legs.size(); ++i) bal -= legs[i];
  const double scale(Max(legs[0][0] + legs[1][0], 1.));
  for (int mu(0); mu < 4; ++mu)
    if (dabs(bal[mu]) > 1e-10*scale) return false;
  return true;
}

/*
  ACRAIC: TCHANNEL_MULTIPHOTON 2. Factorised form of the partition's beta_1
  sum a + sum_j v_j (a = m_Sprod beta_0), the amplitude-level analogue of
  YFS.NLO's REAL_COMBINE product, for Borns near a pole of the reduced frame:
      c_j = <a|a + v_j>/<a|a>,   r_j = v_j - (c_j - 1) a,
      T   = a prod_j c_j + sum_j r_j prod_{l != j} c_l .
  T = a + sum_j v_j to first order (exactly for one photon); the difference
  T - a - sum_j v_j is added.
*/
void Acraic_Base::AddFactorisedRemainder()
{
  const size_t nt(m_b1terms.size());
  if (nt < 2) return;
  const int nh(Amplitude::NHel());
  std::vector<Complex> a(nh);
  double na(0.);
  for (int f(0); f < nh; ++f) {
    a[f] = m_Sprod * m_partborn0.m_A[f];
    na += std::norm(a[f]);
  }
  if (!(na > 0.) || IsBad(na)) return;
  std::vector<Complex> c(nt, Complex(1., 0.));
  for (size_t j(0); j < nt; ++j) {
    Complex proj(0., 0.);
    for (int f(0); f < nh; ++f) proj += std::conj(a[f]) * m_b1terms[j].m_A[f];
    c[j] += proj/na;
  }
  auto product_without = [&](size_t skip) {
    Complex p(1., 0.);
    for (size_t l(0); l < nt; ++l) if (l != skip) p *= c[l];
    return p; };
  const Complex call(product_without(nt));
  for (int f(0); f < nh; ++f) {
    Complex additive(a[f]), factorised(a[f] * call);
    for (size_t j(0); j < nt; ++j) {
      const Complex vj(m_b1terms[j].m_A[f]);
      additive += vj;
      factorised += (vj - (c[j] - 1.) * a[f]) * product_without(j);
    }
    const Complex d(factorised - additive);
    if (IsBad(std::abs(d))) continue;
    m_AmpExpo1.m_A[f]    += d;
    m_AmpBornReal.m_A[f] += d;
  }
}

/*
  The legs for photon k when the event has other photons: Comix does not
  check momentum balance, so the beams are rebuilt back to back on shell at
  (final legs + k)^2, keeping the final legs and k. The missing momentum goes
  into the beams, not the pair, so beta_0 is unchanged and the eikonal nearly
  so (p/(p.k) is invariant under p -> x p). Physical beams at one photon.
*/
bool Acraic_Base::ReducedBeams(const Vec4D &k, Vec4D_Vector &pb) const
{
  pb = m_pacraic;
  Vec4D S(k);
  for (size_t i(2); i < m_pacraic.size(); ++i) S += m_pacraic[i];
  const double S2(S.Abs2());
  const double ma(m_flavs[0].Mass()), mb(m_flavs[1].Mass());
  if (S2 <= sqr(ma + mb)) return false;
  Poincare cms(S);
  Vec4D pa(m_pacraic[0]);
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
  Total eikonal sum_l Q_l theta_l SfactorLeg(p_l) (eq. 6.4); equals
  sum_g m_Sfac[g][j] on the physical legs, recomputed here on reduced legs.
*/
Complex Acraic_Base::TotalEikonal(const Vec4D_Vector &p, const Vec4D &k,
                                int hel)
{
  Complex tot(0., 0.);
  for (size_t l(0); l < m_flavs.size() && l < p.size(); ++l) {
    // ISR only: the final current does not radiate (see m_isronly)
    if (m_isronly && l >= 2) continue;
    const double w(m_flavs[l].Charge() * (l < 2 ? -1. : +1.));
    if (w != 0.) tot += w * SfactorLeg(p[l], k, hel);
  }
  return tot;
}

Complex Acraic_Base::StageEikonal(int stage, const Vec4D_Vector &p,
                                const Vec4D &k, int hel, int iphot)
{
  // every leg of the stage, a resonance built from its daughters in p
  return StageCurrent(stage, iphot, k, hel, p);
}

/*
  beta_1 for one photon from Comix, partition independent:

      beta_1(k) = M_1(k) - s_tot(k) B ,

  s_tot the eikonal of every leg, B the Born on the same legs; no stage
  split, pseudo-flux or (1 - CKine) needed, since

      M_n  ~  sum_wp prod s^wp beta_0(X_wp)  +  sum_j prod_{i!=j} s_tot(k_i)
                                                * beta_1(k_j).

  Adding it on one designated stage of photon j visits each assignment of
  the others once: sum_{wp: wp_j = 0} prod_{i!=j} s^{wp_i} = prod_{i!=j} s_tot(k_i).
*/
bool Acraic_Base::ComixInfraredSubtracted_1_0(const Vec4D &k, int hel,
                                            int iphot)
{
  if (!m_cxbalignok || iphot < 0) return false;
  if (iphot >= (int)m_stage.size() || (int)m_Sfac.size() != m_nstages
      || iphot >= (int)m_Sfac[m_stage[iphot]].size()) return false;
  /*
    M_1 goes on one partition per assignment of the others (photon on the
    reducing stage, never "index 0", which is the non-radiating stage for
    nu nu); the beta_0 subtraction goes in every partition. A collapsed
    (fixed-stage) photon is never enumerated, so it gets M_1 in every
    partition against the total eikonal.
  */
  const bool fixed(m_fixedstage.size() == m_stage.size()
                   && m_fixedstage[iphot] != 0);
  // YFS: IFI false: partitions summed incoherently, each carries its own
  // share of M_1 (AddPlacementVetoes).
  const bool addm1(!m_ifi || fixed || m_stagereduces[m_stage[iphot]]);
  const double rn(RealNorm());
  if (!(rn > 0.)) return false;
  /*
    Per assignment of the others, the two partitions differing only in this
    photon's placement give (hep-ph/0006359)

      prod_{i!=j} s_i^{wp_i} [ M_1^I(k_j; X_{j in I}) + M_1^F(k_j; X_{j in F})
                              - s_I(k_j) beta_0(X_{j in I})
                              - s_F(k_j) beta_0(X_{j in F}) ],

    M_1^{I,F} at the physical spinors with propagators at the partition's X.
    ACRAIC: BETA1_LEGS: 0 (default): physical legs + k_j with per-stage
    propagator shifts (ComixRealShifted). 1: balanced point from
    PartitionLegs (right propagators, wrong spinors), kept for comparison.
  */
  static const beta1legs::code legsmode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                            ["BETA1_LEGS"].Get<beta1legs::code>());
  const bool rebuild(legsmode == beta1legs::balanced || m_redborn);
  Vec4D_Vector pb;
  Vec4D dI;
  PropShifts shifts;
  if (rebuild) { if (!PartitionLegs(iphot, pb)) return false; }
  else {
    pb = m_pacraic;
    shifts = StageShifts(iphot);
    AddExchangeLineShifts(iphot, shifts); // space-like exchange lines
    // YFS: IFI false, or ISR only (every photon on the production stage)
    if (!m_ifi || m_isronly) AddPlacementVetoes(std::vector<int>(1, iphot), shifts);
    dI = PartitionShift(iphot, true);   // the initial stage's photons, for S
  }
  // Softness y = 2k.S/S^2 against the radiating system S (beams less the
  // other initial-state photons); independent of this photon's placement, so
  // both partitions of one beta_1 term are cut together. Below ACRAIC:
  // BETA1_XCUT beta_1 = 0: no M_1, no subtraction.
  Vec4D S(rebuild ? pb[0] + pb[1] : m_pacraic[0] + m_pacraic[1] - dI);
  // A photon on a W decay stage: softness against P_W as produced
  // (daughters + the stage's photons).
  if (WStagesActive() && !rebuild && m_stage[iphot] != m_initstage
      && !fixed) {
    const int g(m_stage[iphot]);
    for (size_t l(0); l < m_stagelegs[g].size(); ++l)
      if (IsResonanceLeg(m_stagelegs[g][l].leg))
        S = StageSystemMomentum(g) + StagePhotonSum(g, -1);
  }
  const double S2(S.Abs2());
  const double y(S2 > 0. ? 2.*(k*S)/S2 : 0.);
  static const double xcut(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                           ["BETA1_XCUT"].Get<double>());
  if (!(y > xcut)) return true;     // m_b2on[iphot] stays 0: beta_1 = 0 here
  /*
    The subtraction is m_partborn0 (the Born InfraredSubtractedME_0_0 added,
    so the (1 - n) cancellation is exact) times the other photons' eikonals
    and this photon's stage eikonal on the legs Comix was handed.

    ACRAIC: TCHANNEL_MULTIPHOTON >= 1: M_1 at OnePhotonScaledLegs, used as
    R_j = M_1/s_j(point) times the physical m_sactu[j] (rescale), so the
    subtraction is exactly s_phys B_0.
  */
  Vec4D km(k);
  bool onept(false);
  if (addm1 && rebuild && TchannelMultiphotonActive()
      && m_sactu.size() == m_allphotons.size()) {
    Vec4D_Vector po;
    Vec4D ko;
    if (OnePhotonScaledLegs(iphot, po, ko)) { pb = po; km = ko; onept = true; }
    else {
      ++m_tmpfail;
      long d(1);
      while (d*10 <= m_tmpfail) d *= 10;
      if (m_tmpfail == d)
        msg_Error()<<METHOD<<"(): ACRAIC: TCHANNEL_MULTIPHOTON: no one-photon "
                   <<"point for "<<m_tmpfail<<" photon(s) so far; those keep "
                   <<"PartitionLegs."<<std::endl;
    }
  }
  Amplitude M1;
  if (addm1) {
    Vec4D_Vector pp(pb);
    pp.push_back(km);
    const bool ok(rebuild ? ComixRealAt(pp, hel, M1, -1.)
                          : ComixRealShifted(pp, hel, M1, shifts));
    if (!ok) return false;
  }
  // the other photons' eikonal product; never sProd/s_j (diverges at zeros of s_j)
  Complex w(1., 0.);
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if ((int)i != iphot)
      w *= (m_sactu.size() == m_allphotons.size() ? m_sactu[i] : m_Sfac[m_stage[i]][i]);
  Complex sj(fixed ? TotalEikonal(pb, km, hel)
                   : StageEikonal(m_stage[iphot], pb, km, hel, iphot));
  Complex rescale(1., 0.);
  if (onept) {
    if (!(std::abs(sj) > 0.)) return true;   // no eikonal at the point
    rescale = m_sactu[iphot]/sj;
    sj = m_sactu[iphot];
  }
  const int nh(Amplitude::NHel());
  // ACRAIC: ORDER 2: store this beta_1's M_1 and eikonal for beta_2's
  // subtraction (Acraic_Beta2.C).
  if (m_order == 2 && iphot < (int)m_b2on.size()) {
    m_b2on[iphot] = 1;
    m_b2sj[iphot] = sj;
    m_b2hasM1[iphot] = addm1 ? 1 : 0;
    if (addm1)
      for (int f(0); f < nh; ++f)
        m_b2M1[iphot].m_A[f] = rescale * m_cxbalign.m_A[f] * M1.m_A[f]/rn;
  }
  double nsub(0.), nm1(0.), nv(0.);
  // TCHANNEL_MULTIPHOTON 2: this photon's term, for AddFactorisedRemainder
  const bool keepterm(onept && TchannelMultiphotonMode() == tchmultiphoton::factorised);
  if (keepterm) m_b1terms.push_back(Amplitude());
  for (int f(0); f < nh; ++f) {
    const Complex sub(w * sj * m_partborn0.m_A[f]);
    Complex v(-sub);
    if (addm1) {
      const Complex m1(rescale * w * m_cxbalign.m_A[f] * M1.m_A[f]/rn);
      v += m1;
      nm1 += std::norm(m1);
      if (iphot < (int)m_realphotM1.size()) m_realphotM1[iphot].m_A[f] += m1;
    }
    nsub += std::norm(sub); nv += std::norm(v);
    m_AmpExpo1.m_A[f]    += v;
    m_AmpBornReal.m_A[f] += v;
    if (iphot < (int)m_realphot.size()) m_realphot[iphot].m_A[f] += v;
    m_snapReal.m_A[f]    += v;
    m_beta10 += v;
    if (keepterm) m_b1terms.back().m_A[f] = v;
  }
  if (m_b1trace) {
    std::string st;
    for (size_t i(0); i < m_stage.size(); ++i) st += ToString(m_stage[i]);
    const Complex sp(m_Sfac[m_stage[iphot]][iphot]);
    // rsoft: M_1 vs total eikonal x Comix Born on the same legs;
    // bratio: partition Born vs that Comix Born.
    auto born_here = [&](const Vec4D_Vector &legs, Amplitude &C) {
      return rebuild ? ComixBornAmplitude(legs, C, NULL, -1., -1.)
                     : ComixBornShifted(legs, C, shifts); };
    auto real_here = [&](const Vec4D_Vector &legs, const Vec4D &kk,
                         Amplitude &M) {
      Vec4D_Vector pp(legs); pp.push_back(kk);
      return rebuild ? ComixRealAt(pp, hel, M, -1.)
                     : ComixRealShifted(pp, hel, M, shifts); };
    const int flip(ComixLegFlips(Amplitude::s_nlegs));
    double rsoft(-1.), bratio(-1.);
    Amplitude Cc;
    if (born_here(pb, Cc)) {
      const Complex stot(TotalEikonal(pb, km, hel));
      double nd(0.), ne(0.), nc(0.), nb(0.);
      for (int f(0); f < nh; ++f) {
        const Complex sb(stot*Cc.m_A[FlipLegs(f, flip)]);
        if (addm1) { nd += std::norm(M1.m_A[f]/rn - sb); ne += std::norm(sb); }
        nc += std::norm(m_cxbalign.m_A[f]*Cc.m_A[FlipLegs(f, flip)]);
        nb += std::norm(m_partborn0.m_A[f]);
      }
      if (addm1 && ne > 0.) rsoft = sqrt(nd/ne);
      if (nc > 0.) bratio = sqrt(nb/nc);
    }
    // The lambda scan needs balancing legs; with shifts only lambda = 1.
    if (addm1) {
      std::string sc;
      for (int e(0); e <= (rebuild ? 3 : 0); ++e) {
        const double lam(pow(10., -(double)e));
        const Vec4D kl(lam*km);
        Vec4D_Vector pl(pb);
        if (rebuild) {
          Vec4D X(pb[0] + pb[1]), Y(X - kl);
          for (size_t i(2); i < m_pacraic.size(); ++i)
            if (i != m_if1 && i != m_if2) Y -= m_pacraic[i];
          if (!LegsAt(X, Y, pl)) continue;
        }
        Amplitude Cl, Ml;
        if (!born_here(pl, Cl) || !real_here(pl, kl, Ml)) continue;
        const Complex stl(TotalEikonal(pl, kl, hel));
        double nd(0.), ne(0.), nm(0.);
        Complex ov(0., 0.);
        for (int f(0); f < nh; ++f) {
          const Complex sb(stl*Cl.m_A[FlipLegs(f, flip)]);
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
             <<" onept="<<onept<<" rescale="<<rescale
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

bool Acraic_Base::ComixBeta1At(const Vec4D &k, int hel,
                             const double propscale, Amplitude &B1)
{
  const double xg(m_s > 0. ? 2.*k[0]/sqrt(m_s) : 0.);
  static const double X_EIK(1e-4);
  // Already eikonal: beta_1 vanishes and the difference is below precision.
  if (!(xg > 10.*X_EIK)) {
    for (int f(0); f < Amplitude::NHel(); ++f) B1.m_A[f] = Complex(0., 0.);
    return true;
  }
  /*
    beta_1(k) = M_1(k) - s(k) beta_0, from the definition. B1 is returned
    Comix-ordered: hence the mask on beta_0 and RealNorm on the eikonal term.
  */
  Amplitude B0, M1;
  Vec4D_Vector pp(m_pacraic);
  pp.push_back(k);
  // Subtract the same Born the partition sum adds: with ACRAIC:
  // BORN_AT_SPRIME (or m_redborn) that is the reduced-leg Born (BornLegsAt).
  static const bool sprime(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                           ["BORN_AT_SPRIME"].Get<bool>());
  bool ok(false);
  if (sprime || m_redborn) {
    Vec4D_Vector pb;
    ok = BornLegsAt(m_PXvec, pb) && ComixBornAmplitude(pb, B0, NULL, -1., -1.);
  }
  if (!ok && !ComixBornAmplitude(m_pacraic, B0, NULL, propscale)) return false;
  if (!ComixRealAt(pp, hel, M1, propscale)) return false;
  const double rn(RealNorm());
  if (!(rn > 0.)) return false;
  const Complex qratio(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
  const Complex sf(Sfactor(m_pacraic[0], m_pacraic[1], k, hel)
                   + qratio*Sfactor(m_pacraic[m_if1], m_pacraic[m_if2], k, hel));
  // the paper's pseudo-flux X^2/(p_c+p_d)^2; exactly 1 for a single ISR photon
  const double pflux(propscale > 0. && m_svarQ > 0. ? propscale/m_svarQ : 1.);
  const int nh(Amplitude::NHel()), msk(ComixLegFlips(Amplitude::s_nlegs));
  for (int f(0); f < nh; ++f)
    B1.m_A[f] = M1.m_A[f] - rn*sf*pflux*B0.m_A[FlipLegs(f, msk)];
  return true;
}

void Acraic_Base::BuildComixPhotonRatios()
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
  // lambda is chosen per photon so the scaled photon sits at x = X_EIK: deep
  // enough to be eikonal, not so deep that spinor products lose precision.
  // Photons already below 10 X_EIK keep the hand-coded amplitude.
  static const double X_EIK(1e-4);

  for (size_t j(0); j < ng; ++j) {
    const Vec4D k0(m_allphotons[j]);
    const int hel(m_PhoHel[j]);

    // Comix at this photon and at the same photon scaled down; phase
    // comparable since Amplitude::SetGauge uses a fixed reference vector.
    Vec4D_Vector bp, bps;
    Vec4D kmap, kmaps;
    Amplitude C, Cs;
    // ACRAIC: COMIX_REAL_ACRAIC_ARGS: ACRAIC's own (non-conserving) spinor
    // arguments, unmapped; needs COMIX: MOMENTUM_PROJECTION: 0. Mapping to a
    // conserving point moves the Born by O(E_gamma), the order of beta_1.
    static const bool acraicargs(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                               ["COMIX_REAL_ACRAIC_ARGS"].Get<bool>());
    if (acraicargs) {
      bp.assign(m_pacraic.begin(), m_pacraic.begin()+4);
      kmap = k0;
      Vec4D_Vector ppc{bp[0], bp[1], bp[2], bp[3], kmap};
      if (!ComixRealAt(ppc, hel, C, m_svarQ)) { ++m_cxppfail; continue; }
    } else
    if (!ComixOnePhotonAmplitude(k0, hel, bp, C, kmap)) { ++m_cxppfail; continue; }
    // Same fermions, only the photon scaled: the Born must stay put. The set
    // is short by (1-lambda) k, accepted with MOMENTUM_PROJECTION off.
    const double xg(2.*kmap[0]/sqrt(m_s));
    if (!(xg > 10.*X_EIK)) { ++m_cxppfail; continue; }  // already eikonal
    const double lam(X_EIK/xg);
    kmaps = lam*kmap;
    bps   = bp;
    Vec4D_Vector pps{bp[0], bp[1], bp[2], bp[3], kmaps};
    if (!ComixRealAt(pps, hel, Cs, acraicargs ? m_svarQ : -1.))
      { ++m_cxppfail; continue; }

    // The eikonal ratio: the polarisation phase depends only on the photon's
    // direction, so it cancels under k -> lambda k.
    const Complex qratio(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
    const Complex Sk (Sfactor(bp[0],  bp[1],  kmap,  hel)
                      + qratio*Sfactor(bp[2],  bp[3],  kmap,  hel));
    const Complex Sks(Sfactor(bps[0], bps[1], kmaps, hel)
                      + qratio*Sfactor(bps[2], bps[3], kmaps, hel));
    if (std::abs(Sks) == 0.) { ++m_cxppfail; continue; }
    const Complex srat(Sk/Sks);

    Amplitude Csub;   // filled below, once Comix's own Born is in hand

    // Per-helicity map to ACRAIC's convention, a_f = e^2 B_hand(f)/B_comix(f):
    // modulus COMIX_REAL_NORM, phase the helicity-dependent spinor convention.
    Amplitude Bc, Bh, H;
    bool okb(false), okh(false);
    if (acraicargs) {
      // The Born at the same (non-conserving) arguments and scale.
      m_sp = m_svarQ;
      MakeProp(false);
      MakePropT(bp);
      okb = ComixBornAmplitude(bp, Bc, NULL, m_svarQ);
      BornAmplitude(bp, Bh, -1., -1., -1);
      okh = HandOnePhotonAmplitude(bp, kmap, hel, H);
      m_sp = sp_save;
      MakeProp(false);
      MakePropT(m_pacraic);
    } else {
      // Mapped route: the Born needs its own conserving 2 -> 2 point (Comix
      // returns NaN on the 2 -> 3 fermion legs alone).
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
      Poincare bq2(Q2);
      bq2.BoostBack(B1);
      bq2.BoostBack(B2);
      Vec4D_Vector b2{B1, B2, bp[2], bp[3]};
      m_sp = s2;
      MakeProp(false);
      MakePropT(b2);
      okb = ComixBornAmplitude(b2, Bc);
      BornAmplitude(b2, Bh, -1., -1., -1);
      okh = HandOnePhotonAmplitude(bp, kmap, hel, H);
      m_sp = sp_save;
      MakeProp(false);
      MakePropT(m_pacraic);
    }
    if (!okb || !okh) { ++m_cxppfail; continue; }

    /*
      beta_1 = M_1(k) - s(k) M_0 with the eikonal term from Comix itself:
      s(lambda k) = s(k)/lambda exactly, so

          s(k) M_0  =  lambda M_1(lambda k)   as lambda -> 0,

      and no gauge reference enters (ACRAIC's m_zeta, m_eta, m_b differ from
      Amplitude::SetGauge's by a little-group phase).
    */
    for (int f = 0; f < nh; ++f)
      Csub.m_A[f] = C.m_A[f] - lam*Cs.m_A[f];
    m_cxcsub[j] = Csub;
    double nrmH(0.), num(0.), den(0.);
    for (int f = 0; f < nh; ++f) nrmH += std::norm(H.m_A[f]);
    nrmH = sqrt(nrmH);
    int nal(0);
    for (int f = 0; f < nh; ++f) {
      const Complex Bcf(Bc.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]);
      const Complex Bhf(m_e*m_e*Bh.m_A[f]);
      if (std::abs(Bcf) == 0. || std::abs(H.m_A[f]) <= 1e-8*nrmH) {
        m_cxratio[j].m_A[f] = Complex(1., 0.);
        continue;
      }
      const Complex a(Bhf/Bcf);                 // Comix -> ACRAIC, per helicity
      const Complex cs(a*Csub.m_A[f]);
      m_cxratio[j].m_A[f] = cs/H.m_A[f];
      num += std::norm(cs - H.m_A[f]);
      den += std::norm(cs) + std::norm(H.m_A[f]);
      ++nal;
    }
    Beta1Check(j, nh, C, Csub, H, Bc, Sk);
    if (nal == 0) { ++m_cxppfail; continue; }
    m_cxratiook[j] = 1;
    if (den > 0.) m_cxppdev += sqrt(num/den);
    ++m_cxppn;
  }
  m_sp = sp_save;
  MakeProp(false);
  MakePropT(m_pacraic);
}


/*!
  Real-subtraction validation scan of beta_1, as NLO_Base::CheckRealSub.
*/
void Acraic_Base::Beta1Scan()
{
  static bool done(false);
  if (done || m_allphotons.empty() || m_pacraic.size() < 4) return;
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
  std::ofstream fc(("AcraicBeta1_comix" + tag + ".txt").c_str());
  std::ofstream fh(("AcraicBeta1_hand"  + tag + ".txt").c_str());
  fc << std::setprecision(10);
  fh << std::setprecision(10);

  Vec4D_Vector bp(m_pacraic.begin(), m_pacraic.begin()+4);
  const int nh(Amplitude::NHel());
  const int hel(m_PhoHel.empty() ? 1 : m_PhoHel[0]);
  const double sp_save(m_sp);

  // the Born both curves are normalised to, once
  m_sp = m_svarQ; MakeProp(false); MakePropT(bp);
  Amplitude Bc, Bh;
  const bool okb(ComixBornAmplitude(bp, Bc, NULL, m_svarQ));
  BornAmplitude(bp, Bh, -1., -1., -1);
  m_sp = sp_save; MakeProp(false); MakePropT(m_pacraic);
  if (!okb) return;
  double nbc(0.), nbh(0.);
  for (int f = 0; f < nh; ++f) {
    nbc += std::norm(Bc.m_A[f]);
    nbh += std::norm(m_e*m_e*Bh.m_A[f]);
  }
  if (!(nbc > 0.) || !(nbh > 0.)) return;
  // Expect |B_comix| = 2 e^2 |B_hand| (sqrt of the initial spin states).
  msg_Debugging()<<"ACRAIC: beta_1 scan Born check: |B_comix|/(2 e^2 |B_hand|) = "
            <<(sqrt(nbc)/(2.*sqrt(nbh)))<<"  (e^2 |B_hand| = "<<sqrt(nbh)
            <<", |B_comix| = "<<sqrt(nbc)<<")\n";
  // Both curves normalised to the same (Comix) Born.
  const double nbref(nbc);

  // Logarithmic scan: ACRAIC: BETA1_SCAN_DECADES decades, BETA1_SCAN_PPD
  // points per decade.
  static const int ndec(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                        ["BETA1_SCAN_DECADES"].SetDefault(8).Get<int>());
  static const int nppd(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                        ["BETA1_SCAN_PPD"].SetDefault(4).Get<int>());
  for (int it(0); it <= ndec*nppd; ++it) {
    const Vec4D k(k0*pow(10., -double(it)/double(nppd)));
    const double xg(2.*k[0]/sqrt(m_s));
    if (xg <= 0.) break;

    // Plot |beta_1(k)|/|s(k)|, which vanishes linearly (beta_1 alone tends
    // to a constant); moduli only, so the polarisation convention drops out.
    const Complex qrat(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
    const Complex Sk(Sfactor(bp[0], bp[1], k, hel)
                     + qrat*Sfactor(bp[2], bp[3], k, hel));
    const double aS(std::abs(Sk));
    if (!(aS > 0.)) continue;

    // Comix: beta_1 = M_1(k) - lambda M_1(lambda k), same fermions and scale
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

    // the hand-coded beta_1 at the same point
    m_sp = m_svarQ; MakeProp(false); MakePropT(bp);
    Amplitude H;
    const bool okh(HandOnePhotonAmplitude(bp, k, hel, H));
    m_sp = sp_save; MakeProp(false); MakePropT(m_pacraic);
    if (okh) {
      double nhh(0.);
      for (int f = 0; f < nh; ++f) nhh += std::norm(H.m_A[f]);
      fh << k[0] << "," << sqrt(nhh/nbref)/aS << std::endl;
    }
  }
  fc.close(); fh.close();
  msg_Debugging()<<"ACRAIC: beta_1 scan written to AcraicBeta1_{comix,hand}"<<tag
            <<".txt\n";
}
