/*!
  \file Ceex_Beta2.C

  beta_2: the O(alpha^2) tree-level double real of CEEX (CEEX: ORDER: 2),
  built from Comix's two-photon spin amplitudes.

  hep-ph/0006359 eq. (95), per partition wp and pair (j,l):

    beta_2^{w_j w_l}(k_j,k_l; X_w) = M_2^{w_j w_l}
          - beta_1^{w_j}(k_j; X_w) s^{w_l}(k_l)
          - beta_1^{w_l}(k_l; X_w) s^{w_j}(k_j)
          - beta_0(X_w) s^{w_j}(k_j) s^{w_l}(k_l),

  multiplied by the other photons' soft factors. KKMC hand-codes the three
  classes (GPS_HiiPlus, GPS_HffPlus, GPS_HifPlus). Here the whole bracket
  comes from amplitudes the partition sum already has:

    - M_2 is ONE Comix two-photon amplitude on the physical legs plus k_j and
      k_l, with the other photons as per-stage propagator shifts
      (StageShifts2). With both photons explicit, a current holding a stage
      gets that stage's other photons, so the one amplitude is the sum over
      the placements of j and l of M_2^{w_j w_l}(X_w), each diagram at its
      own pole - the same "whole amplitude once per assignment of the others"
      rule the Comix beta_1 uses for M_1 (ComixInfraredSubtracted_1_0). It is
      therefore added only in the partition where BOTH photons sit on the
      reducing (initial) stage.
    - beta_1(k_j) is exactly what the beta_1 pass put into THIS partition:
      its M_1 (if this partition carried it) minus its stage eikonal times
      m_partborn0, or nothing if beta_1 was cut. The pass records all three
      (m_b2M1, m_b2sj, m_b2on) so no M_1 is evaluated twice and, as for
      beta_1 against beta_0, the subtractions are the same numbers as the
      additions.
    - beta_0(X_w) is m_partborn0, and s^{w} are the partition's soft factors
      m_sactu (the ones beta_0 and the beta_1 weights carry).

  Summed over the partitions of two photons this gives M_2 identically: the
  n = 2 closure (Beta2Closure). What it tests is the bookkeeping; the physics
  test is the soft limits (Beta2SoftLimit) and the n = 3 KKMC point.

  Replaces the hand-coded beta_2 that stood here (never wired in, and built
  from the GPS U/V matrices, whose photon polarisation is not Comix's).
*/

#include "YFS/CEEX/Ceex_Base.H"
#include "YFS/NLO/Real_Correction.H"

#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"

#include <iostream>
#include <iomanip>

using namespace YFS;
using namespace ATOOLS;

namespace {
  inline size_t FlatIndex2(size_t ferm, int hg, int nlegs)
  { return ferm | ((size_t)hg << nlegs); }
}


void Ceex_Base::ResetBeta1Record()
{
  const size_t n(m_allphotons.size());
  if (m_b2M1.size() != n) m_b2M1.assign(n, Amplitude());
  m_b2hasM1.assign(n, 0);
  m_b2on.assign(n, 0);
  m_b2sj.assign(n, Complex(0., 0.));
}


/*
  Per event: may beta_2 run at all, and which photons can pair?

  Refused (once, with the reason) when the two-photon provider is missing or
  is not this process plus two photons, and for Borns with photons among
  their final legs: the emitted photons' helicity slots among identical
  final photons are the open e+e- -> gamma gamma question, and a wrong slot
  would pass silently.

  A photon can pair when it is ENUMERATED (not collapsed below
  SOFT_PARTITION_CUT or by MAX_PARTITION_PHOTONS: a collapsed photon has no
  stage of its own, only the total eikonal) and above BETA2_XCUT in
  x = 2E/sqrt(s), E in the CEEX frame (the beam CMS), which is KKMC's
  E/E_beam against vcut2.
*/
bool Ceex_Base::PrepareBeta2()
{
  const size_t n(m_allphotons.size());
  m_b2elig.assign(n, 0);
  if (m_order != 2 || n < 2) return false;
  static int state(0);   // 0 unchecked, 1 ok, -1 refused
  if (state == 0) {
    std::string why;
    YFS::Real_Correction *prov(RealProvider(2));
    if (prov == NULL || prov->p_proc == NULL)
      why = "no two-photon real provider (YFS: NLO_Part must contain W, "
            "which builds the process plus two photons)";
    else {
      const Flavour_Vector &pf(prov->p_proc->Flavours());
      const size_t nl(m_flavs.size());
      if (pf.size() != nl + 2) why = "the provider's process has the wrong leg count";
      else {
        for (size_t i(0); i < nl && why.empty(); ++i)
          if (pf[i] != m_flavs[i]) why = "the provider's process has different legs";
        for (size_t i(nl); i < nl + 2 && why.empty(); ++i)
          if (pf[i].Kfcode() != kf_photon) why = "the provider's extra legs are not photons";
      }
    }
    for (size_t i(2); i < m_flavs.size() && why.empty(); ++i)
      if (m_flavs[i].Kfcode() == kf_photon)
        why = "the Born has final-state photons (identical-photon helicity "
              "slots not validated)";
    if (!why.empty()) {
      msg_Error()<<METHOD<<"(): CEEX: ORDER 2 requested but beta_2 is OFF for "
                 <<"this run: "<<why<<". The CEEX column is O(alpha^1)."<<std::endl;
      state = -1;
    }
    else {
      state = 1;
      msg_Info()<<"CEEX: beta_2 from "<<prov->p_proc->Name()<<std::endl;
    }
  }
  if (state < 0) return false;
  static const double xcut(ATOOLS::Settings::GetMainSettings()["CEEX"]
                           ["BETA2_XCUT"].Get<double>());
  const double rs(m_s > 0. ? sqrt(m_s) : 0.);
  if (!(rs > 0.)) return false;
  size_t ne(0);
  for (size_t i(0); i < n; ++i) {
    if (m_fixedstage.size() == n && m_fixedstage[i]) continue;
    if (2.*m_allphotons[i][0]/rs > xcut) { m_b2elig[i] = 1; ++ne; }
  }
  return ne >= 2;
}


bool Ceex_Base::Beta2PairEligible(int j, int l) const
{
  return j >= 0 && l >= 0 && j != l && j < (int)m_b2elig.size()
    && l < (int)m_b2elig.size() && m_b2elig[j] && m_b2elig[l];
}


/*
  Comix's two-photon amplitude, in the alignment the partition sum uses:
  m_cxbalign x (per-photon coupling rescale) x the Comix entry, with the
  fermion flip applied to the index as for the Born and M_1. Photon j is
  the first photon leg, l the second; the amplitude is symmetric under
  swapping them together with their helicities, so which is which only
  fixes the bit. The photon flip bit applies to both (identical particles,
  one convention). RealNorm cancels: M_1 carries it and is divided by it.
*/
bool Ceex_Base::ComixReal2Shifted(const Vec4D_Vector &pp, int hj, int hl,
                                  Amplitude &M2, const PropShifts &shifts,
                                  bool noshift)
{
  YFS::Real_Correction *prov(RealProvider(2));
  if (prov == NULL) return false;
  const int nlg(Amplitude::s_nlegs);
  if (pp.size() != m_flavs.size() + 2) return false;
  const std::vector<METOOLS::Spin_Amplitudes> *amps
    (noshift ? prov->ComixAmplitudes(pp, -1.)
             : prov->ComixAmplitudesShifts(pp, shifts));
  if (amps == NULL || amps->empty()) return false;
  const METOOLS::Spin_Amplitudes &sa((*amps)[0]);
  if ((int)sa.size() != (1 << (nlg + 2))) return false;
  const int fl(m_comixflip), fmaskx((1 << nlg) - 1), pf((fl >> nlg) & 1);
  const int hg(((hj > 0 ? 0 : 1) ^ pf) | (((hl > 0 ? 0 : 1) ^ pf) << 1));
  const Complex cpl(ComixPhotonCoupling(2), 0.);
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f)
    M2.m_A[f] = m_cxbalign.m_A[f] * cpl
      * sa[FlatIndex2((size_t)(f ^ (fl & fmaskx)), hg, nlg)];
  return true;
}


Ceex_Base::PropShifts Ceex_Base::StageShifts2(int j, int l) const
{
  PropShifts sh;
  if (m_stage.size() != m_allphotons.size()) return sh;
  for (size_t g(0); g < m_stagelegs.size(); ++g) {
    const size_t mask(StageShiftMask((int)g));   // as StageShifts
    if (mask == 0) continue;
    sh.push_back(std::make_pair(mask, StagePhotonSum((int)g, j, l)));
  }
  return sh;
}


bool Ceex_Base::PartitionLegs2(int j, int l, Vec4D_Vector &pb) const
{
  const int n((int)m_allphotons.size());
  if (j < 0 || l < 0 || j >= n || l >= n || m_stage.size() != m_allphotons.size())
    return false;
  Vec4D X(m_pceex[0] + m_pceex[1]);
  for (int i(0); i < n; ++i)
    if (i != j && i != l && m_stagereduces[m_stage[i]]) X -= m_allphotons[i];
  Vec4D Y(X - m_allphotons[j] - m_allphotons[l]);
  for (size_t i(2); i < m_pceex.size(); ++i)
    if (i != m_if1 && i != m_if2) Y -= m_pceex[i];
  return LegsAt(X, Y, pb);
}


void Ceex_Base::AddExchangeLineShifts2(int j, int l, PropShifts &sh) const
{
  if (!ExchangeLineShiftsOn() || m_pceex.size() < 4
      || m_flavs.size() != m_pceex.size()) return;
  Vec4D_Vector pb;
  if (!PartitionLegs2(j, l, pb) || pb.size() != m_pceex.size()) return;
  for (size_t i(0); i < m_pceex.size(); ++i) {
    const Vec4D d((i < 2 ? -1. : 1.)*(pb[i] - m_pceex[i]));
    sh.push_back(std::make_pair((((size_t)1) << i)
                                | PHASIC::Process_Base::s_propshiftleg, d));
  }
}


/*
  One pair (j,l) in the current partition, eq. (95) with the other photons'
  soft-factor product W_jl in front (as a product, never as sProd/(s_j s_l),
  which manufactures poles at the zeros of s_j):

    W_jl [ d_j d_l M_2 - s_l beta_1(k_j) - s_j beta_1(k_l) - s_j s_l B_wp ],
    beta_1(k_j) = on_j ? (d_j M_1(k_j) - s'_j B_wp) : 0,

  d = [photon on the reducing stage], s = m_sactu (the partition's soft
  factors), s' the stage eikonal the beta_1 pass subtracted with, B_wp =
  m_partborn0. Must run after the beta_1 loop of the same partition.
*/
bool Ceex_Base::ComixInfraredSubtracted_2_0(int j, int l)
{
  const int n((int)m_allphotons.size());
  if (!m_cxbalignok || j < 0 || l < 0 || j >= n || l >= n
      || m_stage.size() != m_allphotons.size()
      || m_sactu.size() != m_allphotons.size()
      || m_b2on.size() != m_allphotons.size()) return false;
  const bool dj(m_stagereduces[m_stage[j]] != 0), dl(m_stagereduces[m_stage[l]] != 0);
  Complex w(1., 0.);
  for (int i(0); i < n; ++i)
    if (i != j && i != l) w *= m_sactu[i];
  Amplitude M2;
  const bool addm2(dj && dl);
  if (addm2) {
    static const beta1legs::code legsmode(ATOOLS::Settings::GetMainSettings()["CEEX"]
                              ["BETA1_LEGS"].Get<beta1legs::code>());
    const bool rebuild(legsmode == beta1legs::balanced || m_redborn);
    Vec4D_Vector pb;
    PropShifts shifts;
    if (rebuild) { if (!PartitionLegs2(j, l, pb)) { ++m_b2fail; return false; } }
    else {
      pb = m_pceex;
      shifts = StageShifts2(j, l);
      AddExchangeLineShifts2(j, l, shifts);
    }
    Vec4D_Vector pp(pb);
    pp.push_back(m_allphotons[j]);
    pp.push_back(m_allphotons[l]);
    if (!ComixReal2Shifted(pp, m_PhoHel[j], m_PhoHel[l], M2, shifts, rebuild)) {
      ++m_b2fail;
      static bool warned(false);
      if (!warned) {
        warned = true;
        msg_Error()<<METHOD<<"(): Comix two-photon amplitude unavailable; this "
                   <<"pair's beta_2 is dropped. Reported once."<<std::endl;
      }
      return false;
    }
    ++m_b2m2;
  }
  ++m_b2pairs;
  const Complex sj(m_sactu[j]), sl(m_sactu[l]);
  const int nh(Amplitude::NHel());
  double nv(0.), nm2(0.), nss(0.);
  for (int f = 0; f < nh; ++f) {
    const Complex B(m_partborn0.m_A[f]);
    const Complex b1j(m_b2on[j] ? (m_b2hasM1[j] ? m_b2M1[j].m_A[f] : Complex(0., 0.))
                                  - m_b2sj[j]*B : Complex(0., 0.));
    const Complex b1l(m_b2on[l] ? (m_b2hasM1[l] ? m_b2M1[l].m_A[f] : Complex(0., 0.))
                                  - m_b2sj[l]*B : Complex(0., 0.));
    const Complex v(w*((addm2 ? M2.m_A[f] : Complex(0., 0.))
                       - sl*b1j - sj*b1l - sj*sl*B));
    m_AmpBeta2.m_A[f] += v;
    nv += std::norm(v);
    if (addm2) nm2 += std::norm(w*M2.m_A[f]);
    nss += std::norm(w*sj*sl*B);
  }
  static const int tr(ATOOLS::Settings::GetMainSettings()["CEEX"]
                      ["BETA2_TRACE"].Get<int>());
  static int ntr(0);
  if (tr && ntr < tr*64) {
    ++ntr;
    std::string st;
    for (size_t i(0); i < m_stage.size(); ++i) st += ToString(m_stage[i]);
    std::ostringstream o;
    o<<std::setprecision(6)<<"@@@ B2TRACE wp="<<st<<" j="<<j<<" l="<<l
     <<" dj="<<dj<<" dl="<<dl<<" on="<<(int)m_b2on[j]<<(int)m_b2on[l]
     <<" |W s_j s_l B|="<<sqrt(nss)<<" |W M2|="<<sqrt(nm2)
     <<" |beta2|="<<sqrt(nv)<<"\n";
    std::cerr<<o.str();
  }
  return true;
}


/*
  n = 2 closure. With exactly two photons, both able to pair, the O(alpha^2)
  tree amplitude summed over the partitions is M_2 itself, helicity by
  helicity: sum_wp [s s B + s beta_1 + s beta_1 + beta_2] = M_2. Compared
  against Comix's two-photon amplitude on the event's legs with NO shifts
  and no partition machinery at all (plain ComixAmplitudes).
*/
void Ceex_Base::Beta2Closure()
{
  if (m_order != 2 || m_allphotons.size() != 2 || !Beta2PairEligible(0, 1)
      || !m_cxbalignok) return;
  Vec4D_Vector pp(m_pceex);
  pp.push_back(m_allphotons[0]);
  pp.push_back(m_allphotons[1]);
  Amplitude Mex;
  PropShifts none;
  if (!ComixReal2Shifted(pp, m_PhoHel[0], m_PhoHel[1], Mex, none, true)) return;
  double nd(0.), ne(0.), n2(0.), n1(0.), nb(0.);
  const int nh(Amplitude::NHel());
  for (int f = 0; f < nh; ++f) {
    nd += std::norm(m_AmpExpo2.m_A[f] - Mex.m_A[f]);
    ne += std::norm(Mex.m_A[f]);
    n2 += std::norm(m_AmpExpo2.m_A[f]);
    n1 += std::norm(m_AmpExpo1.m_A[f]);
    nb += std::norm(m_AmpBeta2.m_A[f]);
  }
  const double rs(sqrt(m_s));
  const size_t ni(m_isrphotons.size());
  std::ostringstream o;
  o<<std::setprecision(10)<<"@@@ B2CLOS cls="<<(ni == 2 ? "II" : ni == 1 ? "IF" : "FF")
   <<" x0="<<2.*m_allphotons[0][0]/rs<<" x1="<<2.*m_allphotons[1][0]/rs
   <<" h="<<m_PhoHel[0]<<","<<m_PhoHel[1]
   <<" rel="<<(ne > 0. ? sqrt(nd/ne) : -1.)
   <<" rho2="<<n2/4.<<" |M2|^2/4="<<ne/4.<<" rho1="<<n1/4.
   <<" |B2|^2/4="<<nb/4.<<" rho0="<<m_result0<<" nparts="<<m_nparts<<"\n";
  std::cerr<<o.str();
}


/*
  Soft limits of the partition-summed beta_2 at two photons, on BALANCED
  legs: photon j (or l, or both) is scaled by lambda and the radiating pair
  is rebuilt at Y = P - k_j' - k_l' along its own direction (LegsAt), the
  beams unchanged, so every configuration conserves momentum exactly and
  scaling one photon cannot leave an imbalance for the shifts to misread.

  Everything is recomputed on those legs with the same objects the
  partition sum uses: B_wp from ComixBornShifted with StageShifts(-1)
  (m_partborn0 without flux under NO_PSEUDOFLUX 2), M_1 with
  StageShifts(j), M_2 with StageShifts2, stage eikonals on the legs. It
  prints |beta_2|/|M_2| per lambda: the soft theorem makes M_2 -> s M_1 and
  beta_2 finite (next-to-soft), so the ratio must fall linearly.
  s-channel 2 -> 2 without exchange-line shifts only; a diagnostic.
*/
void Ceex_Base::Beta2SoftLimit()
{
  static const int nev(ATOOLS::Settings::GetMainSettings()["CEEX"]
                       ["BETA2_SOFT_TEST"].Get<int>());
  static int done(0);
  if (m_order != 2 || done >= nev || m_allphotons.size() != 2
      || !Beta2PairEligible(0, 1) || !m_cxbalignok) return;
  ++done;
  const std::vector<int> stsave(m_stage);
  const Vec4D_Vector phsave(m_allphotons);
  const Vec4D P(m_pceex[0] + m_pceex[1]);
  const int nh(Amplitude::NHel());
  const int fmaskx(nh - 1);
  const double rn(RealNorm());
  const double rs(sqrt(m_s));
  static const pseudoflux::code pfmode(ATOOLS::Settings::GetMainSettings()["CEEX"]
                          ["NO_PSEUDOFLUX"].Get<pseudoflux::code>());
  const size_t ni(m_isrphotons.size());
  const char *cls(ni == 2 ? "II" : ni == 1 ? "IF" : "FF");
  for (int mode(0); mode < 3; ++mode) {
    std::ostringstream o;
    o<<std::setprecision(6)<<"@@@ B2SOFT cls="<<cls<<" mode="
     <<(mode == 0 ? "k1soft" : mode == 1 ? "k0soft" : "both")
     <<" x0="<<2.*phsave[0][0]/rs<<" x1="<<2.*phsave[1][0]/rs
     <<" pfmode="<<pfmode;
    // each photon's angle to its nearest charged leg, in units of that
    // leg's m/E (the collinear cone), to tell precision from physics
    for (int j(0); j < 2; ++j) {
      double best(1e30), bu(1e30);
      for (size_t i(0); i < m_pceex.size(); ++i) {
        if (m_flavs[i].Charge() == 0.) continue;
        const double c(Vec3D(phsave[j])*Vec3D(m_pceex[i])
                       /(Vec3D(phsave[j]).Abs()*Vec3D(m_pceex[i]).Abs()));
        const double th(acos(Max(-1., Min(1., c))));
        const double u(m_flavs[i].Mass() > 0. ? th/(m_flavs[i].Mass()/m_pceex[i][0]) : th);
        if (th < best) { best = th; bu = u; }
      }
      o<<" th"<<j<<"="<<best<<"("<<bu<<"m/E)";
    }
    for (int e(0); e <= 5; ++e) {
      const double lam(pow(10., -(double)e));
      const double l0(mode == 0 ? 1. : lam), l1(mode == 1 ? 1. : lam);
      const Vec4D k0(l0*phsave[0]), k1(l1*phsave[1]);
      Vec4D Y(P - k0 - k1);
      for (size_t i(2); i < m_pceex.size(); ++i)
        if (i != m_if1 && i != m_if2) Y -= m_pceex[i];
      Vec4D_Vector pl;
      if (!LegsAt(P, Y, pl)) continue;
      m_allphotons[0] = k0; m_allphotons[1] = k1;
      Amplitude B2, M2s, SSB, B1s, SB;
      bool ok(true);
      for (int w0(0); w0 < m_nstages && ok; ++w0)
        for (int w1(0); w1 < m_nstages && ok; ++w1) {
          m_stage[0] = w0; m_stage[1] = w1;
          const bool d0(m_stagereduces[w0] != 0), d1(m_stagereduces[w1] != 0);
          Amplitude C, M10, M11, M2;
          ok = ComixBornShifted(pl, C, StageShifts(-1));
          const Complex s0(StageEikonal(w0, pl, k0, m_PhoHel[0], 0));
          const Complex s1(StageEikonal(w1, pl, k1, m_PhoHel[1], 1));
          if (ok && d0) { Vec4D_Vector pp(pl); pp.push_back(k0);
            ok = ComixRealShifted(pp, m_PhoHel[0], M10, StageShifts(0)); }
          if (ok && d1) { Vec4D_Vector pp(pl); pp.push_back(k1);
            ok = ComixRealShifted(pp, m_PhoHel[1], M11, StageShifts(1)); }
          if (ok && d0 && d1) { Vec4D_Vector pp(pl); pp.push_back(k0); pp.push_back(k1);
            ok = ComixReal2Shifted(pp, m_PhoHel[0], m_PhoHel[1], M2, StageShifts2(0, 1)); }
          if (!ok) break;
          for (int f = 0; f < nh; ++f) {
            const Complex B(m_cxbalign.m_A[f]*C.m_A[f ^ (m_comixflip & fmaskx)]);
            const Complex b10((d0 ? m_cxbalign.m_A[f]*M10.m_A[f]/rn : Complex(0., 0.)) - s0*B);
            const Complex b11((d1 ? m_cxbalign.m_A[f]*M11.m_A[f]/rn : Complex(0., 0.)) - s1*B);
            const Complex m2(d0 && d1 ? M2.m_A[f] : Complex(0., 0.));
            B2.m_A[f]  += m2 - s1*b10 - s0*b11 - s0*s1*B;
            M2s.m_A[f] += m2;
            SSB.m_A[f] += s0*s1*B;
            // the scaled photon's beta_1 with the fixed one's soft factor,
            // and the soft factor times the Born: beta_1/(s_f B) must level
            // off at the (finite) LBK constant as the photon softens
            const Complex bs(mode == 1 ? b10 : b11), sf(mode == 1 ? s1 : s0);
            B1s.m_A[f] += sf*bs;
            SB.m_A[f]  += sf*B;
          }
        }
      if (!ok) { o<<" lam="<<lam<<":fail"; continue; }
      double nb(0.), nm(0.), nss(0.), nb1(0.), nsb(0.);
      for (int f = 0; f < nh; ++f) { nb += std::norm(B2.m_A[f]); nm += std::norm(M2s.m_A[f]);
        nss += std::norm(SSB.m_A[f]); nb1 += std::norm(B1s.m_A[f]); nsb += std::norm(SB.m_A[f]); }
      /*
        lam:A:B:C with A = |beta_2|/|M_2|, B = |beta_2|/|sum s_j s_l B| (the
        Low/LBK criterion: must fall like lambda in all three modes; beta_2
        itself is finite in the single-soft and grows like 1/lambda in the
        double-soft limit), C = |beta_1(k_scaled)|/|s_fixed B| summed over
        partitions (single-soft modes: levels off at the LBK constant).
      */
      o<<" lam="<<lam<<":"<<(nm > 0. ? sqrt(nb/nm) : -1.)
       <<":"<<(nss > 0. ? sqrt(nb/nss) : -1.);
      if (mode < 2) o<<":"<<(nsb > 0. ? sqrt(nb1/nsb) : -1.);
    }
    o<<"\n";
    std::cerr<<o.str();
  }
  m_stage = stsave;
  m_allphotons = phsave;
}


void Ceex_Base::ReportBeta2() const
{
  if (m_order != 2) return;
  msg_Info()<<"CEEX ORDER 2: beta_2 on "<<m_b2events<<" events, "<<m_b2pairs
            <<" pair-partition terms, "<<m_b2m2<<" two-photon amplitudes, "
            <<m_b2fail<<" failures."<<std::endl;
}
