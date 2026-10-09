/*
  The YFS_Handler's diagnostics and probes: per-event print-outs switched on
  by YFS / ACRAIC settings (WEIGHT_PROBE, LO_PROBE, RR_PROBE, ...). None of
  them changes a weight.
*/
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <sstream>
#include <iomanip>
#include "YFS/Main/YFS_Handler.H"
#include "YFS/Main/ISR.H"

using namespace std;
using namespace ATOOLS;
using namespace MODEL;
using namespace YFS;
using namespace PHASIC;
using namespace METOOLS;

namespace {
  /*
    The eikonal of a final-state pair without couplings,
    -(p1/(p1.k) - p2/(p2.k))^2 = 2 p1.p2/(p1.k p2.k) - m1^2/(p1.k)^2 - m2^2/(p2.k)^2.
  */
  double PairEikonal(const Vec4D &k, const Vec4D &p1, const Vec4D &p2)
  {
    return -(p1/(p1*k) - p2/(p2*k)).Abs2();
  }

  /*
    YFS: FSR_DENSITY_GATE (diagnostic, default 0 = off; N > 0 prints the
    first N events with FSR photons). Which eikonal is the density the FSR
    generator actually produced its photons with?

    FSR::GenerateAngles samples the angle from the crude f-bar, which uses the
    PRE-emission velocities; FSR::F then multiplies every photon by the mass
    weight m_f/f-bar, and m_f (KKMC's KarFin dist1) is built from the
    POST-emission pair in its own rest frame. The weight is in the event
    weight (YFS: WEIGHT_MODE Full, the default), so the weighted density is
    m_f: 4 m_f/w'^2 = S(post legs; k), w' the photon energy in the
    post-emission pair frame. One line per FSR photon:
      gpost = 4 m_f/(w'^2 S(post))  - the generator against the eikonal on
              the physical legs (ACRAIC's crude, up to the coupling),
      gpre  = 4 m_f/(w'^2 S(pre))   - against the eikonal on the
              pre-emission legs (YFS.NLO's FSR crude, REAL_SUB_EIK 8),
      pfyy  = (q'+K)^2/q'^2 * FSR::m_yy - ACRAIC's pseudo-flux against the
              generator's rescaling (1 if the flux is the generator's),
      mw    = FSR's mass weight / (m_f/f-bar) (1 if the weight is m_f/f-bar),
      resid = |Q_pre - q1' - q2' - K| (the pre-emission legs, the physical
              legs and the photons in one frame),
    with thk = the photon angle to the nearer physical lepton in units of
    m/E of that lepton and x = 2 k0/sqrt(s).
  */
  void PrintFSRDensityGate(const FSR &fsr, YFS::Dipole &D,
                           const Vec4D_Vector &prefsr, double sqrts,
                           size_t nisr, double m1, double m2)
  {
    const int l(D.Left()), r(D.Right());
    if (l < 0 || r < 0 || l >= (int)prefsr.size() || r >= (int)prefsr.size()) return;
    const Vec4D q1(D.GetNewMomenta(0)), q2(D.GetNewMomenta(1));
    const Vec4D p1(prefsr[l]), p2(prefsr[r]);
    const Vec4D_Vector k(D.GetPhotons());
    Vec4D K;
    for (const Vec4D &g : k) K += g;
    const Vec4D Qp(q1 + q2), res(p1 + p2 - Qp - K);
    const double sq(Qp.Abs2());
    const double pfyy(sq > 0. ? (Qp + K).Abs2()/sq*fsr.m_yy : -1.);
    const double resid(sqrt(sqr(res[0]) + sqr(res[1]) + sqr(res[2]) + sqr(res[3])));
    // the generation entries, by their energy in the post-emission frame
    std::vector<double> wgen;
    for (size_t i(0); i < fsr.m_k0.size(); ++i)
      wgen.push_back(fsr.m_k0[i]*fsr.m_xfact*fsr.m_scalek);
    for (const Vec4D &g : k) {
      const double w((g*Qp)/sqrt(sq));
      size_t best(wgen.size()); double dmin(1e-6);
      for (size_t i(0); i < wgen.size(); ++i)
        if (std::abs(wgen[i]/w - 1.) < dmin) { dmin = std::abs(wgen[i]/w - 1.); best = i; }
      if (best == wgen.size() || best >= fsr.m_dist1.size()) {
        std::cerr<<"@@@ FSRGATE unmatched w="<<w<<std::endl;
        continue;
      }
      const double mf(fsr.m_dist1[best]), fb(fsr.m_dist2[best]);
      const double spost(PairEikonal(g, q1, q2)), spre(PairEikonal(g, p1, p2));
      const double gpost(4.*mf/(w*w*spost)), gpre(4.*mf/(w*w*spre));
      const double mw(best < fsr.m_MassWls.size() && fb > 0. && mf > 0.
                      ? fsr.m_MassWls[best]/(mf/fb) : -1.);
      const double c1(Vec3D(g)*Vec3D(q1)/(Vec3D(g).Abs()*Vec3D(q1).Abs()));
      const double c2(Vec3D(g)*Vec3D(q2)/(Vec3D(g).Abs()*Vec3D(q2).Abs()));
      const bool one(c1 > c2);
      const Vec4D &qn(one ? q1 : q2);
      const double th(acos(Min(1., one ? c1 : c2)));
      const double thk(th/((one ? m1 : m2)/qn[0]));
      std::cerr<<std::setprecision(8)<<"@@@ FSRGATE nisr="<<nisr<<" nfsr="<<k.size()
               <<" x="<<2.*g[0]/sqrts<<" thk="<<thk<<" gpost="<<gpost
               <<" gpre="<<gpre<<" pfyy="<<pfyy<<" mw="<<mw<<" resid="<<resid
               <<" spost="<<spost<<std::endl;
    }
  }
}

/*
  YFS: SIGN_PROBE (diagnostic, default false): one line per event with the
  final weight of the nominal column (ACRAIC when ACRAIC is on) and of the
  squared-level YFS.NLO_EW column on the same event, before the sign clamps,
  with the photon content (visible ISR/FSR multiplicities, hardest energy
  fractions x = 2E/sqrt(s)) and the bare pair mass of legs 2 and 3. For the
  paper's sign table (tab:sign); the Born ME that multiplies both is positive.
*/
void YFS_Handler::ProbeSign(const double w_full)
{
  static const bool sp(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["SIGN_PROBE"].SetDefault(false).Get<bool>());
  if (sp) {
    double xi(0.), xf(0.);
    for (const Vec4D &k : m_ev.m_ISRPhotons) xi = Max(xi, 2.*k[0]/sqrt(m_s));
    for (const Vec4D &k : m_ev.m_FSRPhotons) xf = Max(xf, 2.*k[0]/sqrt(m_s));
    const double mll(m_ev.m_plab.size() > 3 ? (m_ev.m_plab[2] + m_ev.m_plab[3]).Mass() : -1.);
    const double weex(!IsZero(m_ev.m_corr_nominal) ? w_full*m_ev.m_corr_eex/m_ev.m_corr_nominal : 0.);
    std::ostringstream o;   // one write per line: MPI ranks share the stream
    o<<std::setprecision(10)<<"@@@ SIGN nisr="<<m_ev.m_ISRPhotons.size()
     <<" nfsr="<<m_ev.m_FSRPhotons.size()<<" xisr="<<xi<<" xfsr="<<xf
     <<" mll="<<mll<<" wnom="<<w_full<<" weex="<<weex<<"\n";
    std::cerr<<o.str();
  }
}

void YFS_Handler::ProbeLOFinal(const double w_lo, const double w_full,
                               const bool acraic_nom, const double corr_acraic,
                               const double corr_eex)
{
  static const bool lp(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["LO_PROBE"].SetDefault(0).Get<int>() != 0);
  if (lp) {
    std::ostringstream o;
    o<<std::setprecision(10)<<"@@@ LOPROBE2 born="<<m_born<<" wlo="<<w_lo
     <<" wfull="<<w_full<<" corr="<<(acraic_nom ? corr_acraic : corr_eex)
     <<" real="<<m_ev.m_real<<" trial="<<rpa->gen.NumberOfTrials()
     <<" nisr="<<m_ev.m_ISRPhotons.size()<<" nfsr="<<m_ev.m_FSRPhotons.size()<<"\n";
    std::cerr<<o.str();
  }
}

/*
  The FINAL per-column factors, once the IFI real weight and the external
  virtual are in. The ACRAICCMP line in CalculateBeta() prints m_ev.m_real
  before wif is known, which misled a comparison once: with IFI_Real on the
  two columns differ by wif - 1 as well. Part of ACRAIC: WEIGHT_PROBE.
*/
void YFS_Handler::ProbeFinalWeights(const double corr_eex,
                                    const double corr_acraic, const double wif,
                                    const bool acraico2,
                                    const double corr_acraic_o1,
                                    const double corr_acraic_v0,
                                    const double corr_acraic_v1)
{
  static const bool wp(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["WEIGHT_PROBE"].SetDefault(0).Get<int>()!=0);
  if (wp && p_acraic) {
    double ef(0.);
    for (const Vec4D &k : m_ev.m_FSRPhotons) ef = Max(ef, k[0]);
    /*
      Angle of the hardest FSR photon to the nearest final-state lepton,
      in units of the lepton's m/E (the collinear cone), and the two
      eikonals the last photon's beta_1 was divided by (coherent at the
      mapped point, crude on the event).
    */
    double thk(-1.);
    Vec4D kh; for (const Vec4D &k : m_ev.m_FSRPhotons) if (k[0] > kh[0]) kh = k;
    if (kh[0] > 0.)
      for (size_t i(2); i < m_ev.m_plab.size(); ++i) {
        if (!m_flavs[i].IsChargedLepton()) continue;
        const Vec4D &l(m_ev.m_plab[i]);
        const double ct(Vec3D(l)*Vec3D(kh)/(Vec3D(l).Abs()*Vec3D(kh).Abs()));
        const double th(acos(Max(-1., Min(1., ct)))), unit(l.Mass()/l[0]);
        const double u(unit > 0. ? th/unit : th);
        if (thk < 0. || u < thk) thk = u;
      }
    // One write per line: with several MPI ranks on one stream, a line
    // built from many << calls interleaves with other ranks' lines and
    // parses as nonsense (crude weights of 1e71 were seen).
    std::ostringstream o;
    o<<"@@@ WPROBE eex="<<corr_eex<<" acraic="<<corr_acraic
     <<" real="<<m_ev.m_real<<" wif="<<wif
     <<" thk="<<thk
     <<" sloc="<<(p_nlo?p_nlo->m_subloc:0.)
     <<" scr="<<(p_nlo?p_nlo->m_eikeex:0.)
     <<" ff="<<m_ev.m_formfactor<<" yfsw="<<m_ev.m_yfsweight
     <<" nisr="<<m_ev.m_ISRPhotons.size()
     <<" nfsr="<<m_ev.m_FSRPhotons.size()
     <<" xfsr="<<2.*ef/sqrt(m_s)
     <<" mll="<<(p_acraic->SvarQ()>0.?sqrt(p_acraic->SvarQ()):-1.);
    if (acraico2) o<<" acraic1="<<corr_acraic_o1<<" acraicv0="<<corr_acraic_v0
                  <<" acraicv1="<<corr_acraic_v1;
    o<<"\n";
    std::cerr<<o.str();
  }
}

/*
  YFS: LO_PROBE - the factors of the Born-level weight, per event, to compare
  an NLO_Part B run against the YFS.LO stream of an NLO_Part BVR run.
*/
void YFS_Handler::ProbeLO()
{
  static const bool lp(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["LO_PROBE"].SetDefault(0).Get<int>() != 0);
  if (lp) {
    std::ostringstream o;
    o<<std::setprecision(10)<<"@@@ LOPROBE born="<<m_born
     <<" isrw="<<m_isrWeight<<" fsrw="<<m_fsrWeight
     <<" yfsw="<<m_ev.m_yfsweight<<" ff="<<m_ev.m_formfactor
     <<" omv="<<(1.-m_v)<<" nlotype="<<(int)m_nlotype
     <<" nisr="<<m_ev.m_ISRPhotons.size()<<" nfsr="<<m_ev.m_FSRPhotons.size()
     <<" sp="<<m_ev.m_bornMomenta.size()<<"\n";
    std::cerr<<o.str();
  }
}

void YFS_Handler::ProbeLOPre(const double pre_yfsw, const double pre_ff,
                             const double pre_v, const double pre_born)
{
  static const bool lp(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["LO_PROBE"].SetDefault(0).Get<int>() != 0);
  if (lp) {
    std::ostringstream o;
    o<<std::setprecision(12)<<"@@@ LOPRE yfsw="<<m_ev.m_yfsweight/pre_yfsw
     <<" ff="<<m_ev.m_formfactor/pre_ff<<" omv="<<(1.-m_v)/(1.-pre_v)
     <<" born="<<m_born/pre_born<<"\n";
    std::cerr<<o.str();
  }
}

/*
  YFS: RR_PROBE (>= 1): the exactness identities of RR_MODE exact or
  singles_product / RV_MODE remainder, from the SUMS the weight is built
  from. Two-photon events: the real part 1 + (real + rr)/B against the
  two-photon weight the mode targets (NLO_Base::m_rrExactPair: rho_12, or
  rho_12 with the singles' IF subtractions, "@@@ RRID").
  One-photon events: the event's v (virtual/B, what multiplies the real)
  against the v the RV remainder subtracts ("@@@ RVID").
*/
void YFS_Handler::ProbeNNLOEvent()
{
  static const int rrprobe(ATOOLS::Settings::GetMainSettings()["YFS"]["RR_PROBE"]
                           .SetDefault(0).Get<int>());
  if (rrprobe && m_born != 0.) {
    std::ostringstream o;
    o<<std::setprecision(16);
    // per event: the pieces of the squared-level O(alpha^2) weight and the
    // Born-level YFS weight they multiply (yfsweight, form factor and 1 - v
    // are all set before CalculateBeta calls this)
    o<<"@@@ NNLOEV nph="<<p_nlo->m_photons.size()
     <<" nisr="<<m_ev.m_ISRPhotons.size()
     <<" wlo="<<m_ev.m_yfsweight*m_ev.m_formfactor*(1.-m_v)
     <<" born="<<m_born<<" v="<<m_ev.m_nlo_virtual/m_born
     <<" real="<<m_ev.m_nlo_real/m_born<<" rv="<<m_ev.m_nlo_rv/m_born
     <<" rr="<<m_ev.m_nlo_rr/m_born<<" rn="<<m_ev.m_nlo_rn/m_born
     <<" vxr="<<m_ev.m_nlo_vxr/m_born<<"\n";
    if (p_nlo->m_photons.size() == 2 && !IsBad(p_nlo->m_rrExactPair)) {
      const double tot(1. + (m_ev.m_nlo_real + m_ev.m_nlo_rr + m_ev.m_nlo_rn)/m_born);
      o<<"@@@ RRID tot="<<tot<<" exact="<<p_nlo->m_rrExactPair
       <<" rel="<<tot/p_nlo->m_rrExactPair - 1.<<"\n";
    }
    if (p_nlo->m_photons.size() == 1)
      o<<"@@@ RVID vevent="<<m_ev.m_nlo_virtual/m_born<<" vrv="<<p_nlo->VBornEvent()
       <<" rel="<<(m_ev.m_nlo_virtual/m_born - p_nlo->VBornEvent())<<"\n";
    std::cerr<<o.str();
  }
}

/*
  YFS: REAL_STAB diagnostic. Per radiating final-state pair, the residual
  of Q_D(pre-emission) - q_1' - q_2' - K_D, which is what
  NLO_Base::MapMomentaFSRDipole assumes vanishes when it rebuilds the pair
  from the pre-emission legs (MapMomentaFSR starts from the post-emission
  legs plus the photons instead); a non-zero value is the difference
  between the two constructions on a single pair.
*/
void YFS_Handler::ProbeRealStab()
{
  static const bool ds(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>()!=0);
  static long nprint(0);
  if (ds && nprint < 40 && m_ev.m_reallab.size() == m_ev.m_plab.size()) {
    for (auto &D : p_dipoles->GetDipoleFF()) {
      if (!D.IsResonance()) continue;
      Vec4D res(m_ev.m_reallab[D.Left()] + m_ev.m_reallab[D.Right()]
                - m_ev.m_plab[D.Left()] - m_ev.m_plab[D.Right()]);
      for (const Vec4D &k : D.GetPhotons()) res -= k;
      ++nprint;
      std::cerr<<std::setprecision(10)<<"@@@ PREPOST pair=("<<D.Left()<<","<<D.Right()
               <<") nphot="<<D.GetPhotons().size()<<" residual="<<res
               <<" |res|/M="<<(Vec3D(res).Abs()+std::abs(res[0]))/(m_ev.m_reallab[D.Left()]+m_ev.m_reallab[D.Right()]).Mass()
               <<std::endl;
    } }
}

/*
  ACRAIC: REF_PROBE - every momentum of a one-photon event, for an
  independent (offline) evaluation of |M_1|^2/density: the generator's
  Born momenta, the event (m_plab), the s'-frame beams and the photon,
  with the two correction factors, m_born and Comix's Born |M|^2 at
  m_plab and at m_bornMomenta (to identify which point m_born is at).
*/
void YFS_Handler::ProbeReference(const bool haveacraic)
{
  static const bool rp(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["REF_PROBE"].SetDefault(0).Get<int>()!=0);
  if (rp && haveacraic && m_ev.m_ISRPhotons.size() + m_ev.m_FSRPhotons.size() == 1) {
    std::ostringstream o;
    o<<std::setprecision(15);
    auto v4 = [&](const char *lab, const Vec4D &v) {
      o<<" "<<lab<<"="<<v[0]<<","<<v[1]<<","<<v[2]<<","<<v[3]; };
    o<<"@@@ REFP acraic="<<m_ev.m_acraicfactor<<" br="<<m_ev.m_real
     <<" born="<<m_born<<" rho1="<<p_acraic->GetResult()
     <<" rhocr="<<p_acraic->GetRhoCrude()
     <<" me2plab="<<(p_nlo?p_nlo->BornME2At(m_ev.m_plab):-1.)
     <<" me2born="<<(p_nlo?p_nlo->BornME2At(m_ev.m_bornMomenta):-1.);
    for (size_t i(0); i < m_ev.m_bornMomenta.size(); ++i) {
      std::string l("b"+ToString(i)); v4(l.c_str(), m_ev.m_bornMomenta[i]); }
    for (size_t i(0); i < m_ev.m_plab.size(); ++i) {
      std::string l("l"+ToString(i)); v4(l.c_str(), m_ev.m_plab[i]); }
    for (size_t i(0); i < m_ev.m_sprimeBeams.size(); ++i) {
      std::string l("s"+ToString(i)); v4(l.c_str(), m_ev.m_sprimeBeams[i]); }
    for (const Vec4D &k : m_ev.m_ISRPhotons) v4("kI", k);
    for (const Vec4D &k : m_ev.m_FSRPhotons) v4("kF", k);
    o<<"\n";
    std::cerr<<o.str();
  }
}

/*
  The same event's ACRAIC factor next to the fixed-order one it is compared
  with (YFS.ACRAIC / YFS.BR is the event average of their ratio), with the
  hardest ISR and FSR photon energies. Part of ACRAIC: WEIGHT_PROBE.
*/
void YFS_Handler::ProbeAcraicCompare(const bool haveacraic)
{
  static const bool wp(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["WEIGHT_PROBE"].SetDefault(0).Get<int>()!=0);
  if (wp && haveacraic) {
    double ei(0.), ef(0.);
    for (const Vec4D &k : m_ev.m_ISRPhotons) ei = Max(ei, k[0]);
    for (const Vec4D &k : m_ev.m_FSRPhotons) ef = Max(ef, k[0]);
    std::ostringstream o;
    o<<"@@@ ACRAICCMP acraic="<<m_ev.m_acraicfactor<<" br="<<m_ev.m_real
     <<" eex="<<m_ev.m_eex<<" nisr="<<m_ev.m_ISRPhotons.size()
     <<" nfsr="<<m_ev.m_FSRPhotons.size()
     <<" xisr="<<2.*ei/sqrt(m_s)<<" xfsr="<<2.*ef/sqrt(m_s)
     <<" mll="<<(p_acraic->SvarQ()>0.?sqrt(p_acraic->SvarQ()):-1.)
     <<" yfsw="<<m_ev.m_yfsweight
     <<" rho1="<<p_acraic->GetResult()<<" rho0="<<p_acraic->GetResult0()
     <<" rho01="<<p_acraic->GetResult01()
     <<" rhocr="<<p_acraic->GetRhoCrude()<<" born="<<m_born
     <<" v="<<(m_born != 0. ? m_ev.m_nlo_virtual/m_born : 0.)
     <<" dIB="<<p_acraic->ResonanceBornDeltaV()
     <<" G="<<m_bpmc_evtG<<" nalt="<<m_bpmc_evtn
     <<" Gj="<<(p_nlo?p_nlo->m_bpmc_hardG:1.)
     <<" Rj="<<(p_nlo?p_nlo->m_bpmc_hardR:0.)
     <<" subj="<<(p_nlo?p_nlo->m_bpmc_hardsub:0.)
     <<" xh="<<(p_nlo?2.*p_nlo->m_bpmc_hardx/sqrt(m_s):0.)<<"\n";
    std::cerr<<o.str();
  }
}

/*
  What makes a heavy ACRAIC event? The golden stream cannot answer it - it
  carries only rhocrude/result/factor, and adding fields would change the
  md5 the regression depends on. So this is a SEPARATE line: the factor
  alongside the photon multiplicity and the two invariants that drive the
  partition sum, so the tail can be binned against something physical.
*/
void YFS_Handler::ProbeAcraicWeight()
{
  static const bool wp(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["WEIGHT_PROBE"].SetDefault(0).Get<int>()!=0);
  if (wp)
    std::cerr<<"@@@ ACRAICWT factor="<<m_ev.m_acraicfactor
             <<" nphot="<<p_acraic->NPhot()
             <<" mll="<<(p_acraic->SvarQ()>0.?sqrt(p_acraic->SvarQ()):-1.)
             <<" svarY="<<p_acraic->SvarY()
             <<" yfsw="<<m_ev.m_yfsweight
    /*
      rho0/rhocrud is the ONLY ratio here that is bounded: m_AmpExpo0
      accumulates the Born term alone, coherently, and m_rhocrud sums the
      same per-partition objects incoherently, so Cauchy-Schwarz gives
      rho0/rhocrud <= 2^nphot exactly. The full factor r1/r0 is NOT
      bounded - m_AmpExpo1 also carries beta_1 and the virtual, which
      m_rhocrud does not - so a large factor there is not by itself a
      defect. Print both, and test the bound on the one that has it.
    */
             <<" rho0ovcr="<<(p_acraic->GetRhoCrude()>0. ?
                              p_acraic->GetResult0()/p_acraic->GetRhoCrude() : -1.)
             <<" t02="<<p_acraic->T02()<<" t13="<<p_acraic->T13()
             <<" cthk="<<p_acraic->HardestPhotonCosTheta()
             <<std::endl;
}

/*
  Bit-exact reference stream for refactors of the ACRAIC internals
  (env SHERPA_ACRAIC_GOLDEN). The partition sum, the S-factors and the
  amplitude container are all being generalised from the 2 -> 2 /
  ISR-FSR form to arbitrary 2 -> N and arbitrary emission stages; every
  one of those steps has to leave the N = 2 answer IDENTICAL, not close.
  Printed at full precision so a diff is a real test rather than a
  rounded one. Costs nothing when the variable is unset.
*/
void YFS_Handler::ProbeGolden(const double r0, const double r1)
{
  static const bool gold(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["GOLDEN"].Get<int>()!=0);
  if (gold)
    std::cerr<<"@@@ ACRAICGOLD "<<std::setprecision(17)
             <<" rhocrude="<<r0<<" result="<<r1
             <<" factor="<<m_ev.m_acraicfactor<<std::endl;
}

void YFS_Handler::ProbeFSRDensity(YFS::Dipole &D)
{
  static const int gate(ATOOLS::Settings::GetMainSettings()["YFS"]
                        ["FSR_DENSITY_GATE"].SetDefault(0).Get<int>());
  static int ngate(0);
  if (gate > 0 && ngate < gate && !D.GetPhotons().empty()) {
    ++ngate;
    PrintFSRDensityGate(*p_fsr, D, m_ev.m_reallab, sqrt(m_s),
                        m_ev.m_ISRPhotons.size(), m_flavs[D.Left()].Mass(),
                        m_flavs[D.Right()].Mass());
  }
}
