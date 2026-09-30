#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "YFS/NLO/NLO_Base.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include <functional>
#include <map>
#include <cstdlib>
#include <iostream>
#include "YFS/NLO/Virtual.H"
#include "YFS/NLO/VirtualVirtual.H"
#include "YFS/NLO/Photon_Counterterm.H"
#include "MODEL/Main/Model_Base.H"
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include "YFS/NLO/NLO_Base_Internal.H"

using namespace YFS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

double massmin = 2220;
double rcount = 1;
double sumw = 0;

double NLO_Base::PhotonEminNLO() const
{
  static const double fac
    (ATOOLS::Settings::GetMainSettings()["YFS"]["NLO_PHOTON_EMIN"]
     .SetDefault(0.0).Get<double>());
  if (fac<=0.0) return -1.0;
  return fac*0.5*sqrt(m_s)*m_isrcut;
}

/*
  Does the Born have a space-like exchange line (a current with one initial
  leg and part of the final state)? A leg-only propagator shift, which
  COMIX::Amplitude::SetPropShifts applies to such lines and to nothing else,
  either changes |M_0|^2 or leaves it bit for bit. Cached per Born process.
  For YFS: REAL_COMBINE -1.
*/
static bool BornHasExchangeLine(PHASIC::Process_Base *proc, const Vec4D_Vector &p)
{
  static std::map<const PHASIC::Process_Base*, int> cache;
  if (proc == nullptr || p.size() < 4) return false;
  auto it(cache.find(proc));
  if (it != cache.end()) return it->second != 0;
  PHASIC::Process_Base::Prop_Shifts none, leg;
  const double e(1e-3*p[0][0]);
  for (size_t i(0); i < p.size(); ++i)
    leg.push_back(std::make_pair((((size_t)1) << i) | PHASIC::Process_Base::s_propshiftleg,
                                 Vec4D(e, 0.3*e*(i+1), -0.2*e, 0.5*e)));
  std::vector<METOOLS::Spin_Amplitudes> a0, a1;
  double m0(0.), m1(0.);
  if (!proc->BornSpinAmplitudesShifts(p, a0, &m0, none) ||
      !proc->BornSpinAmplitudesShifts(p, a1, &m1, leg)) return false;
  const int has(m0 != m1 ? 1 : 0);
  cache[proc] = has;
  msg_Info()<<"YFS: Born "<<proc->Name()<<(has ? " has" : " has no")
            <<" space-like exchange line; REAL_COMBINE -1 uses the "
            <<(has ? "product" : "sum")<<" over photons."<<std::endl;
  return has != 0;
}

double NLO_Base::CalculateReal() {
  m_wifterms.clear();
  m_subloc = m_eikeex = 0.;   // read by the WPROBE line; stale otherwise
  if (m_coll_real)
    return p_dipoles->CalculateEEX() * m_born;
  if (!m_realtool)
    return 0;
  double real(0);
  m_real_hard1 = 0.;
  m_real_hard2 = 0.;
  m_ifi_prod = 1.;
  m_bpmc_hardx = 0.; m_bpmc_hardG = 1.; m_bpmc_hardR = 0.; m_bpmc_hardsub = 0.;
  static const double trace_thr(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_TRACE"].Get<double>());
  m_realtrace.str(""); m_realtrace.clear();
  double prodw(1.);
  m_rvinfo.clear();
  m_rrphot.clear();
  m_prodfac.assign(m_photons.size(), 1.);
  for (YFS::Photon &g : m_photons) {
    const Vec4D k(g.K());
    // RV_MODE 1: one entry per photon, filled below if the real is evaluated
    if (RVMode() == rvmode::remainder) m_rvinfo.push_back(RVPointInfo());
    m_lastrvinfo = RVPointInfo();
    if (RRMode() == rrmode::exact) m_rrphot.push_back(RRPhotonInfo());   // factor 1 if skipped
    m_lastrrinfo = RRPhotonInfo();
    { static const bool dg(ATOOLS::Settings::GetMainSettings()["YFS"]["PHOTON_DUMP"].Get<int>()!=0);
      if (dg) {
        ATOOLS::Vec4D tot; double eph(0.0);
        for (const YFS::Photon &h : m_photons) { tot+=h.K(); eph+=h.K().E(); }
        double elab(0.0);
        for (size_t j(2);j<m_plab.size();++j) elab+=m_plab[j].E();
        std::cerr<<"@@@ PHOT E="<<k.E()<<" isr="<<(g.IsISR()?1:0)
                 <<" sqrts="<<sqrt(m_s)<<" nph="<<m_photons.size()
                 <<" Ephtot="<<eph<<" Eout="<<elab
                 <<" Etot="<<(eph+elab)<<std::endl;
      } }
    const double phemin(PhotonEminNLO());
    if (phemin>0.0 && k.E()<phemin) { g.m_beta10 = 0.; continue; }
    static constexpr double SOFT_LIMIT_CHECK_TRIGGER_FRAC = 0.2;
    if (m_check_real_sub == CHECK_REAL_SUB_SOFT && (g.IsFSR() || !HasFSR())) {
      if (k.E() < SOFT_LIMIT_CHECK_TRIGGER_FRAC * sqrt(m_s))
        continue;
      CheckRealSub(k, 0);
    }
    static constexpr double COLLINEAR_CHECK_TRIGGER_FRAC = 0.01;
    if (m_check_real_sub == CHECK_REAL_SUB_COLLINEAR && (g.IsFSR() || !HasFSR())) {
      if (k.E() > COLLINEAR_CHECK_TRIGGER_FRAC * sqrt(m_s))
        continue;
      CheckRealCollinearSub(k, 0);
    }
    const double ifi_before(m_ifi_prod);
    double contrib;
    if (g.IsISR() && (m_isr_debug || m_fsr_debug)) {
      contrib = CalculateReal(k);
      double coll = p_dipoles->GetDipoleII().Beta1(k, m_betaorder);
      coll /= p_dipoles->GetDipoleII().Eikonal(k);
      if (contrib != 0)
        m_histograms2d["REAL_COLL_RATIO"]->Insert(k.E(),
                                                  coll * m_born / contrib);
    } else {
      contrib = CalculateReal(k);
    }
    if (RVMode() == rvmode::remainder) m_rvinfo.back() = m_lastrvinfo;
    if (RRMode() == rrmode::exact) m_rrphot.back() = m_lastrrinfo;
    real += contrib;
    // YFS: REAL_COMBINE 1: this photon's factor, its bracket plus the IF
    // interference ratio m_ifi_prod picked up for it (1 without IFI_Real)
    if (m_born != 0.) {
      const double rif(m_ifireal && ifi_before != 0. && !IsBad(m_ifi_prod)
                       ? m_ifi_prod/ifi_before : 1.);
      prodw *= contrib/m_born + rif;
      if (RRMode() == rrmode::exact) m_rrphot.back().factor = contrib/m_born + rif;
      { const size_t ig(&g - &m_photons[0]);
        if (ig < m_prodfac.size()) m_prodfac[ig] = contrib/m_born + rif; }
    }
    if (m_check_real_sub == CHECK_REAL_SUB_SCATTER)
      RecordSubScatter(k, contrib, g.IsISR() ? "realISR" : "realFSR", m_eikeex);
    { static const bool dg2(ATOOLS::Settings::GetMainSettings()["YFS"]["PHOTON_DUMP"].Get<int>()!=0);
      if (dg2) std::cerr<<"@@@ PHC E="<<k.E()<<" isr="<<(g.IsISR()?1:0)
                        <<" contrib="<<contrib
                        <<" failcut="<<(p_real?(p_real->FailCut()?1:0):-1)
                        <<std::endl; }
    g.m_beta10 = contrib;
  }
  HardestBetas(m_photons, [](const YFS::Photon &g) { return g.beta10(); },
               m_real_hard1, m_real_hard2);
  if (m_ifireal && !IsBad(m_ifi_prod)) real += m_born*(m_ifi_prod - 1.);
  /*
    YFS: REAL_COMBINE (name or old integer).
      sum (0): the O(alpha) sum 1 + sum_j delta_j over the event's photons,
        delta_j = beta_1(k_j)/(S~_j B) plus that photon's IF ratio - 1.
      product (1, default): prod_j (1 + delta_j), which agrees with the sum
        at O(alpha) (identical with one photon) and adds the factorised
        beta_2 ~ beta_1 beta_1/beta_0 at O(alpha^2).
      auto (-1): product when the Born has an exchange line, sum otherwise.
  */
  { static const realcombine::code comb(ATOOLS::Settings::GetMainSettings()["YFS"]
        ["REAL_COMBINE"].SetDefault(realcombine::product).Get<realcombine::code>());
    const bool useprod(comb == realcombine::product
                       || (comb == realcombine::automatic
                           && BornHasExchangeLine(p_bornproc, m_plab)));
    // RV_MODE 1's real-virtual is a remainder on top of this product (it
    // supplies neither beta_2 nor v x real), so the product stays on with it
    if (useprod && m_born != 0. && !IsBad(prodw)
        && (!m_rrtool || RRMode() == rrmode::exact) && (!m_realvirt || RVMode() == rvmode::remainder))
      real = m_born*(prodw - 1.); }
  if (trace_thr > 0. && m_born != 0. && std::abs(1. + real/m_born) > trace_thr) {
    Vec4D Q;
    for (size_t i(2); i < m_plab.size(); ++i) Q += m_plab[i];
    std::cerr<<"@@@ RTRACE event nphot="<<m_photons.size()
             <<" BRfactor="<<(1. + real/m_born)
             <<" sqrt_s="<<sqrt(m_s)<<" sqrt_sp="<<Q.Mass()
             <<" born="<<m_born<<" real="<<real<<"\n"
             <<m_realtrace.str()<<"@@@ RTRACE end"<<std::endl;
  }
  return real;
}

// YFS: SUB8_PART (diagnostic): 1 = only the FSR-labelled half of
// REAL_SUB_EIK 8, 2 = only the ISR-labelled half, 0 = both.
static int Sub8Part() {
  static const int p(ATOOLS::Settings::GetMainSettings()["YFS"]["SUB8_PART"]
                     .SetDefault(0).Get<int>());
  return p;
}

static double DipoleFluxShare(const Vec4D &k, const Vec4D &qpost,
                              const Vec4D_Vector &photons)
{
  Vec4D K;
  bool own(false);
  for (const Vec4D &g : photons) { K += g; if (g == k) own = true; }
  const Vec4D Q(qpost + K);
  const double q2(qpost.Abs2()), Q2(Q.Abs2());
  if (!(q2 > 0.) || !(Q2 > 0.)) return 1.;
  if (!own) return (Q + k).Abs2()/Q2;   // an alternative dipole: k added to it
  const double F(Q2/q2);
  double ysum(0.);
  for (const Vec4D &g : photons) ysum += 2.*(g*Q)/Q2;
  const double yk(2.*(k*Q)/Q2);
  return ysum > 0. ? pow(F, yk/ysum) : F;
}

/*
  The FF part of the generator density for photon k: every FF dipole's
  eikonal on its post-emission legs times its flux share (DipoleFluxShare).
  False (and the caller keeps the generation-leg eikonal) under the pole
  scheme, where the radiating pair is the W pair and not the dipole's legs.
*/
static bool GeneratorFSRDensity(const Vec4D &k, YFS::Define_Dipoles &dips,
                                const Vec4D_Vector &postlab, double &sff)
{
  if (dips.PoleActive()) return false;
  double s(0.);
  for (auto &D : dips.GetDipoleFF()) {
    const int l(D.Left()), r(D.Right());
    if (l < 2 || r < 2 || l >= (int)postlab.size() || r >= (int)postlab.size())
      return false;
    s += D.Eikonal(k, postlab[l], postlab[r])
         * DipoleFluxShare(k, postlab[l] + postlab[r], D.GetPhotons());
  }
  if (!(s > 0.) || IsBad(s)) return false;
  sff = s;
  return true;
}

/*
  YFS: SUB8_FSR_LEGS 3 - the same density, but for the emission the
  (n+1)-body point describes rather than the event. REAL_FSR_MAP 2 evaluates
  |M_1(k_j)|^2 with the radiating pair rebuilt at Q_D - k_j (the other photons
  of the pair re-absorbed, the legs along the event's decay axis); the
  single-emission density of that point is S~(point legs; k_j) times the
  point's own flux Q_D^2/(Q_D - k_j)^2. For one photon the point is the event
  and this is GeneratorFSRDensity. With companions the two differ where the
  companions' recoil moves the lepton by a fraction of its dead cone: there
  only the point's eikonal has the collinear structure of the point's
  |M_1(k_j)|^2. kpt is the photon at the point, pt the point's momenta.
*/
static bool PointFSRDensity(const Vec4D &kpt, YFS::Define_Dipoles &dips,
                            const Vec4D_Vector &pt, double &sff)
{
  if (dips.PoleActive()) return false;
  double s(0.);
  for (auto &D : dips.GetDipoleFF()) {
    const int l(D.Left()), r(D.Right());
    if (l < 2 || r < 2 || l >= (int)pt.size() || r >= (int)pt.size()) return false;
    // (q + k)^2/q^2: for the radiating pair q = Q_D - k_j, the point's own
    // flux; for any other pair k added to it (an alternative channel)
    const Vec4D q(pt[l] + pt[r]);
    const double q2(q.Abs2());
    if (!(q2 > 0.)) return false;
    s += D.Eikonal(kpt, pt[l], pt[r]) * (q + kpt).Abs2()/q2;
  }
  if (!(s > 0.) || IsBad(s)) return false;
  sff = s;
  return true;
}

/*
  YFS: ME_PROBE: the photon-lepton angle of g in q, in units of that
  lepton's m/E, minimised over the charged leptons; el receives the energy
  of the lepton that minimises it.
*/
static double DeadConeAngle(const Vec4D_Vector &q, const Vec4D &g,
                            const Flavour_Vector &flavs, double &el)
{
  double best(1e99); el = 0.;
  for (size_t i(2); i < q.size() && i < flavs.size(); ++i) {
    if (!flavs[i].IsChargedLepton()) continue;
    const double ct(Vec3D(q[i])*Vec3D(g)/(Vec3D(q[i]).Abs()*Vec3D(g).Abs()));
    const double th(acos(Max(-1., Min(1., ct)))/(flavs[i].Mass()/q[i][0]));
    if (th < best) { best = th; el = q[i][0]; }
  }
  return best;
}

// S~_II of photon k on the generation legs, the event's pre-emission (Born)
// momenta; 0 without an initial-state dipole
static double GenerationEikonalII(YFS::Define_Dipoles &dips, const Vec4D &k)
{
  if (!dips.HasDipoleII()) return 0.;
  YFS::Dipole &D(dips.GetDipoleII());
  return D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
}

// sum of S~_FF of photon k on the generation legs
static double GenerationEikonalFF(YFS::Define_Dipoles &dips, const Vec4D &k)
{
  double s(0.);
  for (auto &D : dips.GetDipoleFF())
    s += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  return s;
}

// sum of S~_FF of photon k on the legs the dipoles were last built on (for
// the NLO dipoles: the (n+1)-body point's)
static double PointEikonalFF(YFS::Define_Dipoles &dips, const Vec4D &k)
{
  double s(0.);
  for (auto &D : dips.GetDipoleFF())
    s += D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
  return s;
}

// the incoherent crude S~_II + sum S~_FF on the same legs
static double PointCrude(YFS::Define_Dipoles &dips, const Vec4D &k)
{
  double s(0.);
  if (dips.HasDipoleII()) {
    YFS::Dipole &D(dips.GetDipoleII());
    s += D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
  }
  for (auto &D : dips.GetDipoleFF())
    s += D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
  return s;
}

/*
  beta_1(k)/S~ for one photon k of the event,
      tot = (r flux - subloc B/kappa)/density,
  r the real |M_1|^2 at the (n+1)-body point built for k, subloc the
  eikonal it is subtracted with, the density the one the photon was
  generated with (YFS: REAL_SUB_EIK picks both, REAL_FSR_FLUX and Flux_Mode
  the flux, Sub_Mode which eikonal divides). raw returns the numerator only
  (RawBeta). The helpers fill the RealTerms in this order.
*/
double NLO_Base::CalculateReal(Vec4D k, bool raw) {
  RealTerms rt;
  rt.kk = k;
  rt.k = k;
  m_evts += 1;

  msg_Debugging() << METHOD << " raw=" << raw
                  << " k=" << k << " E=" << k.E() << " pt=" << k.PPerp() << "\n";

  MapRealPoint(rt);
  RealME(rt);
  RealFlux(rt);
  RealSubtraction(rt);

  if (!CheckMomentumConservation(rt.p)) {
    msg_Debugging() << METHOD << " momentum conservation failed"
                    << " k.E=" << rt.k.E() << " dip_mass=" << (rt.p[2]+rt.p[3]).Mass() << "\n";
    msg_Error() << "Momentum Conservation fails in " << METHOD << "\n";
    if (m_isr_debug || m_fsr_debug) FillPointHistograms(rt, "");
    return 0;
  }
  if ((rt.p[2] + rt.p[3]).Mass() < massmin)
    massmin = (rt.p[2] + rt.p[3]).Mass();
  if (m_isr_debug || m_fsr_debug) FillPointHistograms(rt, "_pass");
  if (!RealMEUsable(rt)) return 0;

  PointDenominator(rt);
  if (RealSubEik() == realsubeik::multichannel && !PhotonIsFSR(rt.kk))
    MultichannelISRPhoton(rt);
  if (RealSubEik() == realsubeik::multichannel_prefsr && !PhotonIsFSR(rt.kk)
      && p_bornproc != nullptr && Sub8Part() != 1)
    PreFSRMultichannelISRPhoton(rt);
  Sub8Trace(rt);
  AssembleBeta1(rt);

  /*
    YFS: REAL_ALPHA0 (default true since 2026-09-27). The real above is Comix's, with
    the model's alpha on every photon (1/131.9 under G_mu), and the
    subtraction is raised to match it (/m_rescale_alpha), while beta_0 and
    the eikonals carry alpha(0) (USE_MODEL_ALPHA 0). beta_1 is then the hard
    remainder at alpha_model: where the generated channel's S~ B exceeds
    the exact real by far (a Born near its t-channel pole), the one-photon
    factor tends to 1 - alpha_model/alpha(0) = -0.0387, not to 0. CEEX
    rescales its Comix reals to alpha(0) (Ceex_Base::ComixPhotonCoupling).
    true does the same here: beta_1 -> m_rescale_alpha beta_1, the whole real
    correction at alpha(0), as in CEEX. Measured with the same seed, 100k
    events, 0 -> 1: Z-pole mu mu fiducial CEEX/YFS.NLO 1.0056 -> 1.0044,
    m_ff at 60 GeV 1.111 -> 1.079; Bhabha 1.0071 -> 1.0046 and 1.035 ->
    1.007; YFS.NLO rises (its real correction is negative there). false
    restores the alpha_model remainder. No change with USE_MODEL_ALPHA true.
  */
  if (RealAlpha0()) rt.tot *= m_rescale_alpha;

  RecordHardestPhoton(rt);
  RealTrace(rt);
  RealStabilityProbe(rt);
  IFIRealRatio(rt);

  msg_Debugging() << METHOD << " submode=" << m_submode
                  << " r*flux=" << rt.r*rt.flux
                  << " sub=" << rt.subloc * m_born / m_rescale_alpha
                  << " tot=" << rt.tot << "\n";
  if (m_isr_debug)
    m_histograms2d["Real_Flux"]->Insert(
        rt.flux, sqrt(p_dipoles->GetDipoleII().Sprime()));

  if (m_no_subtraction) {
    msg_Debugging() << METHOD << " no_subtraction: returning r/subloc=" << rt.r/rt.subloc << "\n";
    return rt.r / rt.subloc;
  }
  if (IsBad(rt.tot)) ReportBadReal(rt);
  if (m_isr_debug || m_fsr_debug) FillRealHistograms(rt);
  TrackRealAverage(rt.tot);

  if (raw) {
    double rawval = rt.r * rt.flux - rt.subloc * m_born / m_rescale_alpha;
    msg_Debugging() << METHOD << " raw: returning " << rawval << "\n";
    return rawval;
  }

  msg_Debugging() << METHOD << " returning tot=" << rt.tot << "\n";
  StoreRealPointInfo(rt);
  return rt.tot;
}

/*
  The (n+1)-body point for photon k: the NLO dipoles on the event, the
  point from MapMomenta (which moves the photon, rt.k, and sets
  m_map_reduced), the photon appended.
*/
void NLO_Base::MapRealPoint(RealTerms &rt)
{
  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);
  const dipoletype::code fluxtype(p_nlodipoles->WhichResonant(rt.k));

  msg_Debugging() << METHOD << " fluxtype=" << fluxtype << "\n";
  rt.p = m_plab;
  MapMomenta(rt.p, rt.k);

  rt.p.push_back(rt.k);
}

/*
  The real matrix element at the point: masses restored (CheckMasses), the
  NLO dipoles rebuilt on the point's Born legs, r = |M_1|^2 in the
  normalisation of the Born, with the Born-photon symmetry factor and the
  Born-photon multichannel share.
*/
void NLO_Base::RealME(RealTerms &rt)
{
  const double norm = 2. * pow(2 * M_PI, 3);
  Vec4D_Vector &p(rt.p);
  const Vec4D &k(rt.k), &kk(rt.kk);
  /*
    YFS: ME_PROBE - for one final-state photon, the photon-lepton angle in
    units of that lepton's m/E and the lepton energy, in the event (post-
    emission legs, lab photon) and at the point handed to Calc_R, before and
    after CheckMasses; and the real it returns. Dead-cone referee against
    CEEX, which evaluates Comix on the event legs.
  */
  static const bool meprobe(ATOOLS::Settings::GetMainSettings()["YFS"]
                            ["ME_PROBE"].SetDefault(0).Get<int>() != 0);
  double pe_ev(0.), pe_pre(0.), th_ev(-1.), th_pre(-1.);
  const bool probe_this(meprobe && PhotonIsFSR(kk) && m_photons.size() == 1);
  if (probe_this) {
    th_ev  = DeadConeAngle(m_postlab, kk, m_flavs, pe_ev);
    th_pre = DeadConeAngle(p, k, m_flavs, pe_pre);
  }
  const Vec4D_Vector p_before(p);
  CheckMasses(p, 1);

  rt.pp = p;
  rt.pp.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, rt.pp, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, rt.pp, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, rt.pp, m_plab);

  double r = p_real->Calc_R(p) / norm;
  if (probe_this) {
    double pe_post(0.);
    const double th_post(DeadConeAngle(p, p.back(), m_flavs, pe_post));
    double dmax(0.);
    for (size_t i(0); i < p.size(); ++i)
      dmax = Max(dmax, (p[i] - p_before[i]).PSpat() + std::abs((p[i] - p_before[i])[0]));
    std::ostringstream o;
    o<<std::setprecision(10)<<"@@@ MEPROBE x="<<2.*kk.E()/sqrt(m_s)
     <<" th_ev="<<th_ev<<" th_pt="<<th_pre<<" th_stretch="<<th_post
     <<" El_ev="<<pe_ev<<" El_pt="<<pe_pre<<" El_stretch="<<pe_post
     <<" dstretch="<<dmax<<" r/B="<<(m_born != 0. ? r/m_born : 0.)
     <<" kev="<<kk<<" kpt="<<p.back()<<"\n";
    std::cerr<<o.str();
  }
  /*
    Identical photons in the Born final state (YFS: REAL_BORN_PHOTON_SYM).

    The real carries the final-state symmetry factor of the (n+1)-photon
    process, 1/(N_g + 1)!, the Born m_born that of the n-photon one, 1/N_g!,
    so r/m_born is 1/(N_g + 1) of |M_1|^2/|M_0|^2. Define_Dipoles::
    CalculateRealSub divides the coherent eikonal by the same N_g + 1
    (m_N_born_Gamma) so that the SUBTRACTION matches r, but the crude
    eikonal subb that beta_1 is divided by is not divided, so the whole
    beta_1/(S~ B) came out 1/(N_g + 1) of its value: for e+e- -> gamma gamma
    (N_g + 1 = 3) the one-photon Born+real factor was 1 + (w - 1)/3, w the
    exact |M_1|^2 (s'/s)/(e^2 S~ |M_0|^2) at the event, to four digits at
    every photon energy (2026-09-26; w checked against the analytic
    e+e- -> 3 gamma matrix element). Here r and the coherent eikonal taken
    from CalculateRealSub are both multiplied back by N_g + 1, so every
    REAL_SUB_EIK / subtraction mode sees |M_1|^2 and eikonals in one
    normalisation. N_g + 1 = 1 without Born photons: nothing changes there.
    0 restores the old (inconsistent) normalisation.
  */
  r *= BornPhotonSym(1);
  /*
    YFS: REAL_BORN_PHOTON_MULTICHANNEL on: the share of the exact real that
    belongs to the generated assignment of this photon, h_gen/sum_b h_b
    (BornPhotonChannelWeight). Only r is partitioned: each channel keeps its
    own beta_0 and its own subtraction, so the soft bracket (where no other
    channel is open and the factor is exactly 1) is unchanged.
  */
  rt.r_nomc = r;
  r *= BornPhotonChannelWeight(kk);
  rt.r = r;
  m_real = r;
  if (p_real->FailCut()) m_failcut = true;
}

/*
  The flux r is multiplied by: YFS: Flux_Mode, then REAL_FSR_FLUX for a
  final-state photon.
*/
void NLO_Base::RealFlux(RealTerms &rt)
{
  const Vec4D &k(rt.k), &kk(rt.kk);
  double &flux(rt.flux);
  if (m_flux_mode == fluxmode::mapped)
    flux = p_nlodipoles->CalculateFlux(k);
  else if (m_flux_mode == fluxmode::average)
    flux = 0.5 * (p_nlodipoles->CalculateFlux(kk) + p_nlodipoles->CalculateFlux(k));
  else
    flux = p_dipoles->CalculateFlux(kk);

  /*
    Define_Dipoles::CalculateFlux forces fluxtype = initial whenever both ISR
    and FSR are on (the WhichResonant() result in MapRealPoint is never used), so a
    FINAL-state photon received the initial-state flux (Q_X - k)^2/Q_X^2 =
    1 - x, as if it had reduced the beam energy. It has not: the Born scale
    s' is untouched by final-state emission, and the two-body phase space at
    the reduced pair mass differs from the crude one only by the muon
    velocity ratio. e+e- -> mu mu at 0.7 GeV (CMD): the real ME of every FSR
    photon was scaled by 0.70-0.78 before the subtraction, the single-FSR-
    photon events had Born+real/CEEX with median 0.75 and a 90th percentile
    of 6.4 where the two must agree event by event, and the photon spectrum
    grew a 2x bump against KKMC at E_gamma/sqrt(s) = 0.12-0.22.
    Setting flux = 1 (YFS: REAL_FSR_FLUX: no_flux) made the single-FSR-photon
    Born+real/CEEX WORSE (median 0.75 -> 1.67) with the old subtraction.
    With REAL_SUB_EIK multichannel_prefsr (the multichannel crude from the
    pre-FSR system) it is the consistent choice, and the default since
    2026-09-26.

    YFS: REAL_FSR_FLUX: see RealFSRFlux() in NLO_Base_Internal.H.
  */
  const realfsrflux::code fsrflux(RealFSRFlux());
  if (fsrflux != realfsrflux::event_flux && PhotonIsFSR(kk)) {
    if (fsrflux == realfsrflux::no_flux) flux = 1.;
    else if (fsrflux == realfsrflux::flux_squared) flux *= flux;  // (m_ff^2/s')^2
    /*
      own_dipole: the flux of the photon's OWN pair, (Q_D - k)^2/Q_D^2 with Q_D the
      pre-emission momentum of the dipole that radiated k. CalculateFlux
      above takes (Q - k)^2/Q^2 with Q the WHOLE ISR-reduced final state,
      which is the same number for a single resonant pair (Q_D = Q) but
      not with two: the recoil Jacobian the FSR generator produced
      (FSR::RescalePhotons, m_yy) is that of the pair in its own frame,
      and a photon of x = 2E/sqrt(s) = 0.05 at 250 GeV has 1 - x = 0.95
      against a pair-mass ratio of 0.69-0.94 depending on its direction
      relative to the Z boost. See the numbers in the report for
      e+e- -> mu mu tau tau (with REAL_FSR_MAP 2 the fixed-order/CEEX ratio
      on single-FSR-photon events grows with x: 1.00, 1.07, 1.10, 1.10,
      1.28, 1.47 for x in <0.01, 0.01-0.03, 0.03-0.06, 0.06-0.12,
      0.12-0.3, >0.3).
    */
    else if (fsrflux == realfsrflux::own_dipole && m_plab.size() == m_flavs.size()) {
      const YFS::Photon *g(FindPhoton(kk));
      if (g != nullptr && g->Dip() != nullptr) {
        const int l(g->Dip()->Left()), r(g->Dip()->Right());
        if (l >= 2 && r >= 2 && l < (int)m_plab.size() && r < (int)m_plab.size()) {
          const Vec4D Qd(m_plab[l] + m_plab[r]);
          const double q2(Qd.Abs2());
          if (q2 > 0.) flux = (Qd - kk).Abs2()/q2;
        } } }
  }
}

/*
  The subtraction eikonal subloc and the density subb for photon k, and
  the per-photon terms read after CalculateReal() (m_wifterms, m_eikeex,
  m_subloc). Starts from the coherent eikonal at the point and the crude of
  the event; REAL_SUB_EIK then replaces them for its photons, and
  REAL_FSR_FLUX on_crude moves the flux of a final-state photon into subb.
*/
void NLO_Base::RealSubtraction(RealTerms &rt)
{
  // CalculateRealSub is the plain eikonal; the symmetry factor is on r above
  rt.subloc = p_nlodipoles->CalculateRealSub(rt.k);
  rt.subb   = p_dipoles->CalculateRealSubEEX(rt.kk);
  rt.flux0 = rt.flux; rt.subloc0 = rt.subloc; rt.subb0 = rt.subb;   // SUB8_TRACE
  /*
    Which eikonal beta_1 subtracts (YFS: REAL_SUB_EIK, name or old integer;
    default multichannel_prefsr (8) since 2026-09-26, see below). Each name
    is followed by its old integer:
      coherent_point (0): the coherent one at the mapped, post-emission
        point (the old default).
      crude (1): the event's crude, S~_II + S~_FF on the pre-emission legs -
        the density the photon was generated with, so the weight is
        r flux/(S~ B) exactly.
      coherent_born (2): the coherent eikonal on the pre-emission (Born) legs
        of the event, the legs the form factor exponent is built on.
      crude_fsr (3), post_emission_fsr (4), coherent_born_fsr (5),
      assignment_born_fsr (6), multichannel (7), multichannel_prefsr (8):
        described at their branches below and at the SubEik* helpers.
  */
  static const realsubeik::code subeik(RealSubEik());
  if (subeik == realsubeik::crude) rt.subloc = rt.subb;
  else if (subeik == realsubeik::coherent_born) rt.subloc = p_dipoles->CalculateRealSub(rt.kk);
  /*
    crude_fsr (3): the event's crude for FINAL-state photons only. With the
    post-emission map (REAL_FSR_MAP >= 1) the coherent eikonal at the
    (n+1)-body point is 0.2-0.6 of the event's crude for a hard FSR photon,
    so 1 + (r flux - S_pt B)/S~_ev carries an extra (1 - S_pt/S~_ev) that
    the exact one-photon weight r flux/(S~_ev B) does not have. Subtracting
    the crude the photon was generated with removes it; ISR photons keep
    their reduced-point subtraction and denominator (m_map_reduced).
    Measured 2026-09-25, MODE: FSR, Z-pole mu mu, one photon, 3 < theta/(m/E)
    < 30: FO/CEEX 1.00/1.00/1.00/1.00/0.97 for x < 0.01 .. > 0.4 with this
    and REAL_FSR_FLUX 1, against 1.00/1.00/1.01/1.09/1.34 before.
  */
  else if (subeik == realsubeik::crude_fsr && PhotonIsFSR(rt.kk)) rt.subloc = rt.subb;
  else if (subeik == realsubeik::coherent_born_fsr && PhotonIsFSR(rt.kk))
    SubEikCoherentBornFSR(rt);
  else if (subeik == realsubeik::assignment_born_fsr && PhotonIsFSR(rt.kk))
    SubEikAssignmentBornFSR(rt);
  else if (subeik == realsubeik::multichannel && PhotonIsFSR(rt.kk))
    SubEikMultichannelFSR(rt);
  else if (subeik == realsubeik::multichannel_prefsr && PhotonIsFSR(rt.kk)
           && p_bornproc != nullptr && Sub8Part() != 2)
    SubEikPreFSRMultichannelFSR(rt);
  else if (subeik == realsubeik::post_emission_fsr && PhotonIsFSR(rt.kk)
           && m_postlab.size() == m_plab.size())
    SubEikPostEmissionFSR(rt);
  if (RealFSRFlux() == realfsrflux::on_crude && PhotonIsFSR(rt.kk)
      && m_postlab.size() == m_plab.size())
    FluxOnFSRCrude(rt);
  m_wifterms.push_back(rt.subb != 0. ? rt.subloc/rt.subb : 1.);
  m_eikeex = rt.subb;
  m_subloc = rt.subloc;

  msg_Debugging() << METHOD << " r=" << rt.r << " flux=" << rt.flux
                  << " (mode=" << m_flux_mode << ")"
                  << " subloc=" << rt.subloc << " subb=" << rt.subb
                  << " born=" << m_born << " alpha=" << m_rescale_alpha << "\n";
}

/*
  REAL_SUB_EIK post_emission_fsr (4): for FINAL-state photons, subtraction AND
  denominator on the
  POST-emission legs. The FSR generator's crude weight (FSR::Weight,
  m_wt2/m_yy, KKMC's KarFin) already carries the Jacobian from the
  generation variables to the physical momenta, so the density to divide
  by is the eikonal on the physical legs - which is what CEEX's rho_crude
  uses (m_Sprod on m_pceex). Dividing by the eikonal on the PRE-emission
  Born momenta instead gives, for a soft photon whose parent lepton was
  kicked by a hard companion, a bracket (r_j - S~_pre B)/S~_pre of O(1)
  inside the collinear cone (the lepton moved by >> m/E): measured with
  mode 3, Z-pole mu mu at the peak, nfsr >= 3 all soft, 10th percentile
  of the fixed-order factor -5 while CEEX has 0.97. With the
  post-emission legs the point built by MapMomentaFSR (post legs, other
  photons re-absorbed by rescaling, directions kept) has the same
  collinear structure as the crude, and the bracket vanishes as it must.
  The one-photon identity r/(S~_post B) = CEEX is exact by construction.
*/
void NLO_Base::SubEikPostEmissionFSR(RealTerms &rt)
{
  const Vec4D &kk(rt.kk);
  double s(GenerationEikonalII(*p_dipoles, kk));
  for (auto &D : p_dipoles->GetDipoleFF()) {
    const int l(D.Left()), r(D.Right());
    if (l >= 2 && r >= 2 && l < (int)m_postlab.size() && r < (int)m_postlab.size())
      s += D.Eikonal(kk, m_postlab[l], m_postlab[r]); }
  if (s > 0. && !IsBad(s)) { rt.subloc = rt.subb = s; }
}

/*
  REAL_SUB_EIK coherent_born_fsr (5): for FINAL-state photons, the COHERENT
  eikonal on the pre-emission
  legs of the event, |J_II + J_FF|^2 on the Born momenta. Mode 3 (the
  incoherent crude) leaves each photon's soft initial-final interference
  2 Re(J_II J_FF^*) in its bracket, additively and with either sign,
  while m_ifi_prod supplies the same interference multiplicatively:
  traced at the Z pole, events with three soft wide-angle FSR photons
  had brackets -0.6, -0.9, -0.99 each (r/(S~_crude B) = 0.35, 0.06,
  0.007, i.e. destructive interference) and a factor of -3, 10th
  percentile -5. With the coherent eikonal the soft brackets vanish and
  the one-photon weight is still exact: 1 + (r - S_coh B)/S~_cr plus the
  product's (S_coh/S~_cr - 1) is r/(S~_cr B). Mode 0 subtracts the
  coherent eikonal too, but at the MAPPED point, which for a hard photon
  is 0.2-0.6 of the event's - the (1 - S_pt/S~_ev) excess. ISR photons
  keep their reduced-point pair (m_map_reduced).
*/
void NLO_Base::SubEikCoherentBornFSR(RealTerms &rt)
{
  // The same two pieces m_ifi_prod is built from (IFIRealRatio): the incoherent
  // crude on the generation legs plus the initial-final interference on
  // those legs. CalculateRealSub(kk) is NOT this - it sits on the
  // post-emission dipole momenta and reproduced mode 4 exactly.
  const double ifg(p_dipoles->CalculateRealSubIF(rt.kk));
  if (!IsBad(ifg) && !IsZero(rt.subb)) rt.subloc = rt.subb + ifg;
}

/*
  REAL_SUB_EIK assignment_born_fsr (6): mode 3 with the crude carrying the
  Born of each ASSIGNMENT. A hard
  photon labelled FSR at wide angle to its pair has a tiny FF crude, but
  the exact ME does not know the label: at 250 GeV such a photon is really
  the initial-state one that brings an off-peak pair onto the Z peak, so
  r is on-peak-large while m_born (pre-emission pair, off peak) is tiny
  and the weight r/(S~ B) reached 1e4 (clean probe: FO 12919, CEEX 0.30,
  x_fsr 0.39, S~_crude 1e-5). CEEX divides by sum_assignments |s|^2 B at
  that assignment's reduced s', so it gives O(1). Here: subtraction and
  denominator become S~_II B_I/B + S~_FF, B_I the Born at the scaled
  ISR-reduced point of this photon, S~_FF and S~_II on the generation
  legs. At the Z pole B_I/B ~ 1 and this reduces to mode 3.
*/
void NLO_Base::SubEikAssignmentBornFSR(RealTerms &rt)
{
  const Vec4D &kk(rt.kk);
  const double sII(GenerationEikonalII(*p_dipoles, kk));
  const double sFF(GenerationEikonalFF(*p_dipoles, kk));
  double br(1.);
  if (sII + sFF > 0. && ReducedBornISR(kk, br)) {
    const double s6(sII*br + sFF);
    if (s6 > 0. && !IsBad(s6)) { rt.subloc = rt.subb = s6; }
  } else rt.subloc = rt.subb;
}

/*
  REAL_SUB_EIK multichannel (7), FSR-labelled photon: the MULTICHANNEL
  crude. The generator reaches a given final state
  with a photon labelled ISR or FSR, each channel with its own Born and
  Jacobian: g_I = S~_II B_I / flux_I, g_F = S~_FF B_F, the two one-photon
  identities measured separately (ISR: r flux/(S~_II B) exact; FSR:
  r/(S~_FF B_pre) exact). The unbiased weight divides r by g_I + g_F
  WHATEVER the label - mode 6 did so only for FSR-labelled photons and
  counted the ISR region twice (Z pole, 66-76 GeV: 1.20 of KKMC) and
  dropped flux_I. For an FSR-labelled photon B = B_F (the event's Born),
  so subb = S~_FF + S~_II (B_I/B)/flux_I, B_I the Born at (P - k)^2 and
  flux_I = (P - k)^2/P^2; subloc adds the IF interference on the same
  legs (as mode 5) so soft brackets vanish, the product m_ifi_prod
  supplying it multiplicatively. No flux on r (REAL_FSR_FLUX 1). The
  ISR-labelled half is MultichannelISRPhoton.
*/
void NLO_Base::SubEikMultichannelFSR(RealTerms &rt)
{
  const Vec4D &kk(rt.kk);
  double &subloc(rt.subloc), &subb(rt.subb);
  const double sII(GenerationEikonalII(*p_dipoles, kk));
  const double sFF(GenerationEikonalFF(*p_dipoles, kk));
  double bI(1.), s2(1.);
  const double ifg(p_dipoles->CalculateRealSubIF(kk));
  if (sFF > 0. && ReducedBornAlt(kk, -1, bI, s2) && s2 > 0.) {
    const double s7(sFF + sII*bI/s2);
    if (s7 > 0. && !IsBad(s7)) { subb = s7; subloc = s7 + (IsBad(ifg) ? 0. : ifg); }
  } else { subloc = subb + (IsBad(ifg) ? 0. : ifg); }
  rt.flux = 1.;
}

/*
  REAL_SUB_EIK multichannel_prefsr (8), FSR-labelled photon: mode 7 with
  the alternative channel built from the right system.
  Mode 7 took the invariant of the other label from the FULL beams,
  (P -+ k)^2, and normalised to the Born at m_bornMomenta (the Born
  kinematics at the full s), neither of which is the system the event's
  Born m_born sits at once ISR photons are present. Here both come from
  the event's PRE-FSR final state Q_pre (NLO m_plab): an FSR-labelled
  photon's ISR alternative has its Born at (Q_pre - k)^2 with the ISR
  flux (Q_pre - k)^2/Q_pre^2, an ISR-labelled photon's FSR alternative
  at (Q_pre + k)^2; the Born ratios are taken between two points of the
  same construction (PreFSRBornRatio). The subtraction is scaled with the
  denominator, subloc -> subloc*(new subb)/(old subb), so that together
  with m_ifi_prod the one-photon weight is exactly r/(g_I + g_F) and the
  soft bracket is unchanged (the scale -> 1 as k -> 0). The FSR-labelled
  base is mode 5 (subloc = S~ + S_IF on the generation legs), flux 1.
  The ISR-labelled half is PreFSRMultichannelISRPhoton.
*/
void NLO_Base::SubEikPreFSRMultichannelFSR(RealTerms &rt)
{
  const Vec4D &kk(rt.kk);
  double &subloc(rt.subloc), &subb(rt.subb);
  const double sII(GenerationEikonalII(*p_dipoles, kk));
  /*
    YFS: SUB8_FSR_LEGS (2026-09-27; name or old integer; default
    generation (0), the pre-emission eikonal below).
    post_emission (1): the FSR channel density on the POST-emission legs
    times the pair
    flux, S~_FF(q_1', q_2'; k) (q_D + k)^2/q_D^2, q_D the post-emission
    pair of the radiating dipole - what the FSR generator produces
    (FSR::F's mass weight is the eikonal of the post-emission legs) and
    what CEEX's rho_crude uses (m_Sprod on the physical legs times
    m_pflux). The pre-emission eikonal below has the same size outside
    the dead cone but a narrower dead cone (m/E_pre instead of m/E_post):
    Z-pole mu mu, one hard FSR photon at theta < 2 m/E, it is up to 5x
    the generator density, and YFS.NLO fell to 0.2-0.8 of the exact
    |M_1|^2 there while CEEX stayed at 1.00. See
    NOTES-collinear-fsr-referee-2026-09-27.md.
  */
  static const sub8fsrlegs::code fsrlegs(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["SUB8_FSR_LEGS"].SetDefault(sub8fsrlegs::generation).Get<sub8fsrlegs::code>());
  /*
    event_density (2) and point_density (3)
    (NOTES-deadcone-crude-2026-09-28.md): the FSR channel density
    the generator really has, post-emission eikonal times pair flux,
    so that the one-photon weight is r/(g B), the density CEEX divides by
    (hard half of SUB8_FSRSUB 4: g (1 + S_IF/crude), on the new g).
    Unlike 1, the SOFT half keeps the generation-leg scale (gsoft below):
    it is the point's coherent current, the soft limit of r, and must not
    be rescaled by a post/pre eikonal ratio, which leaves soft photons
    next to a lepton kicked by a companion with O(1) brackets
    (post_emission doubled YFS.NLO's error at the Z pole).
    event_density (2): the EVENT's density, GeneratorFSRDensity - S~ on the event's
       post-emission legs times the dipole's flux F_D shared among its
       photons. Right at one photon; with companions it describes a
       different emission from the one |M_1|^2 is evaluated for (the
       REAL_FSR_MAP 2 point re-absorbs the companions, which moves the
       lepton by more than its dead cone for electrons): Z-pole e e,
       YFS.NLO fiducial error x4, CEEX/NLO 1.034 at s'/s 0.85-0.93.
    point_density (3, recommended): the POINT's density,
       PointFSRDensity - S~ on the legs of the (n+1)-body point times its
       own single-emission flux. Identical to event_density at one photon; with companions it is the density
       of the emission the point describes, and it reduces to the old
       generation-leg crude where the companions dominate the recoil.
  */
  const double sFFgen(GenerationEikonalFF(*p_dipoles, kk));
  const bool densityscale(fsrlegs == sub8fsrlegs::event_density
                          || fsrlegs == sub8fsrlegs::point_density);
  const double sFF(Sub8FSRChannelDensity(rt, fsrlegs, sFFgen));
  const double ifg(p_dipoles->CalculateRealSubIF(kk));
  const Vec4D Q(PreFSRSystem()), R(Q - kk);
  const double s2(Q.Abs2() > 0. ? R.Abs2()/Q.Abs2() : 0.);
  double bI(0.);
  if (!(s2 > 0.) || !PreFSRBornRatio(R, bI)) bI = 0.;   // no ISR channel
  const double isrchannel(s2 > 0. ? sII*bI/s2 : 0.);
  const double g(sFF + isrchannel);
  // the scale of the soft (point) subtraction: see SUB8_FSR_LEGS
  // event_density, point_density
  const double gsoft(densityscale ? sFFgen + isrchannel : g);
  /*
    Which subtraction (YFS: SUB8_FSRSUB, name or old integer):
    generation (0): mode 5's, on the generation (pre-FSR) legs, scaled
       with g. The
       one-photon weight is then exactly r/(g_I + g_F) (MODE: FSR cone
       identity 1.00), but r's soft limit sits on the legs of the point
       it is evaluated at, not on the generation legs: with a hard ISR
       photon (250 GeV radiative return) the FSR point's beams are the
       generator's reduced beams, tilted by the ISR transverse momentum,
       and next to a hard FSR companion its leptons are the kicked ones.
       Soft FSR photons there kept O(1) brackets (mean +0.18 per photon
       below x = 1e-3, single photons up to 250), YFS.NLO +6% at 250 GeV
       from this half alone.
    point (1): the coherent current at that point (the default's subloc),
       scaled with g: soft brackets vanish, the identity for hard photons
       does not hold (Z pole +4.4% on Born+real).
    blend (4, default): point for soft photons, generation for hard ones,
       w = 1/(1 + (y/y0)^2), y = 2 k.Q_pre/Q_pre^2, y0 = YFS: SUB8_Y0
       (0.01). 250 GeV FSR half 1.062 -> 0.992 of CEEX, Z pole and the
       cone identity unchanged.
  */
  static const sub8fsrsub::code fsub(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["SUB8_FSRSUB"].SetDefault(sub8fsrsub::blend).Get<sub8fsrsub::code>());
  if (g > 0. && !IsBad(g) && subb > 0.) {
    const double cru(subb);
    const double sgen(g*(1. + (IsBad(ifg) ? 0. : ifg/cru)));
    if (fsub == sub8fsrsub::point) subloc *= gsoft/cru;
    else if (fsub == sub8fsrsub::blend) {
      /*
        blend: the point's coherent subtraction for soft photons, mode 5's
        for hard ones, w = 1/(1 + (y/y0)^2), y = 2 k.Q_pre/Q_pre^2.
      */
      static const double y0(ATOOLS::Settings::GetMainSettings()["YFS"]
                             ["SUB8_Y0"].SetDefault(0.01).Get<double>());
      const double q2(Q.Abs2()), y(q2 > 0. ? 2.*(kk*Q)/q2 : 1.);
      const double w(1./(1. + sqr(y/y0)));
      subloc = w*subloc*gsoft/cru + (1. - w)*sgen;
    }
    else subloc = sgen;
    subb   = g;
    rt.flux = 1.;
  }
  m_sub8[0] = sII; m_sub8[1] = sFF; m_sub8[2] = bI; m_sub8[3] = s2;
}

/*
  The FSR channel density S~_FF of REAL_SUB_EIK multichannel_prefsr on the
  legs SUB8_FSR_LEGS selects; sFFgen, the generation-leg eikonal, where
  the chosen density is not available.
*/
double NLO_Base::Sub8FSRChannelDensity(const RealTerms &rt,
                                       sub8fsrlegs::code fsrlegs, double sFFgen)
{
  const Vec4D &kk(rt.kk);
  double sFF(0.);
  if (fsrlegs == sub8fsrlegs::event_density) {
    if (!GeneratorFSRDensity(kk, *p_dipoles, m_postlab, sFF)) sFF = sFFgen;
  }
  else if (fsrlegs == sub8fsrlegs::point_density) {
    if (!PointFSRDensity(rt.k, *p_dipoles, rt.pp, sFF)) sFF = sFFgen;
  }
  else if (fsrlegs == sub8fsrlegs::post_emission && m_postlab.size() == m_plab.size()) {
    for (auto &D : p_dipoles->GetDipoleFF()) {
      const int l(D.Left()), r(D.Right());
      if (l >= 2 && r >= 2 && l < (int)m_postlab.size() && r < (int)m_postlab.size()) {
        const Vec4D qd(m_postlab[l] + m_postlab[r]);
        const double qd2(qd.Abs2());
        const double F(qd2 > 0. ? (qd + kk).Abs2()/qd2 : 1.);
        sFF += F*D.Eikonal(kk, m_postlab[l], m_postlab[r]);
      }
      else sFF += D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1));
    }
  }
  else sFF = sFFgen;
  return sFF;
}

/*
  REAL_FSR_FLUX: on_crude (3) - the crude a FINAL-state photon is divided by carries
  F = (q + k)^2/q^2 on its final-state part, q the post-emission pair, and
  r carries no flux. This is the crude CEEX divides by (validated against
  KKMC): on single-FSR-photon events K rho_crude(CEEX)/((S~_II + F S~_FF) B)
  has median 0.96 where K rho_crude/(S~ B) has 1.15 and a tail to 2.1.
*/
void NLO_Base::FluxOnFSRCrude(RealTerms &rt)
{
  const Vec4D &kk(rt.kk);
  Vec4D q; for (size_t i = 2; i < m_postlab.size(); ++i) q += m_postlab[i];
  const double q2(q.Abs2());
  if (q2 > 0.) {
    const double F((q + kk).Abs2()/q2);
    const double sII(GenerationEikonalII(*p_dipoles, kk));
    const double sFF(GenerationEikonalFF(*p_dipoles, kk));
    if (sII + F*sFF > 0.) { rt.subb = sII + F*sFF; rt.flux = 1.; m_eikeex = rt.subb; }
  }
}

// ISR_DEBUG / FSR_DEBUG: the photon and the point's first pair, for every
// point (suffix "") or for those that conserve momentum ("_pass")
void NLO_Base::FillPointHistograms(const RealTerms &rt, const std::string &suffix)
{
  m_histograms1d["k_E" + suffix]->Insert(rt.k.E());
  m_histograms1d["k_pt" + suffix]->Insert(rt.k.PPerp());
  m_histograms1d["dip_mass" + suffix]->Insert((rt.p[2] + rt.p[3]).Mass());
}

// false (and CalculateReal returns 0) for a vanishing or non-finite real
bool NLO_Base::RealMEUsable(const RealTerms &rt)
{
  if (IsZero(rt.r)) {
    msg_Debugging() << METHOD << " r=0, returning 0\n";
    m_zero_real_amp++;
    return false;
  }
  if (IsBad(rt.r) || IsBad(rt.flux)) {
    msg_Debugging() << METHOD << " bad point: r=" << rt.r << " flux=" << rt.flux << "\n";
    msg_Error() << "Bad point for YFS Real\n"
                << "  Real ME : " << rt.r << "\n"
                << "  Flux    : " << rt.flux << "\n";
    return false;
  }
  return true;
}

/*
  The CRUDE eikonal of the mapped point: S~_II + S~_FF, each dipole on the
  momenta of the point, incoherently, like CalculateRealSubEEX does for the
  event. Only used when the mapped point is not the event (m_map_reduced).

  beta_1 at the reduced point carries the photon's collinear structure at
  a different overall scale from the event (the 1/x^2 of the scaled photon,
  electron-mass terms), so the bounded object is the residual
  r flux - S~coh B divided by an eikonal OF THE POINT, and the crude sum is
  the one to divide by: it is positive, and |J_II + J_FF|^2 <= 2 (|J_II|^2
  + |J_FF|^2) so wherever residual/S~coh is bounded so is residual/crude.
  Dividing by the coherent S~coh instead (the first version of this
  branch) put the coherent ZEROS of the II-FF interference pattern into the
  denominator: the full |M_1|^2 does not vanish there, its hard remainder
  does not, and the ratio did not stay bounded. Measured on e+e- -> u ubar
  at the Z pole, 100k events: YFS.BR 4904 +- 6.9% against 4481 +- 0.16%
  with the legacy map, the heavy events all having S~coh(point) 20-60x
  below the event's crude (a 5.8 GeV photon: 3.3e-6 against 1.2e-4, x =
  0.93). nu nu has no FF current and was not affected (S~coh = crude).

  Crude(point) ~ crude(event)/x^2 up to the bounded Doppler change of the
  final-state angles, so residual/crude(point) is the event-normalised
  beta_1/S~crude - and the n = 1 formula (r flux - subloc B)/subb, with
  subb the event's crude, is the same object with x = 1.
*/
void NLO_Base::PointDenominator(RealTerms &rt)
{
  rt.subb_loc = rt.subb;
  if (m_map_reduced) {
    rt.subb_loc = PointCrude(*p_nlodipoles, rt.k);
    if (IsZero(rt.subb_loc) || IsBad(rt.subb_loc)) rt.subb_loc = rt.subb;
  }
  // the denominator as it would be without REAL_SUB_EIK's changes to it
  rt.denom0 = m_map_reduced ? rt.subb_loc : rt.subb0;
}

/*
  REAL_SUB_EIK multichannel (7), ISR-labelled photon: the final-state channel term of the
  crude carries its own Born and the ISR flux, f = flux B_F/B with B_F the
  Born at (P + k)^2 (the photon returned to beams and pair). At one photon
  r flux/(B (S~_II + f S~_FF)) = r/(g_I + g_F). The same (f - 1) S~_FF is
  added to the subtraction, which keeps the soft bracket at its old value.
*/
void NLO_Base::MultichannelISRPhoton(RealTerms &rt)
{
  const Vec4D &k(rt.k), &kk(rt.kk);
  double &subb_loc(rt.subb_loc), &subloc(rt.subloc), &subb(rt.subb);
  const double r(rt.r), flux(rt.flux);
  double bF(1.), s2(1.);
  if (ReducedBornAlt(kk, +1, bF, s2)) {
    const double f(flux*bF);
    const double sFFpt(m_map_reduced ? PointEikonalFF(*p_nlodipoles, k)
                                     : GenerationEikonalFF(*p_dipoles, kk));
    static const int t7(ATOOLS::Settings::GetMainSettings()["YFS"]
                        ["SUB7_TRACE"].SetDefault(0).Get<int>());
    static long n7(0);
    const double sb0(subb_loc), sl0(subloc);
    if (!IsBad(f) && !IsBad(sFFpt)) {
      subb_loc += (f - 1.)*sFFpt;
      subloc   += (f - 1.)*sFFpt;
      if (!m_map_reduced) subb += (f - 1.)*sFFpt;
    }
    if (t7 && n7 < t7) { ++n7;
      Vec4D P(m_bornMomenta[0] + m_bornMomenta[1]);
      std::ostringstream o;
      o<<"@@@ SUB7 x="<<2.*kk.E()/sqrt(m_s)<<" flux="<<flux<<" bF="<<bF<<" s2="<<s2
       <<" f="<<f<<" sFFpt="<<sFFpt<<" subb_loc "<<sb0<<"->"<<subb_loc
       <<" subloc "<<sl0<<"->"<<subloc<<" reduced="<<(m_map_reduced?1:0)
       <<" sqrtP2="<<P.Mass()<<" r*flux/B="<<r*flux/m_born
       <<" bracket="<<(r*flux/m_born - subloc)/subb_loc<<"\n";
      std::cerr<<o.str(); }
  }
}

/*
  REAL_SUB_EIK multichannel_prefsr (8), ISR-labelled photon: the FSR channel term of the crude,
  S~_FF with the Born at (Q_pre + k)^2 in the ISR photon's own units
  (times its flux), the whole subtraction scaled with the denominator.
*/
void NLO_Base::PreFSRMultichannelISRPhoton(RealTerms &rt)
{
  const Vec4D &k(rt.k), &kk(rt.kk);
  double &subb_loc(rt.subb_loc), &subloc(rt.subloc), &subb(rt.subb);
  /*
    The FSR channel's invariant is taken at the POINT r is evaluated at,
    (final legs of the point + its photon)^2: the ratio this bracket
    forms is r/(S~ B) there, and at n >= 2 the scaled point
    (MapMomentaScaled) boosts the pair rigidly, which moves both the
    photon-lepton angles (S~_FF of the point) and (Q' + k')^2. With the
    event's (Q_pre + k)^2 but the point's S~_FF the correction removed a
    spurious final-state collinear enhancement of the point from the
    denominator while r kept its share of it (250 GeV, radiative return,
    S~_FF(point) up to 60% of the crude at the point against 4% in the
    event; YFS.NLO +8% in 86-96 GeV). At one photon the point is the
    event and nothing changes.
  */
  const Vec4D Q(PreFSRSystem());
  Vec4D R(Q + kk);
  if (m_map_reduced && rt.pp.size() == m_flavs.size() && Q.Abs2() > 0.) {
    Vec4D Rp(k);
    for (size_t i(2); i < rt.pp.size(); ++i) Rp += rt.pp[i];
    const double r2(Rp.Abs2());
    if (r2 > 0.) R = sqrt(r2/Q.Abs2())*Q;
  }
  double &bF(rt.sub8_bF), &f(rt.sub8_f);
  bF = 0.;
  if (!PreFSRBornRatio(R, bF)) bF = 0.;
  f = rt.flux*bF;
  const double sFFpt(m_map_reduced ? PointEikonalFF(*p_nlodipoles, k)
                                   : GenerationEikonalFF(*p_dipoles, kk));
  // YFS: SUB8_ISR_REDUCED (default true): also rescale when the point is
  // the reduced one (m_map_reduced)
  static const bool isrred(ATOOLS::Settings::GetMainSettings()["YFS"]
                           ["SUB8_ISR_REDUCED"].SetDefault(true).Get<bool>());
  const double nb(subb_loc + (f - 1.)*sFFpt);
  if ((isrred || !m_map_reduced) && !IsBad(nb) && nb > 0. && subb_loc > 0.) {
    const double scale(nb/subb_loc);
    rt.sub8_scale = scale;
    subb_loc = nb;
    subloc  *= scale;
    if (!m_map_reduced) subb *= scale;
  }
}

/*
  YFS: SUB8_TRACE n: "@@@ SUB8" for the first n photons under REAL_SUB_EIK
  multichannel_prefsr, both labels, with the terms before (flux0,
  subloc0, subb0) and after the scheme; clears the m_sub8 scratch the FSR
  half left for it.
*/
void NLO_Base::Sub8Trace(const RealTerms &rt)
{
  const bool se8(RealSubEik() == realsubeik::multichannel_prefsr);
  static const int t8(ATOOLS::Settings::GetMainSettings()["YFS"]
                      ["SUB8_TRACE"].SetDefault(0).Get<int>());
  static long n8(0);
  if (se8 && t8 && n8 < t8) { ++n8;
    const Vec4D &kk(rt.kk);
    const bool fsr(PhotonIsFSR(kk));
    const double r(rt.r), flux(rt.flux), subloc(rt.subloc), subb(rt.subb),
                 subb_loc(rt.subb_loc), flux0(rt.flux0), subloc0(rt.subloc0),
                 subb0(rt.subb0), bF(rt.sub8_bF), f(rt.sub8_f), scale(rt.sub8_scale);
    double me2ref(-1.);
    Vec4D_Vector pr;
    if (PreFSRBornPoint(PreFSRSystem(), pr)) me2ref = BornME2At(pr);
    std::ostringstream o;
    o<<"@@@ SUB8 fsr="<<(fsr?1:0)<<" x="<<2.*kk.E()/sqrt(m_s)
     <<" nph="<<m_photons.size()<<" mpre="<<PreFSRSystem().Mass()
     <<" flux="<<flux<<" r*flux/B="<<r*flux/m_born
     <<" subloc="<<subloc<<" subb="<<(m_map_reduced?subb_loc:subb)
     <<" bracket="<<(r*flux/m_born - subloc)/(m_map_reduced?subb_loc:subb)
     <<" isr:bF="<<bF<<" f="<<f<<" scale="<<scale
     <<" fsr:sII="<<m_sub8[0]<<" sFF="<<m_sub8[1]<<" bI="<<m_sub8[2]<<" s2="<<m_sub8[3]
     <<" Bref/m_born="<<(m_born>0.?me2ref/m_born:-1.)
     <<" bracket0="<<(r*flux0/m_born - subloc0)/(m_map_reduced?subb_loc/scale:subb0)
     <<" flux0="<<flux0<<" subloc0="<<subloc0<<" subb0="<<subb0
     <<" reduced="<<(m_map_reduced?1:0);
    // the photon's dead-cone angle to the post-emission leptons, and the
    // hardest ISR photon and other FSR photon
    double thl(1e9);
    for (size_t i(2); i < m_postlab.size() && i < m_flavs.size(); ++i) {
      if (!m_flavs[i].IsChargedLepton()) continue;
      const Vec4D &lq(m_postlab[i]);
      const double ct(Vec3D(lq)*Vec3D(kk)/(Vec3D(lq).Abs()*Vec3D(kk).Abs()));
      thl = Min(thl, acos(Max(-1., Min(1., ct)))/(lq.Mass()/lq[0]));
    }
    double xisr(0.), xfsr(0.);
    for (const YFS::Photon &gph : m_photons) {
      const double xx(2.*gph.K()[0]/sqrt(m_s));
      if (gph.IsFSR()) { if (gph.K() != kk) xfsr = Max(xfsr, xx); }
      else xisr = Max(xisr, xx); }
    o<<" thl="<<thl<<" xisrmax="<<xisr<<" xfsrother="<<xfsr;
    // which legs the FF eikonal sits on: generation (GetBornMomenta),
    // pre-FSR lab (m_plab), post-FSR lab (m_postlab)
    double sg(0.), sl(0.), sp(0.), dmax(0.);
    for (auto &D : p_dipoles->GetDipoleFF()) {
      const int l(D.Left()), rr(D.Right());
      sg += D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1));
      if (l < (int)m_plab.size() && rr < (int)m_plab.size()) {
        sl += D.Eikonal(kk, m_plab[l], m_plab[rr]);
        dmax = Max(dmax, (D.GetBornMomenta(0) - m_plab[l]).PSpat()
                         + std::abs((D.GetBornMomenta(0) - m_plab[l])[0])); }
      if (l < (int)m_postlab.size() && rr < (int)m_postlab.size())
        sp += D.Eikonal(kk, m_postlab[l], m_postlab[rr]);
    }
    o<<" sFFgen="<<sg<<" sFFpre="<<sl<<" sFFpost="<<sp<<" |dBorn|="<<dmax;
    o<<"\n";
    std::cerr<<o.str();
  }
  m_sub8[0] = m_sub8[1] = m_sub8[2] = m_sub8[3] = -1.;
}

/*
  beta_1/S~ = (r flux - subloc B/kappa)/density, the density by Sub_Mode:
  the subtraction eikonal (local), the crude at the point (global), or the
  event's crude with nothing subtracted (off). The terms are kept for the
  real-virtual and the double real (m_last*), the *0 ones as they were
  before REAL_SUB_EIK.
*/
void NLO_Base::AssembleBeta1(RealTerms &rt)
{
  const double r(rt.r), flux(rt.flux), subloc(rt.subloc), subb(rt.subb),
               subb_loc(rt.subb_loc), subloc0(rt.subloc0), subb0(rt.subb0);
  m_lastflux = flux; m_lastsubloc = subloc;
  m_lastdenom = (m_submode == submode::local ? subloc
                 : m_submode == submode::global ? subb_loc : subb);
  m_lastsubloc0 = subloc0;
  m_lastdenom0 = (m_submode == submode::local ? subloc0
                  : m_submode == submode::global ? rt.denom0 : subb0);
  if (m_submode == submode::local)
    rt.tot = (r * flux - subloc * m_born / m_rescale_alpha) / subloc;
  else if (m_submode == submode::global)
    rt.tot = (r * flux - subloc * m_born / m_rescale_alpha) / subb_loc;
  else if (m_submode == submode::off)
    rt.tot = (r * flux) / subb;
  else
    msg_Error() << METHOD << " unknown YFS subtraction mode " << m_submode << "\n";
}

// WEIGHT_PROBE bookkeeping for the hardest photon of the event: the
// exact-over-crude ratio before the multichannel share, the subtraction
// over the denominator (the "1 - sub" floor), and G.
void NLO_Base::RecordHardestPhoton(const RealTerms &rt)
{
  if (rt.kk.E() >= m_bpmc_hardx) {
    const double den(m_submode == submode::local ? rt.subloc : rt.subb_loc);
    m_bpmc_hardx = rt.kk.E();
    m_bpmc_hardG = m_bpmc_lastG;
    m_bpmc_hardR = (m_born != 0. && den != 0.) ? rt.r_nomc*rt.flux/(m_born*den) : 0.;
    m_bpmc_hardsub = (den != 0.) ? rt.subloc/(m_rescale_alpha*den) : 0.;
  }
}

/*
  YFS: REAL_TRACE: one line per photon into m_realtrace, which
  CalculateReal() prints for events whose Born+real factor exceeds the
  threshold.
*/
void NLO_Base::RealTrace(const RealTerms &rt)
{
  static const double tr(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_TRACE"].Get<double>());
  if (tr > 0.) {
    const Vec4D &k(rt.k), &kk(rt.kk);
    const Vec4D_Vector &p(rt.p);
    const double r(rt.r), flux(rt.flux), subloc(rt.subloc), subb(rt.subb),
                 subb_loc(rt.subb_loc), tot(rt.tot);
    // The lab photon (kk), the mapped photon (k), the reduced beams of the
    // (n+1)-body point, the real ME and the two eikonals it is compared to.
    Vec4D Q;
    for (size_t i(2); i < m_plab.size(); ++i) Q += m_plab[i];
    double cmin(2.);
    for (size_t i(0); i < 2; ++i) {
      const double ct(Vec3D(m_plab[i]) * Vec3D(kk) / (Vec3D(m_plab[i]).Abs() * Vec3D(kk).Abs()));
      cmin = std::min(cmin, 1. - std::fabs(ct));
    }
    const double S(subloc * m_born / m_rescale_alpha);
    const double sII(GenerationEikonalII(*p_dipoles, kk));
    const double sFF(GenerationEikonalFF(*p_dipoles, kk));
    m_realtrace<<std::setprecision(6)
               <<"  gam Elab="<<kk.E()<<" x="<<2.*kk.E()/sqrt(m_s)
               <<" SII="<<sII<<" SFF="<<sFF<<" fsr="<<(PhotonIsFSR(kk)?1:0)
               <<" 1-|cos|="<<cmin
               <<" Emap="<<k.E()
               <<" sqrt_sj="<<(p[0]+p[1]).Mass()
               <<" sqrt_sp="<<Q.Mass()
               <<" M_ff="<<(p[2]+p[3]).Mass()
               <<" r="<<r<<" flux="<<flux<<" rflux="<<r*flux
               <<" S~loc*B="<<S<<" S~loc="<<subloc<<" S~crude="<<subb
               <<" S~crude_pt="<<subb_loc<<" reduced="<<(m_map_reduced?1:0)
               <<" Bpp/B="<<(m_born>0.?BornME2At(rt.pp)/m_born:-1.)
               <<" born="<<m_born
               <<" beta1/S~="<<tot<<" beta1/(S~B)="<<(m_born!=0.?tot/m_born:0.)
               <<" failcut="<<(p_real->FailCut()?1:0)
               <<"\n     klab="<<kk<<" kmap="<<k
               <<"\n     Pa="<<m_bornMomenta[0]<<" Pb="<<m_bornMomenta[1]
               <<" IIborn0="<<p_dipoles->GetDipoleII().GetBornMomenta(0)
               <<"\n     pa_j="<<p[0]<<" pb_j="<<p[1]<<" Qlab="<<Q
               <<"\n";
  }
}

// YFS: REAL_STAB: "@@@ RSTAB", how deep the cancellation in r flux - S B
// goes and how far the mapped photon moved
void NLO_Base::RealStabilityProbe(const RealTerms &rt)
{
  static const bool ds(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>()!=0);
  if (ds) {
    const Vec4D &k(rt.k), &kk(rt.kk);
    const double r(rt.r), flux(rt.flux), subloc(rt.subloc), subb(rt.subb), tot(rt.tot);
    const double S(subloc * m_born / m_rescale_alpha);
    const double num(r * flux - S);
    const double subb_k(p_dipoles->CalculateRealSubEEX(k));
    const double subloc_kk(p_nlodipoles->CalculateRealSub(kk));
    const double dk((k-kk).PSpat()/(kk.PSpat()>0.?kk.PSpat():1.));
    const double coskk(Vec3D(k)*Vec3D(kk)
                       /(Vec3D(k).Abs()*Vec3D(kk).Abs()));
    double cmin(2.), pkmin(1e300);
    for (size_t i(0); i < m_plab.size(); ++i) {
      if (m_flavs[i].Charge() == 0.) continue;
      const Vec4D &q(m_plab[i]);
      const double ct(Vec3D(q) * Vec3D(kk) / (Vec3D(q).Abs() * Vec3D(kk).Abs()));
      cmin = std::min(cmin, 1. - std::fabs(ct));
      pkmin = std::min(pkmin, (q * kk) / (q.E() * kk.E()));
    }
    std::cerr << "@@@ RSTAB x=" << (2. * kk.E() / sqrt(m_s))
              << " 1-|cos|=" << cmin
              << " pk/EE=" << pkmin
              << " Rflux=" << r * flux
              << " S=" << S
              << " num=" << num
              << " depth=" << (num != 0. ? std::fabs(S / num) : -1.)
              << " tot=" << tot << " born=" << m_born
              << " subloc=" << subloc << " subb=" << subb
              << " subb_k=" << subb_k << " subloc_kk=" << subloc_kk
              << " dk=" << dk << " 1-cos_kkk=" << (1.-coskk)
              << std::endl;
  }
}

/*
  YFS: IFI_Real, Sub_Mode global: the photon's initial-final interference
  as the factor 1 + S_IF/S~_crude on the generation legs, multiplied into
  m_ifi_prod for photons above the IF dipoles' omega, with its statistics
  (m_ifi_*).
*/
void NLO_Base::IFIRealRatio(const RealTerms &rt)
{
  const Vec4D &kk(rt.kk);
  const double subloc(rt.subloc), subb(rt.subb);
  const bool ifi_above = (kk.E() > p_dipoles->IFIOmega());
  if (m_ifireal && ifi_above && m_submode == submode::global &&
      !IsZero(subb) && !IsBad(subloc) && !IsBad(subb)) {
    const double cru_gen = p_dipoles->CalculateRealSubEEX(kk);
    const double if_gen  = p_dipoles->CalculateRealSubIF(kk);
    const double ratio = (IsZero(cru_gen) || IsBad(cru_gen) || IsBad(if_gen))
                       ? 1. : 1. + if_gen/cru_gen;
    if (!IsBad(ratio)) {
      m_ifi_prod *= ratio;
      ++m_ifi_n; m_ifi_sum += ratio; m_ifi_sum2 += ratio*ratio;
      const int ib = std::min(4, (int)(10.*kk.E()/sqrt(m_s)));
      if (ib >= 0) {
        ++m_ifi_x_n[ib];
        m_ifi_x_r[ib] += ratio;
        const double ex = subloc/(m_rescale_alpha*subb);
        m_ifi_x_e[ib] += IsBad(ex) ? ratio : ex;
      }
      if (ratio < m_ifi_min) m_ifi_min = ratio;
      if (ratio > m_ifi_max) m_ifi_max = ratio;
    }
    msg_Debugging() << METHOD << " IFI_Real ratio=" << ratio << "\n";
  }
}

void NLO_Base::ReportBadReal(const RealTerms &rt)
{
  msg_Debugging() << METHOD << " tot is NaN/Inf"
                  << " r*flux=" << rt.r*rt.flux
                  << " subloc*born=" << rt.subloc*m_born
                  << " subb=" << rt.subb << "\n";
  msg_Error() << "NLO real is NaN\n"
              << "  R        : " << rt.r << "\n"
              << "  Local  S : " << rt.subloc * m_born << "\n"
              << "  Global S : " << rt.subb << "\n";
}

// ISR_DEBUG / FSR_DEBUG: the IF eikonal, the real and beta_1 over the photon
void NLO_Base::FillRealHistograms(const RealTerms &rt)
{
  const Vec4D &k(rt.k);
  const Vec4D_Vector &p(rt.p);
  m_histograms2d["IFI_EIKONAL"]->Insert(k.Y(), k.PPerp(),
                                        p_nlodipoles->CalculateRealSubIF(k));
  m_histograms2d["REAL_SUB"]->Insert((p[0] + p[1]).Mass(), k.E(), rt.tot / m_born);
  m_histograms2d["REAL"]->Insert(k.E(), k.Theta(), rt.r);
  m_histograms2d["REAL_SUB"]->Insert(k.E(), k.Theta(), rt.tot);
}

// the running average of beta_1/S~ over all photons, for the debug line on
// jumps of more than 10% after the first 1000
void NLO_Base::TrackRealAverage(double tot)
{
  sumw += tot;
  rcount += 1;
  double avg = sumw / rcount;
  if (rcount == 1000)
    m_ravg = avg;
  if (rcount > 1000) {
    double diff = fabs(1. - m_ravg / avg) * 100;
    if (diff > 10) {
      msg_Debugging() << METHOD << " large weight jump: " << diff << "%"
                      << " prev=" << m_ravg << " curr=" << avg
                      << " n=" << rcount << "\n";
      m_ravg = avg;
    }
  }
}

/*
  RV_MODE 1 and RR_MODE 1: the point and the photon's exact-over-crude
  ratio, in the same units and coupling as tot (REAL_ALPHA0's kappa
  included), for the real-virtual remainder and the double real.
*/
void NLO_Base::StoreRealPointInfo(const RealTerms &rt)
{
  const double r(rt.r), flux(rt.flux), tot(rt.tot);
  if (RVMode() == rvmode::remainder && m_born != 0. && m_lastdenom != 0.) {
    m_lastrvinfo.ok = true;
    m_lastrvinfo.p = rt.p;
    m_lastrvinfo.rho = (RealAlpha0() ? m_rescale_alpha : 1.)*r*flux/(m_lastdenom*m_born);
  }
  if (RRMode() == rrmode::exact && m_born != 0. && m_lastdenom != 0. && !IsBad(tot)) {
    m_lastrrinfo.ok = true;
    m_lastrrinfo.delta = tot/m_born;
    m_lastrrinfo.flux = flux;
    m_lastrrinfo.denom = m_lastdenom;
    m_lastrrinfo.subloc = m_lastsubloc/m_lastdenom;
    m_lastrrinfo.rho = (RealAlpha0() ? m_rescale_alpha : 1.)*r*flux/(m_lastdenom*m_born);
    m_lastrrinfo.pt = rt.p;
  }
}

double NLO_Base::CrudeOnLegs(const Vec4D_Vector &pt, const Vec4D &k)
{
  Vec4D_Vector legs(pt.begin(), pt.begin() + m_flavs.size());
  p_nlodipoles->MakeDipolesII(m_flavs, legs, legs);
  p_nlodipoles->MakeDipoles(m_flavs, legs, legs);
  return PointCrude(*p_nlodipoles, k);
}

namespace {
  inline size_t PopCount(unsigned m) {
    size_t n(0);
    for (; m; m &= m-1) ++n;
    return n;
  }
  inline size_t LowestBit(unsigned m) {
    size_t i(0);
    for (; !(m&1u); m >>= 1) ++i;
    return i;
  }
}

bool NLO_Base::RealMEForSubset(unsigned mask, double &me,
                               std::vector<double> &subloc)
{
  const size_t n(PopCount(mask));
  YFS::Real_Correction *prov(RealProvider(n));
  if (prov == NULL) {
    msg_Error()<<METHOD<<"(): no "<<n<<"-photon real ME provider. The subset "
               <<"recursion is general, the matrix element is not."<<std::endl;
    return false;
  }

  Vec4D_Vector ks;
  std::vector<size_t> idx;
  for (size_t i(0); i < m_photons.size(); ++i)
    if (mask & (1u<<i)) { ks.push_back(m_photons[i].K()); idx.push_back(i); }

  Vec4D_Vector p(m_plab);
  MapMomenta(p, ks);                   // boosts ks in place, any multiplicity
  for (const Vec4D &k : ks) p.push_back(k);

  // dipoles of the mapped HARD system (the photons popped back off)
  Vec4D_Vector hard(p);
  for (size_t j(0); j < n; ++j) hard.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, hard, m_plab);
  p_nlodipoles->MakeDipoles  (m_flavs, hard, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, hard, m_plab);

  subloc.assign(m_photons.size(), 0.);
  for (size_t j(0); j < n; ++j)
    subloc[idx[j]] = p_nlodipoles->CalculateRealSub(ks[j]);

  double flux;
  if (m_flux_mode == fluxmode::mapped) {
    Vec4D ksum;
    for (const Vec4D &k : ks) ksum += k;
    flux = p_nlodipoles->CalculateFlux(ksum);
  } else {
    flux = 1.;
    for (const Vec4D &k : ks) flux *= p_dipoles->CalculateFlux(k);
  }

  if (!CheckMomentumConservation(p)) {
    msg_Debugging()<<METHOD<<"(): momentum conservation failed at n="<<n<<"\n";
    return false;
  }

  // (2pi)^{3n} per photon, as the one- and two-photon norms had explicitly
  const double norm(2. * pow(2*M_PI, 3*n));
  const double r(prov->Calc_R(p) / norm * BornPhotonSym(n));
  if (prov->FailCut()) { m_failcut = true; return false; }
  if (IsBad(r) || IsBad(flux)) {
    msg_Debugging()<<METHOD<<"(): bad point r="<<r<<" flux="<<flux<<"\n";
    return false;
  }

  me = r * flux;
  return true;
}

double NLO_Base::RawBeta(unsigned mask)
{
  if (mask == 0) return m_born / m_rescale_alpha;

  std::unordered_map<unsigned, double>::const_iterator it(m_rawbeta.find(mask));
  if (it != m_rawbeta.end()) return it->second;

  double beta(0.);
  if (PopCount(mask) == 1) {
    beta = CalculateReal(m_photons[LowestBit(mask)].K(), /*raw*/true);
  } else {
    double me(0.);
    std::vector<double> subloc;
    if (!RealMEForSubset(mask, me, subloc)) { m_rawbeta[mask] = 0.; return 0.; }
    beta = me;
    // every non-empty subset S of `mask`: eikonal product over S times the
    // residual of what is left. s = (s-1)&mask walks the subsets of mask.
    for (unsigned s(mask); s; s = (s-1) & mask) {
      double eik(1.);
      for (unsigned t(s); t; t &= t-1) eik *= subloc[LowestBit(t)];
      beta -= eik * RawBeta(mask ^ s);
    }
  }
  m_rawbeta[mask] = beta;
  return beta;
}

double NLO_Base::CalculateRealN(unsigned mask)
{
  if (mask == 0) return 0.;
  if (m_photons.size() > 8*sizeof(unsigned)-1) {
    msg_Error()<<METHOD<<"(): "<<m_photons.size()<<" photons exceeds the subset "
               <<"mask width; no fixed-order correction for this event."<<std::endl;
    return 0.;
  }
  const size_t n(PopCount(mask));

  /*
    Cut against photons that the resummation has ALREADY accounted for.

    A photon below the YFS infrared cutoff is unresolved by construction: it
    belongs to the exponentiated form factor, not to a fixed-order residual.
    Computing beta_n for it is a cancellation with no physical content, carried
    out exactly where the numerics are worst - measured on the zpole card, the
    points where the double real disagreed with OpenLoops by factors of 10^7
    are precisely these, with median photon energy 7e-07 GeV against 1.3e-03
    for the points that agree.

    The guard is per photon, so it needs nothing to generalise: a set is
    dropped if ANY of its photons is unresolved, which is what the pair loop
    already did with ||. It is applied HERE as well as in the subset loop so
    that no caller of CalculateRealN can bypass it - the loop only avoids
    doing the work.

    Controlled by NLO_PHOTON_EMIN (multiples of the cutoff energy, 0 = off).
  */
  const double phemin(PhotonEminNLO());
  if (phemin > 0.) {
    for (size_t i(0); i < m_photons.size(); ++i) {
      if ((mask & (1u<<i)) && m_photons[i].K().E() < phemin) {
        const Vec4D &k(m_photons[i].K());
        msg_Debugging()<<METHOD<<"(): photon E="<<k.E()<<" below the resummed "
                       <<"threshold "<<phemin<<", no fixed-order correction\n";
        m_softRR++;
        return 0.;
      }
    }
  }

  // The photons were GENERATED against the crude eikonal, so the weight the
  // caller wants is the residual divided by that crude product - evaluated in
  // the lab, from the dipoles as generated, not from any mapping.
  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles  (m_flavs, m_plab, m_plab);
  double crude(1.);
  for (size_t i(0); i < m_photons.size(); ++i)
    if (mask & (1u<<i)) crude *= p_dipoles->CalculateRealSubEEX(m_photons[i].K());
  m_rr_eik = crude;

  const double beta(RawBeta(mask));

  if (IsZero(crude) || IsBad(crude)) {
    msg_Debugging()<<METHOD<<"(): crude eikonal product "<<crude<<", returning 0\n";
    return 0.;
  }
  const double tot(beta / crude);
  if (IsBad(tot)) {
    msg_Error()<<METHOD<<"(): beta_"<<n<<" is NaN: beta="<<beta
               <<" crude="<<crude<<std::endl;
    return 0.;
  }
  msg_Debugging()<<METHOD<<"(): n="<<n<<" beta="<<beta<<" crude="<<crude
                 <<" tot="<<tot<<" (residuals cached this event: "
                 <<m_rawbeta.size()<<")\n";
  return tot;
}

/*!
  Sum beta_n over every n-photon subset of this event's photons.

  This is the generalisation of the nested i<j loop in CalculateRealReal():
  n=2 reproduces its pair enumeration, n=3 gives triples with no new loop. The
  infrared guard is applied per photon, so it generalises for free - a subset
  is dropped if ANY of its photons is below threshold, which is what the pair
  loop already did with ||.

  Cost is C(m,n) subsets - with 20 photons, 190 pairs but 1140 triples - but
  the residuals themselves are shared: m_rawbeta is keyed on a mask over
  m_photons and lives for the event, so each beta_1 is computed once and each
  beta_2 is reused by every triple built on it. What is still missing is a
  TRUNCATION policy: every subset is enumerated however soft its photons are,
  beyond the resummed-photon cut. HardestBetas is the seed for that.
*/
double NLO_Base::CalculateRealMultiplicity(size_t n)
{
  if (n == 0) return 0.;
  if (RealProvider(n) == NULL) return 0.;
  const size_t m(m_photons.size());
  if (m < n) return 0.;

  if (m > 8*sizeof(unsigned)-1) {
    msg_Error()<<METHOD<<"(): "<<m<<" photons exceeds the subset mask width."
               <<std::endl;
    return 0.;
  }
  // Photons already covered by the resummation get no fixed-order correction;
  // see the note in CalculateRealN. Applied here too so the enumeration can
  // skip the work, not only refuse the result.
  const double phemin(PhotonEminNLO());
  double sum(0.);

  // lexicographic n-subsets of {0..m-1}, as bitmasks over m_photons so that a
  // residual cached by one subset is found by every other subset containing it
  std::vector<size_t> c(n);
  for (size_t i(0); i < n; ++i) c[i] = i;
  while (true) {
    unsigned mask(0);
    bool soft(false);
    for (size_t j(0); j < n; ++j) {
      if (phemin > 0. && m_photons[c[j]].K().E() < phemin) soft = true;
      mask |= (1u << c[j]);
    }
    if (!soft) sum += CalculateRealN(mask);

    size_t i(n);
    while (i > 0 && c[i-1] == m - n + (i-1)) --i;
    if (i == 0) break;
    ++c[i-1];
    for (size_t j(i); j < n; ++j) c[j] = c[j-1] + 1;
  }
  return sum;
}

double NLO_Base::BornME2At(const Vec4D_Vector &p) {
  if (p_bornproc == nullptr) return -1.;
  std::vector<METOOLS::Spin_Amplitudes> amps;
  double me2(0.);
  if (!p_bornproc->BornSpinAmplitudes(p, amps, &me2)) return -1.;
  return (IsBad(me2) || me2 < 0.) ? -1. : me2;
}

/*
  The Born of the photon's INITIAL-state assignment relative to the event's:
  B(s' with k removed as an ISR photon) / B(m_bornMomenta). The point is the
  scaled reduction (KKMC's beta_1 convention) with the label ignored.
*/
bool NLO_Base::ReducedBornISR(const Vec4D &k, double &bratio) {
  bratio = 1.;
  static long nfail[4] = {0,0,0,0};
  auto fail = [&](int i, const char *why) {
    if (nfail[i]++ == 0)
      msg_Error()<<"NLO_Base::ReducedBornISR(): "<<why<<" (reported once)."<<std::endl;
    return false; };
  if (p_bornproc == nullptr) return fail(0, "no Born process pointer");
  if (m_me2born < 0.) m_me2born = BornME2At(m_bornMomenta);
  if (!(m_me2born > 0.)) return fail(1, "Born |M|^2 at the event momenta unavailable");
  /*
    The Born point of the ISR assignment: invariant s'' = (P - k)^2 with P
    the event's (ISR-reduced) beam sum, the same scattering angles. Built in
    the rest frame of R = P - k: beams back to back along the event's beam
    axis, final legs with their event directions rescaled by one xi to the
    total energy sqrt(s''); then boosted onto R. MapMomentaScaled is not
    usable here, its balance check assumes an ISR-labelled photon.
  */
  Vec4D_Vector b(m_bornMomenta);
  if (b.size() < 4 || m_flavs.size() != b.size()) return fail(2, "no Born momenta");
  const Vec4D P(b[0] + b[1]), R(P - k);
  const double s2(R.Abs2());
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (!(s2 > sqr(m1 + m2)) || !(R[0] > 0.)) return fail(2, "reduced s'' below threshold");
  Poincare toP(P);
  Vec4D n4(b[0]); toP.Boost(n4);
  Vec3D n(n4);
  if (!(n.Abs() > 0.)) return fail(2, "degenerate beam axis");
  n = n/n.Abs();
  std::vector<Vec4D> q; std::vector<double> mm2, pp2; double msum(0.);
  for (size_t i(2); i < b.size(); ++i) {
    Vec4D qi(b[i]); toP.Boost(qi); q.push_back(qi);
    const double mi(m_flavs[i].Mass());
    mm2.push_back(mi*mi); pp2.push_back(Vec3D(qi).Sqr()); msum += mi;
  }
  const double M2(sqrt(s2));
  if (msum >= M2) return fail(2, "reduced s'' below the final-state mass");
  auto etot = [&](double xi) { double e(0.);
    for (size_t j(0); j < mm2.size(); ++j) e += sqrt(mm2[j] + xi*xi*pp2[j]);
    return e; };
  double lo(0.), hi(1.);
  while (etot(hi) < M2) hi *= 2.;
  for (int it(0); it < 200; ++it) { const double mid(0.5*(lo+hi)); (etot(mid) < M2 ? lo : hi) = mid; }
  const double xi(0.5*(lo+hi));
  Vec4D_Vector p(b.size());
  const double lam(0.5*sqrt(Lambda(s2, m1*m1, m2*m2)/s2));
  p[0] = Vec4D(sqrt(lam*lam + m1*m1),  lam*n);
  p[1] = Vec4D(sqrt(lam*lam + m2*m2), -lam*n);
  for (size_t i(2); i < b.size(); ++i) {
    const Vec3D v(xi*Vec3D(q[i-2]));
    p[i] = Vec4D(sqrt(mm2[i-2] + v.Sqr()), v);
  }
  Poincare fromR(R);
  for (Vec4D &pi : p) fromR.BoostBack(pi);
  const double me2(BornME2At(p));
  if (!(me2 > 0.)) return fail(3, "Born |M|^2 at the reduced point unavailable");
  bratio = me2/m_me2born;
  return !IsBad(bratio);
}

bool NLO_Base::ReducedBornAlt(const Vec4D &k, int sign, double &bratio,
                              double &s2ratio) {
  s2ratio = 1.;
  bratio = 1.;
  static long nfailalt[4] = {0,0,0,0};
  auto fail = [&](int i, const char *why) {
    if (nfailalt[i]++ == 0)
      msg_Error()<<"NLO_Base::ReducedBornAlt(): "<<why<<" (reported once)."<<std::endl;
    return false; };
  if (p_bornproc == nullptr) return fail(0, "no Born process pointer");
  if (m_me2born < 0.) m_me2born = BornME2At(m_bornMomenta);
  if (!(m_me2born > 0.)) return fail(1, "Born |M|^2 at the event momenta unavailable");
  /*
    The Born point of the ISR assignment: invariant s'' = (P - k)^2 with P
    the event's (ISR-reduced) beam sum, the same scattering angles. Built in
    the rest frame of R = P - k: beams back to back along the event's beam
    axis, final legs with their event directions rescaled by one xi to the
    total energy sqrt(s''); then boosted onto R. MapMomentaScaled is not
    usable here, its balance check assumes an ISR-labelled photon.
  */
  Vec4D_Vector b(m_bornMomenta);
  if (b.size() < 4 || m_flavs.size() != b.size()) return fail(2, "no Born momenta");
  const Vec4D P(b[0] + b[1]), R(sign < 0 ? P - k : P + k);
  const double s2(R.Abs2());
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (!(s2 > sqr(m1 + m2)) || !(R[0] > 0.)) return fail(2, "reduced s'' below threshold");
  Poincare toP(P);
  Vec4D n4(b[0]); toP.Boost(n4);
  Vec3D n(n4);
  if (!(n.Abs() > 0.)) return fail(2, "degenerate beam axis");
  n = n/n.Abs();
  std::vector<Vec4D> q; std::vector<double> mm2, pp2; double msum(0.);
  for (size_t i(2); i < b.size(); ++i) {
    Vec4D qi(b[i]); toP.Boost(qi); q.push_back(qi);
    const double mi(m_flavs[i].Mass());
    mm2.push_back(mi*mi); pp2.push_back(Vec3D(qi).Sqr()); msum += mi;
  }
  const double M2(sqrt(s2));
  if (msum >= M2) return fail(2, "reduced s'' below the final-state mass");
  auto etot = [&](double xi) { double e(0.);
    for (size_t j(0); j < mm2.size(); ++j) e += sqrt(mm2[j] + xi*xi*pp2[j]);
    return e; };
  double lo(0.), hi(1.);
  while (etot(hi) < M2) hi *= 2.;
  for (int it(0); it < 200; ++it) { const double mid(0.5*(lo+hi)); (etot(mid) < M2 ? lo : hi) = mid; }
  const double xi(0.5*(lo+hi));
  Vec4D_Vector p(b.size());
  const double lam(0.5*sqrt(Lambda(s2, m1*m1, m2*m2)/s2));
  p[0] = Vec4D(sqrt(lam*lam + m1*m1),  lam*n);
  p[1] = Vec4D(sqrt(lam*lam + m2*m2), -lam*n);
  for (size_t i(2); i < b.size(); ++i) {
    const Vec3D v(xi*Vec3D(q[i-2]));
    p[i] = Vec4D(sqrt(mm2[i-2] + v.Sqr()), v);
  }
  Poincare fromR(R);
  for (Vec4D &pi : p) fromR.BoostBack(pi);
  const double me2(BornME2At(p));
  if (!(me2 > 0.)) return fail(3, "Born |M|^2 at the reduced point unavailable");
  bratio = me2/m_me2born;
  s2ratio = (P.Abs2() > 0. ? s2/P.Abs2() : 1.);
  return !IsBad(bratio);
}

Vec4D NLO_Base::PreFSRSystem() const {
  Vec4D Q;
  for (size_t i(2); i < m_plab.size(); ++i) Q += m_plab[i];
  return Q;
}

bool NLO_Base::PreFSRBornPoint(const Vec4D &R, Vec4D_Vector &p) const {
  if (m_plab.size() < 4 || m_flavs.size() != m_plab.size()) return false;
  const Vec4D Q(PreFSRSystem());
  const double s2(R.Abs2());
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (!(s2 > sqr(m1 + m2)) || !(R[0] > 0.) || !(Q.Abs2() > 0.)) return false;
  // final legs: their momenta in Q_pre's rest frame, one common rescaling
  Poincare toQ(Q);
  std::vector<Vec3D> q; std::vector<double> mm2; double msum(0.);
  for (size_t i(2); i < m_plab.size(); ++i) {
    Vec4D qi(m_plab[i]); toQ.Boost(qi); q.push_back(Vec3D(qi));
    mm2.push_back(sqr(m_flavs[i].Mass())); msum += m_flavs[i].Mass();
  }
  const double M(sqrt(s2));
  if (msum >= M) return false;
  auto etot = [&](double xi) { double e(0.);
    for (size_t j(0); j < q.size(); ++j) e += sqrt(mm2[j] + xi*xi*q[j].Sqr());
    return e; };
  double lo(0.), hi(1.);
  while (etot(hi) < M) hi *= 2.;
  for (int it(0); it < 200; ++it) { const double mid(0.5*(lo+hi)); (etot(mid) < M ? lo : hi) = mid; }
  const double xi(0.5*(lo+hi));
  p.assign(m_plab.size(), Vec4D());
  // beams along the lab beam axis as a pure boost carries it into R's frame
  const double sgn(m_bornMomenta.size() > 0 && m_bornMomenta[0][3] < 0. ? -1. : 1.);
  const double lam(0.5*sqrt(Lambda(s2, m1*m1, m2*m2)/s2));
  p[0] = Vec4D(sqrt(lam*lam + m1*m1), 0., 0.,  sgn*lam);
  p[1] = Vec4D(sqrt(lam*lam + m2*m2), 0., 0., -sgn*lam);
  for (size_t i(2); i < m_plab.size(); ++i) {
    const Vec3D v(xi*q[i-2]);
    p[i] = Vec4D(sqrt(mm2[i-2] + v.Sqr()), v);
  }
  Poincare fromR(R);
  for (Vec4D &pi : p) fromR.BoostBack(pi);
  return true;
}

bool NLO_Base::PreFSRBornRatio(const Vec4D &R, double &bratio) {
  bratio = 1.;
  static long nfail[3] = {0,0,0};
  auto fail = [&](int i, const char *why) {
    if (nfail[i]++ == 0)
      msg_Error()<<"NLO_Base::PreFSRBornRatio(): "<<why<<" (reported once)."<<std::endl;
    return false; };
  if (p_bornproc == nullptr) return fail(0, "no Born process pointer");
  Vec4D_Vector p;
  if (m_me2pre < 0.) {
    if (!PreFSRBornPoint(PreFSRSystem(), p)) return fail(1, "no pre-FSR Born point");
    m_me2pre = BornME2At(p);
  }
  if (!(m_me2pre > 0.)) return fail(1, "Born |M|^2 at the pre-FSR point unavailable");
  if (!PreFSRBornPoint(R, p)) return false;      // below threshold: no channel
  const double me2(BornME2At(p));
  if (!(me2 >= 0.)) return fail(2, "Born |M|^2 at the alternative point unavailable");
  bratio = me2/m_me2pre;
  return !IsBad(bratio);
}

double NLO_Base::BornPhotonSym(size_t nextra) const {
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["REAL_BORN_PHOTON_SYM"].SetDefault(true).Get<bool>());
  if (!on) return 1.;
  size_t nb(0);
  for (const Flavour &f : m_flavs) if (f.IsPhoton()) ++nb;
  /*
    C(nb + nextra, nextra), not (nb + nextra)!/nb!. YFS sums its photons as
    an unordered set (pairs i < j), so only the mixing of the extra photons
    with the Born's photons needs a factor. The one-photon value is nb + 1
    either way (gamma gamma: 3, validated); the double real of an s-channel
    Born was 2, and the soft limit of beta_2/(S~_1 S~_2) then tended to +B
    per photon pair instead of 0 (Z-pole mu mu, median 1.65 for x < 1e-4;
    YFS.NLO+RR 3192 pb against YFS.NLO 1185 pb). Measured 2026-09-26.
  */
  double sym(1.);
  for (size_t j(1); j <= nextra; ++j) sym *= double(nb + j)/double(j);
  return sym;
}

/*
  Born-photon multichannel (YFS: REAL_BORN_PHOTON_MULTICHANNEL, name or old
  integer: off (0); on (1, default), the weight below applied; report (2),
  G computed and reported for WEIGHT_PROBE but not applied).

  e+e- -> gamma gamma (and any Born with final-state photons): the generator
  picks the Born pair, then adds ISR photons. Once an ISR photon passes the
  generation cuts itself, the same gamma gamma gamma final state is also
  reached with THAT photon in the Born pair and one of the Born photons as
  the ISR photon, each assignment with its own crude density
      h = prod_{ISR} S~_II(k) * B(pair at its reduced point) * s/s'
  (the density the one-photon weight r flux/(S~ B) divides by, flux = s'/s).
  Each channel on its own already integrates the exact real to
  |M_{n+1}|^2 (BornPhotonSym), so with several channels open the real is
  counted once per channel. The standard multichannel weight divides the
  exact real by the SUM of the densities of all open channels, i.e. the
  generated channel's real is multiplied by h_gen/sum_b h_b = 1/G. Any
  partition of unity is unbiased; this one is the generator's density up to
  the Born integrator's own importance sampling, which cancels channel by
  channel. G = 1 exactly when no other assignment passes the cuts (soft
  photons can never be Born photons), so the soft limit, and every process
  without Born photons, is untouched.

  Measured on the Z-pole gamma gamma card (BR, Comix real, 1.6M weighted
  events per run), 2026-09-27, NOTES-gammagamma-multichannel-2026-09-27.md:
  the diphoton mass m_gg in 63.8-71.1 GeV - Born-eligible third photon with
  the generation cut E > 0.20 sqrt(s), not with 0.25 - had EFRAC 0.20/0.25
  = 1.078 +- 0.008 (YFS.NLO), 1.170 +- 0.010 (YFS.CEEX); with this switch
  CEEX 1.010 +- 0.010 and sigma_fid 26.561 in both runs. YFS.NLO then
  shows 0.910 +- 0.007: the beta_0 of every open channel stays, and its
  alpha-scheme floor (YFS: REAL_ALPHA0) scales with the channel count; with
  REAL_ALPHA0 1 as well, 0.986 +- 0.008. One-photon events: the CEEX
  factor times G equals YFS.NLO's exact/crude times alpha(0)/alpha_model
  (0.962747) to 1e-6 in 99.8% of events, G > 2 included.
*/
namespace {
  // The generator's reduced Born point of a final state: beams back to back
  // along the lab axis in the final state's rest frame (pure boost), the
  // final legs as they are.
  bool GeneratorPointOf(const ATOOLS::Flavour_Vector &fl, const Vec4D_Vector &fin,
                        double zsign, Vec4D_Vector &p) {
    Vec4D Q;
    for (const Vec4D &q : fin) Q += q;
    const double s2(Q.Abs2());
    const double m1(fl[0].Mass()), m2(fl[1].Mass());
    if (!(s2 > sqr(m1 + m2)) || !(Q[0] > 0.)) return false;
    const double lam(0.5*sqrt(Max(0., sqr(s2 - m1*m1 - m2*m2) - 4.*m1*m1*m2*m2)/s2));
    Vec4D b0(sqrt(lam*lam + m1*m1), 0., 0.,  zsign*lam);
    Vec4D b1(sqrt(lam*lam + m2*m2), 0., 0., -zsign*lam);
    Poincare fromQ(Q);
    fromQ.BoostBack(b0); fromQ.BoostBack(b1);
    p.clear(); p.push_back(b0); p.push_back(b1);
    for (const Vec4D &q : fin) p.push_back(q);
    return true;
  }
}

double NLO_Base::BornPhotonChannelSum(const Vec4D_Vector &lab,
                                      const std::vector<Vec4D> &isr,
                                      YFS::Dipole &dII, int jswap, int *nalt) {
  if (nalt) *nalt = 0;
  if (p_bornproc == nullptr || lab.size() != m_flavs.size() || isr.empty())
    return 1.;
  std::vector<size_t> slot;                     // Born photon legs
  for (size_t i(2); i < m_flavs.size(); ++i)
    if (m_flavs[i].IsPhoton()) slot.push_back(i);
  if (slot.empty()) return 1.;
  const double zs(m_bornMomenta.size() > 0 && m_bornMomenta[0][3] < 0. ? -1. : 1.);
  const double s(m_bpmc_s > 0. ? m_bpmc_s
                 : (m_bornMomenta.size() > 1 ? (m_bornMomenta[0] + m_bornMomenta[1]).Abs2()
                                             : m_s));
  const Vec4D b0(dII.GetBornMomenta(0)), b1(dII.GetBornMomenta(1));
  auto eik = [&](const Vec4D &k) { return dII.Eikonal(k, b0, b1); };
  PHASIC::Combined_Selector *sel(p_bornproc->Selector());
  const int res0(sel ? sel->Result() : 1);
  Vec4D_Vector fin0(lab.begin() + 2, lab.end());
  /*
    h of the assignment "the photons 'pick' are the Born photons (in the order
    of slot), everything else is ISR", relative factor only: the S~ of the
    photons that leave the Born divided by those that enter it, times the
    Born and the flux. Returns -1 when the assignment is not a channel.
  */
  // candidates: Born photons first, then the ISR photons
  std::vector<Vec4D> cand;
  for (size_t i : slot) cand.push_back(lab[i]);
  const size_t nb(slot.size());
  for (const Vec4D &k : isr) cand.push_back(k);
  auto hof = [&](const std::vector<size_t> &pick, bool trig) -> double {
    Vec4D_Vector fin(fin0);
    for (size_t a(0); a < nb; ++a) fin[slot[a] - 2] = cand[pick[a]];
    Vec4D_Vector p;
    if (!GeneratorPointOf(m_flavs, fin, zs, p)) return -1.;
    const double sp((p[0] + p[1]).Abs2());
    if (trig) {
      if (s > 0. && 1. - sp/s > m_bpmc_vmax) return -1.;
      if (sel && !sel->Trigger(p)) return -1.;
    }
    const double me2(BornME2At(p));
    if (!(me2 > 0.)) return -1.;
    double h(me2*s/sp);
    // S~ of every photon of this assignment that is ISR; the generated
    // assignment's ISR photons are the common reference, so only photons
    // that changed role matter: Born photons that became ISR multiply,
    // ISR photons that became Born divide.
    std::vector<bool> inborn(cand.size(), false);
    for (size_t a(0); a < nb; ++a) inborn[pick[a]] = true;
    for (size_t c(0); c < cand.size(); ++c) {
      const bool born0(c < nb);
      if (born0 && !inborn[c]) h *= eik(cand[c]);
      if (!born0 && inborn[c]) h /= eik(cand[c]);
    }
    return (IsBad(h) || h < 0.) ? -1. : h;
  };
  std::vector<size_t> gen(nb);
  for (size_t a(0); a < nb; ++a) gen[a] = a;
  const double hgen(hof(gen, false));
  // YFS: BPMC_TRACE n - the first n events: the Born of the generated
  // assignment at the reconstructed generator point against m_born (the
  // generator's own Born; must be 1), and G.
  static const int bptr(ATOOLS::Settings::GetMainSettings()["YFS"]
                        ["BPMC_TRACE"].SetDefault(0).Get<int>());
  static int nbptr(0);
  if (bptr > 0 && nbptr < bptr) {
    Vec4D_Vector p;
    if (GeneratorPointOf(m_flavs, fin0, zs, p) && m_born > 0.) {
      ++nbptr;
      std::cerr<<"@@@ BPMC Bgen/m_born="<<BornME2At(p)/m_born
               <<" sqrt_sp="<<(p[0]+p[1]).Mass()<<" jswap="<<jswap<<std::endl;
    }
  }
  double sum(0.);
  int n(0);
  if (hgen > 0.) {
    if (jswap >= 0) {
      if ((size_t)jswap < isr.size())
        for (size_t a(0); a < nb; ++a) {
          std::vector<size_t> pick(gen);
          pick[a] = nb + jswap;
          const double h(hof(pick, true));
          if (h > 0.) { sum += h; ++n; }
        }
    } else {
      // ISR photons that pass as a Born photon in at least one single swap
      std::vector<size_t> ok;
      for (size_t j(0); j < isr.size(); ++j)
        for (size_t a(0); a < nb; ++a) {
          std::vector<size_t> pick(gen);
          pick[a] = nb + j;
          Vec4D_Vector fin(fin0), p;
          for (size_t b(0); b < nb; ++b) fin[slot[b] - 2] = cand[pick[b]];
          if (GeneratorPointOf(m_flavs, fin, zs, p) && (!sel || sel->Trigger(p))) {
            ok.push_back(nb + j); break; }
        }
      if (!ok.empty()) {
        // every ordered choice of nb distinct photons from Born + ok, with
        // the Born photon legs kept in slot order (identical photons: one
        // ordering per set, the first in ascending candidate index)
        std::vector<size_t> pool(gen);
        pool.insert(pool.end(), ok.begin(), ok.end());
        std::vector<size_t> idx(nb);
        std::function<void(size_t, size_t)> rec = [&](size_t d, size_t from) {
          if (d == nb) {
            std::vector<size_t> pick(nb);
            for (size_t a(0); a < nb; ++a) pick[a] = pool[idx[a]];
            if (pick == gen) return;
            const double h(hof(pick, true));
            if (h > 0.) { sum += h; ++n; }
            return;
          }
          for (size_t i(from); i < pool.size(); ++i) { idx[d] = i; rec(d + 1, i + 1); }
        };
        rec(0, 0);
      }
    }
  }
  if (sel) sel->SetResult(res0);
  if (nalt) *nalt = n;
  if (!(hgen > 0.) || n == 0) return 1.;
  const double G(1. + sum/hgen);
  return IsBad(G) ? 1. : G;
}

double NLO_Base::BornPhotonChannelWeight(const Vec4D &k) {
  m_bpmc_lastG = 1.; m_bpmc_lastn = 0;
  static const bornphotonmc::code mode(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["REAL_BORN_PHOTON_MULTICHANNEL"].SetDefault(bornphotonmc::on)
      .Get<bornphotonmc::code>());
  if (mode == bornphotonmc::off || !p_dipoles || !p_dipoles->HasDipoleII()) return 1.;
  std::vector<Vec4D> isr;
  int j(-1);
  for (const YFS::Photon &g : m_photons) {
    if (g.IsFSR()) continue;
    if (g.K() == k) j = (int)isr.size();
    isr.push_back(g.K());
  }
  if (j < 0) return 1.;
  const double G(BornPhotonChannelSum(m_plab, isr, p_dipoles->GetDipoleII(), j,
                                      &m_bpmc_lastn));
  m_bpmc_lastG = G;
  return mode == bornphotonmc::on ? 1./G : 1.;   // report: compute and report only
}

const YFS::Photon *NLO_Base::FindPhoton(const Vec4D &k) const {
  for (const YFS::Photon &g : m_photons)
    if (g.K() == k) return &g;
  return nullptr;
}

bool NLO_Base::PhotonIsFSR(const Vec4D &k) const {
  const YFS::Photon *g(FindPhoton(k));
  return g != nullptr && g->IsFSR();
}
