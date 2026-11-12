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
  Whether the Born has a space-like exchange line (a current with one initial
  leg and part of the final state), the only lines COMIX::Amplitude::
  SetPropShifts moves. Cached per Born process. For YFS: REAL_COMBINE -1.
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
    leg.push_back(std::make_pair(LegUnion(LegBit(i), PHASIC::Process_Base::s_propshiftleg),
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
  m_real_ifi_soft = 0.;
  if (m_coll_real)
    return p_dipoles->CalculateEEX() * m_born;
  if (!m_realtool)
    return 0;
  double real(0);
  m_real_hard1 = 0.;
  m_real_hard2 = 0.;
  m_ifi_prod = 1.;
  double prodw(1.);
  m_rvinfo.clear();
  m_rrphot.clear();
  m_prodfac.assign(m_photons.size(), 1.);
  for (YFS::Photon &g : m_photons) {
    const Vec4D k(g.K());
    // RV_MODE 1: one entry per photon, filled below if the real is evaluated
    if (RVIsRemainder()) m_rvinfo.push_back(RVPointInfo());
    m_lastrvinfo = RVPointInfo();
    if (RRFromSingles()) m_rrphot.push_back(RRPhotonInfo());   // factor 1 if skipped
    m_lastrrinfo = RRPhotonInfo();
    // the photons IFIRealRatio would reweight, if CalculateReal(k) reaches it
    const bool ifi_eligible(m_ifireal && m_submode == submode::global
                            && k.E() > p_dipoles->IFIOmega());
    if (ifi_eligible) ++m_ifi_eligible;
    const double phemin(PhotonEminNLO());
    if (phemin>0.0 && k.E()<phemin) {
      if (ifi_eligible) ++m_ifi_skip_emin;
      g.m_beta10 = 0.;
      continue;
    }
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
    const double contrib(CalculateReal(k));
    if (ifi_eligible && !m_ifi_applied) CountIFISkip();
    if (RVIsRemainder()) m_rvinfo.back() = m_lastrvinfo;
    if (RRFromSingles()) m_rrphot.back() = m_lastrrinfo;
    real += contrib;
    // YFS: REAL_COMBINE 1: this photon's factor, its bracket plus the IF
    // interference ratio m_ifi_prod picked up for it (1 without IFI_Real)
    if (m_born != 0.) {
      const double rif(m_ifireal && ifi_before != 0. && !IsBad(m_ifi_prod)
                       ? m_ifi_prod/ifi_before : 1.);
      prodw *= contrib/m_born + rif;
      if (RRFromSingles()) m_rrphot.back().factor = contrib/m_born + rif;
      const size_t ig(&g - &m_photons[0]);
      if (ig < m_prodfac.size()) m_prodfac[ig] = contrib/m_born + rif;
      if (m_acraic_compare) (g.IsISR() ? m_dbg_facisr : m_dbg_facfsr) = contrib/m_born + rif;
    }
    if (m_check_real_sub == CHECK_REAL_SUB_SCATTER)
      RecordSubScatter(k, contrib, g.IsISR() ? "realISR" : "realFSR", m_eikeex);
    g.m_beta10 = contrib;
  }
  HardestBetas(m_photons, [](const YFS::Photon &g) { return g.beta10(); },
               m_real_hard1, m_real_hard2);
  if (m_ifireal && !IsBad(m_ifi_prod)) m_real_ifi_soft = m_born*(m_ifi_prod - 1.);
  real += m_real_ifi_soft;
  /*
    YFS: REAL_COMBINE: sum (0) is 1 + sum_j delta_j, product (1, default)
    prod_j (1 + delta_j), delta_j = beta_1(k_j)/(S~_j B) plus that photon's IF
    ratio - 1. The two agree at O(alpha), the product adds the factorised
    beta_2 ~ beta_1 beta_1/beta_0. auto (-1): product with an exchange line.
  */
  const bool useprod(RealCombineMode() == realcombine::product
                     || (RealCombineMode() == realcombine::automatic
                         && BornHasExchangeLine(p_bornproc, m_plab)));
  // RV_MODE 1's real-virtual is a remainder on top of this product (it
  // supplies neither beta_2 nor v x real), so the product stays on with it
  if (useprod && m_born != 0. && !IsBad(prodw)
      && (!m_rrtool || RRFromSingles()) && (!m_realvirt || RVIsRemainder()))
    real = m_born*(prodw - 1.);
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
  YFS: SUB8_FSR_LEGS 3: the density of the emission the (n+1)-body point
  describes. REAL_FSR_MAP 2 rebuilds the radiating pair at Q_D - k_j, so the
  density is S~(point legs; k_j) times the point's flux Q_D^2/(Q_D - k_j)^2,
  equal to GeneratorFSRDensity for one photon. kpt is the photon at the
  point, pt the point's momenta.
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
  m_real_exit = real_exit::completed;
  m_ifi_applied = false;

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
    m_real_exit = real_exit::momentum_conservation;
    return 0;
  }
  if ((rt.p[2] + rt.p[3]).Mass() < massmin)
    massmin = (rt.p[2] + rt.p[3]).Mass();
  if (!RealMEUsable(rt)) {
    m_real_exit = real_exit::unusable_me;
    return 0;
  }

  PointDenominator(rt);
  if (RealSubEik() == realsubeik::multichannel && !PhotonIsFSR(rt.kk))
    MultichannelISRPhoton(rt);
  if (RealSubEik() == realsubeik::multichannel_prefsr && !PhotonIsFSR(rt.kk)
      && p_bornproc != nullptr && Sub8Part() != 1)
    PreFSRMultichannelISRPhoton(rt);
  AssembleBeta1(rt);

  /*
    YFS: REAL_ALPHA0 (default true): the real is Comix's at the model's alpha
    and the subtraction is raised to match (/m_rescale_alpha), while beta_0
    and the eikonals carry alpha(0). true rescales beta_1 by m_rescale_alpha,
    the whole real correction at alpha(0) as in ACRAIC
    (Acraic_Base::ComixPhotonCoupling). false keeps the alpha_model remainder.
  */
  if (RealAlpha0()) rt.tot *= m_rescale_alpha;

  IFIRealRatio(rt);

  msg_Debugging() << METHOD << " submode=" << m_submode
                  << " r*flux=" << rt.r*rt.flux
                  << " sub=" << rt.subloc * m_born / m_rescale_alpha
                  << " tot=" << rt.tot << "\n";

  if (m_no_subtraction) {
    msg_Debugging() << METHOD << " no_subtraction: returning r/subloc=" << rt.r/rt.subloc << "\n";
    return rt.r / rt.subloc;
  }
  if (IsBad(rt.tot)) ReportBadReal(rt);
  TrackRealAverage(rt.tot);

  if (raw) {
    double rawval = rt.r * rt.flux - rt.subloc * m_born / m_rescale_alpha;
    msg_Debugging() << METHOD << " raw: returning " << rawval << "\n";
    return rawval;
  }

  if (m_acraic_compare) {
    m_dbg_r = rt.r; m_dbg_flux = rt.flux; m_dbg_subloc = rt.subloc;
    m_dbg_subb = rt.subb; m_dbg_tot = rt.tot;
    const Vec4D dk(rt.k - rt.kk);
    m_dbg_dk = Max(Max(dabs(dk[0]),dabs(dk[1])),Max(dabs(dk[2]),dabs(dk[3])));
    m_dbg_dp = 0.;
    for (size_t i(0); i < m_reallab.size() && i < rt.p.size(); ++i) {
      const Vec4D d(rt.p[i] - m_reallab[i]);
      m_dbg_dp = Max(m_dbg_dp, Max(Max(dabs(d[0]),dabs(d[1])),Max(dabs(d[2]),dabs(d[3]))));
    }
    if (PhotonIsFSR(rt.kk) && rt.p.size() > 3) {
      Vec4D X(rt.p.back());
      for (size_t i(2); i + 1 < rt.p.size(); ++i) X += rt.p[i];
      m_dbg_x2pt = X.Abs2();
      m_dbg_x2gen = PreFSRSystem().Abs2();
    }
    else if (!PhotonIsFSR(rt.kk) && rt.p.size() > 3) {
      Vec4D F;
      for (size_t i(2); i + 1 < rt.p.size(); ++i) F += rt.p[i];
      m_dbg_x2isr = F.Abs2();
      if (!(m_dbg_x2gen > 0.)) m_dbg_x2gen = PreFSRSystem().Abs2();
    }
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
  const Vec4D &kk(rt.kk);
  CheckMasses(p, 1);

  rt.pp = p;
  rt.pp.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, rt.pp, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, rt.pp, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, rt.pp, m_plab);

  double r = p_real->Calc_R(p) / norm;
  /*
    YFS: REAL_BORN_PHOTON_SYM. The real carries 1/(N_g+1)! for the (n+1)-photon
    final state, the Born 1/N_g!, so r/m_born is 1/(N_g+1) of |M_1|^2/|M_0|^2,
    while the crude subb is not divided. Multiplying r back by N_g+1 puts
    every subtraction mode in one normalisation (1 without Born photons).
  */
  r *= BornPhotonSym(1);
  /*
    YFS: REAL_BORN_PHOTON_MULTICHANNEL: r is reduced to the share of this
    photon's generated assignment, h_gen/sum_b h_b (BornPhotonChannelWeight).
    Each channel keeps its own beta_0 and subtraction.
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
    Define_Dipoles::CalculateFlux gives every photon the initial-state flux
    (Q_X - k)^2/Q_X^2 when ISR and FSR are both on, which a final-state photon
    must not have: the Born scale s' is untouched by FSR. REAL_FSR_FLUX
    selects the replacement (RealFSRFlux() in NLO_Base_Internal.H).
  */
  const realfsrflux::code fsrflux(RealFSRFlux());
  if (fsrflux != realfsrflux::event_flux && PhotonIsFSR(kk)) {
    if (fsrflux == realfsrflux::no_flux) flux = 1.;
    else if (fsrflux == realfsrflux::flux_squared) flux *= flux;  // (m_ff^2/s')^2
    /*
      own_dipole: the flux (Q_D - k)^2/Q_D^2 of the photon's OWN pair, Q_D
      its pre-emission momentum. CalculateFlux uses the whole ISR-reduced
      final state, the same for one resonant pair but not for two, while the
      recoil Jacobian of the FSR generator (FSR::RescalePhotons, m_yy) is
      that of the pair in its own frame.
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
  rt.flux0 = rt.flux; rt.subloc0 = rt.subloc; rt.subb0 = rt.subb;
  /*
    Which eikonal beta_1 subtracts (YFS: REAL_SUB_EIK, default
    multichannel_prefsr (8)):
      coherent_point (0): the coherent one at the mapped post-emission point.
      crude (1): the event's crude S~_II + S~_FF on the pre-emission legs, the
        generation density, so the weight is r flux/(S~ B).
      coherent_born (2): the coherent eikonal on the pre-emission legs, those
        of the form factor exponent.
      3 to 8: see the branches below and the SubEik* helpers.
  */
  static const realsubeik::code subeik(RealSubEik());
  if (subeik == realsubeik::crude) rt.subloc = rt.subb;
  else if (subeik == realsubeik::coherent_born) rt.subloc = p_dipoles->CalculateRealSub(rt.kk);
  /*
    crude_fsr (3): the event's crude for FINAL-state photons only. For a hard
    FSR photon the coherent eikonal at the (n+1)-body point is smaller than the
    event's crude, which adds a spurious (1 - S_pt/S~_ev) to the bracket.
    ISR photons keep their reduced-point subtraction (m_map_reduced).
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
  REAL_SUB_EIK post_emission_fsr (4): for FINAL-state photons, subtraction and
  denominator on the POST-emission legs. FSR::Weight (m_wt2/m_yy) already
  carries the Jacobian to the physical momenta, so the generation density is
  the eikonal on the physical legs, as in ACRAIC's rho_crude. On the
  pre-emission legs a soft photon next to a kicked lepton gets an O(1)
  bracket inside the collinear cone. One photon: r/(S~_post B) = ACRAIC.
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
  eikonal |J_II + J_FF|^2 on the pre-emission legs. The incoherent crude (3)
  leaves the soft interference 2 Re(J_II J_FF^*) in the bracket while
  m_ifi_prod supplies it multiplicatively. With the coherent eikonal the soft
  brackets vanish and the one-photon weight stays r/(S~_cr B).
*/
void NLO_Base::SubEikCoherentBornFSR(RealTerms &rt)
{
  // The pieces m_ifi_prod is built from (IFIRealRatio): crude plus IF
  // interference on the generation legs. Not CalculateRealSub(kk), which
  // sits on the post-emission dipole momenta.
  const double ifg(p_dipoles->CalculateRealSubIF(rt.kk));
  if (!IsBad(ifg) && !IsZero(rt.subb)) rt.subloc = rt.subb + ifg;
}

/*
  REAL_SUB_EIK assignment_born_fsr (6): mode 3 with the crude carrying the
  Born of each ASSIGNMENT. A hard FSR-labelled photon at wide angle may really
  be the ISR photon bringing an off-peak pair onto the Z peak, so r is large
  while m_born is tiny. Subtraction and denominator become S~_II B_I/B +
  S~_FF, B_I the Born at the scaled ISR-reduced point of this photon.
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
  REAL_SUB_EIK multichannel (7), FSR-labelled photon: the generator reaches a
  final state with the photon labelled ISR or FSR, g_I = S~_II B_I/flux_I,
  g_F = S~_FF B_F, and the unbiased weight divides r by g_I + g_F whatever
  the label. For an FSR-labelled photon B = B_F, so subb = S~_FF + S~_II
  (B_I/B)/flux_I, B_I the Born at (P - k)^2, flux_I = (P - k)^2/P^2. subloc
  adds the IF interference on the same legs (as mode 5). No flux on r. The
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
  REAL_SUB_EIK multichannel_prefsr (8), FSR-labelled photon: mode 7 with the
  alternative channel built from the event's PRE-FSR final state Q_pre (NLO
  m_plab). An FSR-labelled photon's ISR alternative has its Born at
  (Q_pre - k)^2 with flux (Q_pre - k)^2/Q_pre^2, an ISR-labelled photon's FSR
  alternative at (Q_pre + k)^2, Born ratios taken between points of the same
  construction (PreFSRBornRatio). subloc is scaled with the denominator, so
  with m_ifi_prod the one-photon weight is r/(g_I + g_F) and the soft bracket
  is unchanged. The base is mode 5, flux 1. The ISR-labelled half is
  PreFSRMultichannelISRPhoton.
*/
void NLO_Base::SubEikPreFSRMultichannelFSR(RealTerms &rt)
{
  const Vec4D &kk(rt.kk);
  double &subloc(rt.subloc), &subb(rt.subb);
  const double sII(GenerationEikonalII(*p_dipoles, kk));
  /*
    YFS: SUB8_FSR_LEGS (default point_density (3)). generation (0): the
    pre-emission eikonal. post_emission (1): the FSR density on the
    POST-emission legs times the pair flux, S~_FF(q_1', q_2'; k)
    (q_D + k)^2/q_D^2, as the FSR generator produces and ACRAIC's rho_crude
    uses. The pre-emission eikonal has a narrower dead cone (m/E_pre).
  */
  static const sub8fsrlegs::code fsrlegs(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["SUB8_FSR_LEGS"].SetDefault(sub8fsrlegs::point_density).Get<sub8fsrlegs::code>());
  /*
    event_density (2) and point_density (3): the FSR density the generator
    really has, post-emission eikonal times pair flux, so the one-photon
    weight is r/(g B). Unlike 1, the SOFT half keeps the generation-leg scale
    (gsoft below), the point's coherent current, which a post/pre eikonal
    ratio would spoil for a lepton kicked by a companion.
    event_density: GeneratorFSRDensity, S~ on the event's post-emission legs
       times the dipole's flux F_D shared among its photons. With companions
       it describes a different emission from the one |M_1|^2 is evaluated for.
    point_density: PointFSRDensity, S~ on the legs of the (n+1)-body point
       times its own single-emission flux. Equal to event_density at one photon.
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
  double gsoft(densityscale ? sFFgen + isrchannel : g);
  /*
    YFS: SUB8_SOFT_PAIR_FLUX (default true): g carries each dipole's pair flux
    (q_D + k)^2/q_D^2 = 1 + O(x), gsoft does not, leaving soft brackets at
    -x r_IF. true gives gsoft the same flux, dipole by dipole, on the
    generation-leg eikonals.
  */
  static const bool softflux(ATOOLS::Settings::GetMainSettings()["YFS"]
                             ["SUB8_SOFT_PAIR_FLUX"].SetDefault(true).Get<bool>());
  if (softflux && densityscale && rt.pp.size() == m_flavs.size()) {
    double sgf(0.);
    for (auto &D : p_dipoles->GetDipoleFF()) {
      const int l(D.Left()), r(D.Right());
      double F(1.);
      if (l >= 2 && r >= 2 && l < (int)rt.pp.size() && r < (int)rt.pp.size()) {
        const Vec4D q(rt.pp[l] + rt.pp[r]);
        if (q.Abs2() > 0.) F = (q + rt.k).Abs2()/q.Abs2();
      }
      sgf += F*D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1));
    }
    if (sgf > 0. && !IsBad(sgf)) gsoft = sgf + isrchannel;
  }
  /*
    Which subtraction (YFS: SUB8_FSRSUB, default blend (4)):
    generation (0): mode 5's on the generation legs, scaled with g. Gives the
       one-photon weight r/(g_I + g_F), but r's soft limit sits on the point's
       legs, so soft FSR photons next to a hard ISR or FSR companion keep
       O(1) brackets.
    point (1): the coherent current at the point, scaled with g. Soft brackets
       vanish, the hard-photon identity does not hold.
    blend (4): point for soft photons, generation for hard, with
       w = 1/(1 + (y/y0)^2), y = 2 k.Q_pre/Q_pre^2, y0 = YFS: SUB8_Y0.
  */
  static const sub8fsrsub::code fsub(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["SUB8_FSRSUB"].SetDefault(sub8fsrsub::blend).Get<sub8fsrsub::code>());
  if (g > 0. && !IsBad(g) && subb > 0.) {
    const double cru(subb);
    const double sgen(g*(1. + (IsBad(ifg) ? 0. : ifg/cru)));
    if (fsub == sub8fsrsub::point) subloc *= gsoft/cru;
    else if (fsub == sub8fsrsub::blend) {
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
  REAL_FSR_FLUX on_crude (3): the crude a FINAL-state photon is divided by
  carries F = (q + k)^2/q^2 on its final-state part, q the post-emission pair,
  and r carries no flux. This is the crude ACRAIC divides by.
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
  momenta of the point, incoherently, like CalculateRealSubEEX for the event.
  Only used when the mapped point is not the event (m_map_reduced).

  beta_1 at the reduced point has the photon's collinear structure at a
  different scale from the event, so the residual r flux - S~coh B is divided
  by the crude of the point: positive, and |J_II + J_FF|^2 <= 2 (|J_II|^2 +
  |J_FF|^2), so it stays bounded where residual/S~coh does. The coherent
  S~coh would put the zeros of the II-FF interference into the denominator.
  Crude(point) ~ crude(event)/x^2, so residual/crude(point) is the
  event-normalised beta_1/S~crude, and the n = 1 formula with subb the
  event's crude is the same object at x = 1.
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
  const double flux(rt.flux);
  double bF(1.), s2(1.);
  if (ReducedBornAlt(kk, +1, bF, s2)) {
    const double f(flux*bF);
    const double sFFpt(m_map_reduced ? PointEikonalFF(*p_nlodipoles, k)
                                     : GenerationEikonalFF(*p_dipoles, kk));
    if (!IsBad(f) && !IsBad(sFFpt)) {
      subb_loc += (f - 1.)*sFFpt;
      subloc   += (f - 1.)*sFFpt;
      if (!m_map_reduced) subb += (f - 1.)*sFFpt;
    }
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
    (final legs of the point + its photon)^2: at n >= 2 the scaled point
    (MapMomentaScaled) boosts the pair rigidly, moving both S~_FF of the point
    and (Q' + k')^2, and r keeps its share of that collinear structure. At one
    photon the point is the event.
  */
  const Vec4D Q(PreFSRSystem());
  Vec4D R(Q + kk);
  if (m_map_reduced && rt.pp.size() == m_flavs.size() && Q.Abs2() > 0.) {
    Vec4D Rp(k);
    for (size_t i(2); i < rt.pp.size(); ++i) Rp += rt.pp[i];
    const double r2(Rp.Abs2());
    if (r2 > 0.) R = sqrt(r2/Q.Abs2())*Q;
  }
  double bF(0.), f(0.);
  if (!PreFSRBornRatio(R, bF)) bF = 0.;
  f = rt.flux*bF;
  /*
    YFS: SUB8_ISR_FSR_CHANNEL_FLUX (default true): the FSR generator produces
    a photon with the post-emission eikonal times the pair flux (q + k)^2/q^2,
    as ACRAIC's rho_crude. true multiplies the FSR channel term by R^2/Q^2,
    the pair flux of the FSR interpretation (pre-FSR system R = Q + k, pair
    Q), on top of the ISR photon's own flux.
  */
  static const bool fsrchflux(ATOOLS::Settings::GetMainSettings()["YFS"]
                              ["SUB8_ISR_FSR_CHANNEL_FLUX"].SetDefault(true).Get<bool>());
  if (fsrchflux && Q.Abs2() > 0. && R.Abs2() > 0.) f *= R.Abs2()/Q.Abs2();
  const double sFFpt(m_map_reduced ? PointEikonalFF(*p_nlodipoles, k)
                                   : GenerationEikonalFF(*p_dipoles, kk));
  // YFS: SUB8_ISR_REDUCED (default true): also rescale when the point is
  // the reduced one (m_map_reduced)
  static const bool isrred(ATOOLS::Settings::GetMainSettings()["YFS"]
                           ["SUB8_ISR_REDUCED"].SetDefault(true).Get<bool>());
  const double nb(subb_loc + (f - 1.)*sFFpt);
  if ((isrred || !m_map_reduced) && !IsBad(nb) && nb > 0. && subb_loc > 0.) {
    const double scale(nb/subb_loc);
    subb_loc = nb;
    subloc  *= scale;
    if (!m_map_reduced) subb *= scale;
  }
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
      m_ifi_applied = true;
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

void NLO_Base::CountIFISkip()
{
  switch (m_real_exit) {
  case real_exit::momentum_conservation: ++m_ifi_skip_momcons;  break;
  case real_exit::unusable_me:           ++m_ifi_skip_unusable; break;
  case real_exit::completed:             ++m_ifi_skip_other;    break;
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
  if (RVIsRemainder() && m_born != 0. && m_lastdenom != 0.) {
    m_lastrvinfo.ok = true;
    m_lastrvinfo.p = rt.p;
    m_lastrvinfo.rho = (RealAlpha0() ? m_rescale_alpha : 1.)*r*flux/(m_lastdenom*m_born);
  }
  if (RRFromSingles() && m_born != 0. && m_lastdenom != 0. && !IsBad(tot)) {
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

bool NLO_Base::RealMEForSubset(unsigned mask, double &me,
                               std::vector<double> &subloc)
{
  const size_t n(LegCount(mask));
  YFS::Real_Correction *prov(RealProvider(n));
  if (prov == NULL) {
    msg_Error()<<METHOD<<"(): no "<<n<<"-photon real ME provider. The subset "
               <<"recursion is general, the matrix element is not."<<std::endl;
    return false;
  }

  Vec4D_Vector ks;
  std::vector<size_t> idx;
  for (size_t i(0); i < m_photons.size(); ++i)
    if (HasLeg(mask, i)) { ks.push_back(m_photons[i].K()); idx.push_back(i); }

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
  if (LegCount(mask) == 1) {
    beta = CalculateReal(m_photons[LowestLeg(mask)].K(), /*raw*/true);
  } else {
    double me(0.);
    std::vector<double> subloc;
    if (!RealMEForSubset(mask, me, subloc)) { m_rawbeta[mask] = 0.; return 0.; }
    beta = me;
    // every non-empty subset S of `mask`: eikonal product over S times the
    // residual of what is left.
    for (unsigned s(mask); s; s = NextSubset(s, mask)) {
      double eik(1.);
      for (unsigned t(s); t; t = WithoutLowestLeg(t)) eik *= subloc[LowestLeg(t)];
      beta -= eik * RawBeta(RemoveLegs(mask, s));
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
  const size_t n(LegCount(mask));

  /*
    Cut against photons the resummation has ALREADY accounted for: a photon
    below the YFS infrared cutoff belongs to the exponentiated form factor,
    and beta_n for it is a numerically unstable cancellation with no physical
    content. A set is dropped if ANY of its photons is unresolved. Applied
    here as well as in the subset loop so that no caller bypasses it.
    Controlled by NLO_PHOTON_EMIN (multiples of the cutoff energy, 0 = off).
  */
  const double phemin(PhotonEminNLO());
  if (phemin > 0.) {
    for (size_t i(0); i < m_photons.size(); ++i) {
      if (HasLeg(mask, i) && m_photons[i].K().E() < phemin) {
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
    if (HasLeg(mask, i)) crude *= p_dipoles->CalculateRealSubEEX(m_photons[i].K());
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
  Sum beta_n over every n-photon subset of this event's photons (n=2 is the
  pair loop of CalculateRealReal()). A subset is dropped if ANY of its
  photons is below the resummed threshold. Cost is C(m,n) subsets, but
  m_rawbeta is keyed on a mask over m_photons and lives for the event, so
  each residual is computed once. HardestBetas is the seed for a truncation
  policy beyond the resummed-photon cut.
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
      mask = AddLeg(mask, c[j]);
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
    The Born point of the ISR assignment: invariant s'' = (P - k)^2, P the
    event's beam sum, same scattering angles. Built in the rest frame of
    R = P - k with beams back to back along the event's beam axis and the
    final legs rescaled by one xi to energy sqrt(s''), then boosted onto R.
    (MapMomentaScaled's balance check assumes an ISR-labelled photon.)
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
    The Born point of the alternative assignment: invariant s'' = R^2,
    R = P - k (sign < 0) or P + k, P the event's beam sum, same scattering
    angles. Built in the rest frame of R with beams back to back along the
    event's beam axis and the final legs rescaled by one xi to energy
    sqrt(s''), then boosted onto R.
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
  if (!GeneratorReaches(p)) return false;        // outside the generated region: no channel
  const double me2(BornME2At(p));
  if (!(me2 >= 0.)) return fail(2, "Born |M|^2 at the alternative point unavailable");
  bratio = me2/m_me2pre;
  return !IsBad(bratio);
}

/*
  Whether the generator can produce an event through the Born point p: its
  ISR v within the generated range (m_bpmc_vmax, from the selector's smallest
  invariant mass) and p inside the Born selectors. A channel the generator
  never runs has no density and must not enter the multichannel sum g.
  YFS: SUB8_REACHABLE_CHANNELS false keeps every channel.
*/
bool NLO_Base::GeneratorReaches(const Vec4D_Vector &p) {
  static const bool reachable(ATOOLS::Settings::GetMainSettings()["YFS"]
                              ["SUB8_REACHABLE_CHANNELS"].SetDefault(true).Get<bool>());
  if (!reachable || p.size() < 2) return true;
  const double s(m_bpmc_s > 0. ? m_bpmc_s : m_s);
  if (s > 0. && 1. - (p[0] + p[1]).Abs2()/s > m_bpmc_vmax) return false;
  PHASIC::Combined_Selector *sel(p_bornproc ? p_bornproc->Selector() : nullptr);
  if (sel == nullptr) return true;
  const int res0(sel->Result());
  const bool pass(sel->Trigger(p));
  sel->SetResult(res0);
  return pass;
}

double NLO_Base::BornPhotonSym(size_t nextra) const {
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["REAL_BORN_PHOTON_SYM"].SetDefault(true).Get<bool>());
  if (!on) return 1.;
  size_t nb(0);
  for (const Flavour &f : m_flavs) if (f.IsPhoton()) ++nb;
  /*
    C(nb + nextra, nextra), not (nb + nextra)!/nb!: YFS sums its photons as an
    unordered set, so only the mixing of the extra photons with the Born's
    photons needs a factor. The one-photon value is nb + 1.
  */
  double sym(1.);
  for (size_t j(1); j <= nextra; ++j) sym *= double(nb + j)/double(j);
  return sym;
}

/*
  Born-photon multichannel (YFS: REAL_BORN_PHOTON_MULTICHANNEL: off (0), on
  (1, default), report (2) computes G for WEIGHT_PROBE without applying it).

  For a Born with final-state photons (e+e- -> gamma gamma) the generator
  picks the Born pair, then adds ISR photons. Once an ISR photon passes the
  generation cuts, the same final state is also reached with that photon in
  the Born pair and a Born photon as the ISR one, each assignment with its
  own crude density
      h = prod_{ISR} S~_II(k) * B(pair at its reduced point) * s/s'.
  Each channel alone integrates the exact real to |M_{n+1}|^2 (BornPhotonSym),
  so with several open the real is counted once per channel. The multichannel
  weight divides the real by the sum of the densities of all open channels,
  i.e. multiplies the generated channel's real by h_gen/sum_b h_b = 1/G.
  G = 1 when no other assignment passes the cuts, so soft photons and
  processes without Born photons are untouched.
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

// YFS: BORN_PHOTON_CHANNELS: isr (ISR photons, S~_II) or stage (every photon, the caller's density)
namespace {
  YFS::bornphotonchannels::code BornPhotonChannels() {
    static const YFS::bornphotonchannels::code mode(
        ATOOLS::Settings::GetMainSettings()["YFS"]["BORN_PHOTON_CHANNELS"]
        .SetDefault(YFS::bornphotonchannels::isr).Get<YFS::bornphotonchannels::code>());
    return mode;
  }
}

double NLO_Base::BornPhotonChannelSum(const Vec4D_Vector &lab,
                                      const std::vector<Vec4D> &isr,
                                      YFS::Dipole &dII, int jswap, int *nalt,
                                      const ChannelDensity *density) {
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
  const bool allchannels(BornPhotonChannels() != bornphotonchannels::isr && p_dipoles);
  auto eik = [&](const Vec4D &k) {
    double e(dII.Eikonal(k, b0, b1));
    if (allchannels)
      for (auto &D : p_dipoles->GetDipoleFF()) {
        if (!D.IsResonance()) continue;
        const size_t l(D.Left()), r(D.Right());
        if (l < lab.size() && r < lab.size()) e += D.Eikonal(k, lab[l], lab[r]);
      }
    return e;
  };
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
    // the caller's density of this assignment (ACRAIC's crude)
    if (density) {
      const double h((*density)(slot, pick));
      return (IsBad(h) || !(h > 0.)) ? -1. : h;
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
  // the generated assignment's density, taken only once another assignment
  // is a channel (with a density it is one more crude per event)
  double hgen(-1.);
  bool hgenknown(false);
  auto gendensity = [&]() {
    if (!hgenknown) { hgen = hof(gen, false); hgenknown = true; }
    return hgen;
  };
  double sum(0.);
  int n(0);
  if (jswap >= 0) {
    if ((size_t)jswap < isr.size() && gendensity() > 0.)
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
    if (!ok.empty() && gendensity() > 0.) {
      // ordered choices of nb distinct photons, one ordering per set
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
  if (sel) sel->SetResult(res0);
  if (nalt) *nalt = n;
  if (!(hgen > 0.) || n == 0) return 1.;
  const double G(1. + sum/hgen);
  return IsBad(G) ? 1. : G;
}

double NLO_Base::BornPhotonChannelWeight(const Vec4D &k) {
  m_bpmc_lastn = 0;
  static const bornphotonmc::code mode(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["REAL_BORN_PHOTON_MULTICHANNEL"].SetDefault(bornphotonmc::on)
      .Get<bornphotonmc::code>());
  if (mode == bornphotonmc::off || !p_dipoles || !p_dipoles->HasDipoleII()) return 1.;
  // YFS: BORN_PHOTON_CHANNELS stage (no ACRAIC crude here): FSR-labelled
  // photons are alternatives too, on the post-emission legs
  const bool allchannels(BornPhotonChannels() != bornphotonchannels::isr);
  std::vector<Vec4D> isr;
  int j(-1);
  for (const YFS::Photon &g : m_photons) {
    if (g.IsFSR() && !allchannels) continue;
    if (g.K() == k) j = (int)isr.size();
    isr.push_back(g.K());
  }
  if (j < 0) return 1.;
  const double G(BornPhotonChannelSum(allchannels ? m_postlab : m_plab, isr,
                                      p_dipoles->GetDipoleII(), j, &m_bpmc_lastn));
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
