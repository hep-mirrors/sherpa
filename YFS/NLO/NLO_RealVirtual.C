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

/*
  YFS: RV_PROBE (diagnostic, default 0): one line per real-virtual photon on
  std::cerr, "@@@ RVPROBE", with the loop over the tree lt = V_fin/T of the
  last loop call, the YFS dim-reg subtraction on the event's Born legs (the
  dipoles CalculateRealVirtual builds, Born momenta m_plab) and on the
  (n+1)-body point's own charged legs, the Born's v = V_sub/B, and the
  one-loop pole coefficients of each. dv_ev and dv_pt are v_{n+1} - v with
  the two subtractions: the object that must vanish in the soft limit.
*/
static bool RVProbeOn() {
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PROBE"]
                       .SetDefault(0).Get<int>() != 0);
  return on;
}

bool NLO_Base::RVRemainder() const { return m_realvirt && RVMode() == rvmode::remainder; }

/*
  RV_MODE 1: v_{n+1} at the (n+1)-body point p (photon last), in the units of
  v = V_sub/B: the provider's loop over its own tree at the Born virtual's
  mu, minus the Born virtual's subtraction function on the point's charged
  legs. The subtraction is built after the loop call, so that a dim-reg
  B~ reads this call's IR scale. Poles: the loop's 1/eps coefficient and the
  subtraction's must cancel on these legs; failures beyond 1e-6 relative are
  counted (m_rvPoleFail) and reported at the end.
*/
/*
  YFS: RV_LOOP_FRAME - the frame RV_MODE remainder hands the (n+1)-point loop
  to OpenLoops in (name or old integer).
  lab (0): the real's point as it is (the legacy RV's frame). Beams on z
     up to the rounding of the rotate-back (|p_T|/E ~ 1e-33): the
     hp_mode 1 offset of 0.54% on e+e- -> mu mu gamma, and photons within
     ~3e-6 rad of a beam lose precision on the axis.
  canonical (1, default): CanonicalBeamFrame, the beams' rest frame turned
     by a fixed generic rotation.
  tilted (2): the beams' rest frame with the beams exactly on z, then tilted
     by 1e-3 rad.
  beam_axis (3): the beams' rest frame with the beams exactly on z (no
     tilt): stable for every nu nu gamma point seen.
  No single frame is right everywhere for photons deep in the electron's dead
  cone; YFS: RV_LOOP_CHECK true compares two (see RealVirtualFactor).
*/

static Vec4D_Vector TiltedBeamFrame(const Vec4D_Vector &p)
{
  Vec4D_Vector q(BeamRestFrameOnZ(p));
  static const double c(cos(1e-3)), s(sin(1e-3));
  for (Vec4D &v : q) v = Vec4D(v[0], c*v[1] + s*v[3], v[2], -s*v[1] + c*v[3]);
  return q;
}

static Vec4D_Vector RVLoopPoint(const Vec4D_Vector &p, rvloopframe::code frame)
{
  switch (frame) {
  case rvloopframe::lab:       return p;
  case rvloopframe::tilted:    return TiltedBeamFrame(p);
  case rvloopframe::beam_axis: return BeamRestFrameOnZ(p);
  case rvloopframe::canonical: break;
  }
  return CanonicalBeamFrame(p);
}

/*
  YFS: RV_LOOP_CHECK (default false). true: for a photon within
  RV_LOOP_CHECK_ANGLE (default 1e-4 rad, ~10 m_e/E at the Z pole) of a beam,
  evaluate V_fin/T a second time, in frame tilted if RV_LOOP_FRAME is
  beam_axis and in frame beam_axis otherwise,
  and drop the photon's real-virtual (counted, reported at the end) when the
  two differ by more than RV_LOOP_CHECK_TOL (default 1e-3). The remainder of
  such photons is O(1e-3 .. 1e-2) of |M_1|^2 at most; OpenLoops' value there
  can be O(1e4) in one frame.
*/
static bool RVLoopCheckNeeded(const Vec4D_Vector &p)
{
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_LOOP_CHECK"]
                       .SetDefault(false).Get<bool>());
  static const double amax(ATOOLS::Settings::GetMainSettings()["YFS"]
                           ["RV_LOOP_CHECK_ANGLE"].SetDefault(1e-4).Get<double>());
  if (!on || p.size() < 3) return false;
  const Vec4D &k(p.back());
  for (size_t i(0); i < 2; ++i) {
    const double c(Vec3D(p[i])*Vec3D(k)/(Vec3D(p[i]).Abs()*Vec3D(k).Abs()));
    if (acos(Max(-1., Min(1., c))) < amax) return true;
  }
  return false;
}

bool NLO_Base::RealVirtualFactor(const Vec4D_Vector &pin, double &v1) {
  if (pin.size() != m_flavs.size() + 1) return false;
  const Vec4D_Vector p(RVLoopPoint(pin, RVLoopFrame()));
  double mur2(sqr(rpa->gen.Ecms()));
  if (m_looptool && p_virt->p_loop_me->IRscale() > 0.)
    mur2 = sqr(p_virt->p_loop_me->IRscale());
  double lt(0.);
  if (RVLoopCheckNeeded(pin)) {
    static const double tol(ATOOLS::Settings::GetMainSettings()["YFS"]
                            ["RV_LOOP_CHECK_TOL"].SetDefault(1e-3).Get<double>());
    double lt2(0.);
    const Vec4D_Vector q(RVLoopPoint(pin, RVLoopFrame() == rvloopframe::beam_axis
                                          ? rvloopframe::tilted : rvloopframe::beam_axis));
    if (!p_realvirt->LoopOverTree(pin, q, mur2, lt2)) return false;
    ++m_rvLoopChecked;
    if (!p_realvirt->LoopOverTree(pin, p, mur2, lt)) return false;
    if (!(std::abs(lt - lt2) <= tol)) { ++m_rvLoopUnstable; return false; }
  }
  else if (!p_realvirt->LoopOverTree(pin, p, mur2, lt)) return false;
  Vec4D_Vector pp(p);
  pp.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, pp, pp);
  p_nlodipoles->MakeDipoles(m_flavs, pp, pp);
  p_nlodipoles->MakeDipolesIF(m_flavs, pp, pp);
  p_nlodipoles->p_yfsFormFact->p_virt = p_realvirt->p_loop_me.get();
  const double sub(p_nlodipoles->CalculateVirtualSub());
  if (IsBad(sub)) return false;
  v1 = lt - sub/m_rescale_alpha;
  m_rvLastLt = lt;
  const double e1l(p_realvirt->p_loop_me->ME_E1()*p_realvirt->m_factor);
  const double e1s(p_nlodipoles->Get_E1()/m_rescale_alpha);
  if (m_dim_reg && std::abs(e1l + e1s) > 1e-6*std::max(std::abs(e1l), 1e-3))
    ++m_rvPoleFail;
  if (RVProbeOn()) {
    std::ostringstream o;
    o<<std::setprecision(10)<<"@@@ RVREM lt="<<lt<<" sub="<<sub/m_rescale_alpha
     <<" v1="<<v1<<" E1loop="<<e1l<<" E1sub="<<e1s<<" mu="<<sqrt(mur2);
    // YFS: RV_PROBE 2: also the evaluated (canonical-frame) point and the
    // event's Born point, at full precision, for the pyol referee
    static const int lvl(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PROBE"].Get<int>());
    if (lvl >= 2) {
      o<<std::setprecision(17)<<" PC=";
      for (const Vec4D &q : p) o<<q<<" ";
      o<<" PB=";
      for (const Vec4D &q : CanonicalBeamFrame(m_plab)) o<<q<<" ";
      o<<std::setprecision(10);
    }
    o<<"\n";   // one record per line (it ran into the next record before)
    std::cerr<<o.str();
  }
  return !IsBad(v1);
}

/*
  An (n+1)-body point with a soft photon of energy fraction x added to the
  Born point born: in the beams' rest frame the final state is shrunk in its
  own rest frame to the invariant mass (P - k)^2 (3-momenta scaled, masses
  kept) and boosted onto P - k; the beams are unchanged; photon last.
*/
static bool SoftPhotonPoint(const Vec4D_Vector &born, const double x,
                            Vec4D_Vector &pt)
{
  if (born.size() < 4) return false;
  const Vec4D_Vector q(BeamRestFrameOnZ(born));
  const Vec4D P(q[0] + q[1]);
  const double rs(P.Mass());
  const Vec3D n(Vec3D(0.36, 0.48, 0.8)/Vec3D(0.36, 0.48, 0.8).Abs());
  const double e(0.5*x*rs);
  const Vec4D k(e, e*n);
  const Vec4D Qp(P - k);
  const double Mp(Qp.Mass());
  // final state in its own rest frame (q is in P's rest frame already)
  Vec4D Qf;
  for (size_t i(2); i < q.size(); ++i) Qf += q[i];
  Poincare restf(Qf);
  Vec4D_Vector f;
  std::vector<double> m2;
  for (size_t i(2); i < q.size(); ++i) {
    Vec4D v(q[i]);
    restf.Boost(v);
    f.push_back(v);
    m2.push_back(Max(0., v.Abs2()));
  }
  double xi(1.);
  for (int it(0); it < 50; ++it) {
    double E(0.), dE(0.);
    for (size_t i(0); i < f.size(); ++i) {
      const double p2(Vec3D(f[i]).Sqr()), ei(sqrt(xi*xi*p2 + m2[i]));
      E += ei;
      dE += xi*p2/ei;
    }
    if (dE <= 0.) return false;
    const double step((E - Mp)/dE);
    xi -= step;
    if (std::abs(step) < 1e-15) break;
  }
  pt.assign(q.begin(), q.begin() + 2);
  Poincare toq(Qp);
  for (size_t i(0); i < f.size(); ++i) {
    const Vec3D p3(xi*Vec3D(f[i]));
    Vec4D v(sqrt(p3.Sqr() + m2[i]), p3);
    toq.BoostBack(v);
    pt.push_back(v);
  }
  pt.push_back(k);
  return true;
}

/*
  YFS: RV_PHOTON_CT (default none) - the emitted photon's charge renormalisation
  in the real-virtual. With REAL_ALPHA0 1 the real couples the photon with
  alpha(0), but the loop provider renormalises it in the input scheme; in a
  G_mu card (OpenLoops, Z pole) v_{n+1} - v -> -0.0075 as the photon goes
  soft. The counterterm c is added to every v_{n+1} (one external photon).
  (name or old integer)
  none (0): no counterterm (right for an alpha(0) card, c = 0 there).
  calibrated (1): once per run, c = -(v_{n+1} - v_B) at points with a photon
     of x = 1e-4 and 2e-4 added to the first event's Born (SoftPhotonPoint),
     extrapolated linearly to x = 0. It forces the soft limit to 0 by
     construction, so a soft-limit test is then no test, and it is taken at
     one point of one process.
  analytic (2): AnalyticPhotonCounterterm, Photon_Counterterm.C: the on-shell
     photon's counterterm 2 Re[dZe(alpha(0)) - dZe(scheme)] alpha/(4 pi) from
     the model's masses and widths, as OpenLoops computes it for its pdg 2002
     photon with Sherpa's settings (ew_scheme 2, EW_REN_SCHEME, CMS):
     0 for alpha(0), 1 - alpha(0)/alpha for alpha(M_Z), Delta r with the
     light-fermion Delta alpha replaced by 1 - alpha(0)/alpha for G_mu
     (Z pole: 0.0075211 against 0.0075209 from 1). Independent of the
     kinematics, so the soft limit stays a test.
  OpenLoops' own on-shell photon (OL_RV_PHOTON_SCHEME On_shell, pdg 2002)
  gives the same constant at every point (pyol: 0.00752 at FSR x = 0.1,
  0.01 and ISR x = 0.01), but once a 2002 process has been evaluated
  OpenLoops multiplies every later tree with a photon by alpha(0)/alpha, the
  Real_Generator's real included (-8% on sigma_fid); calibrated and analytic
  keep the 22 registration.
*/
static rvphotonct::code PhotonCTMode() {
  static const rvphotonct::code m(ATOOLS::Settings::GetMainSettings()["YFS"]
    ["RV_PHOTON_CT"].SetDefault(rvphotonct::none).Get<rvphotonct::code>());
  return m;
}

double NLO_Base::AnalyticPhotonCounterterm() const {
  YFS::Photon_CT_Input in;
  in.MZ = Flavour(kf_Z).Mass();      in.WZ = Flavour(kf_Z).Width();
  in.MW = Flavour(kf_Wplus).Mass();  in.WW = Flavour(kf_Wplus).Width();
  in.MH = Flavour(kf_h0).Mass();     in.WH = Flavour(kf_h0).Width();
  in.MT = Flavour(kf_t).Mass();      in.WT = Flavour(kf_t).Width();
  in.me = Flavour(kf_e).Mass();      in.mmu = Flavour(kf_mu).Mass();
  in.mtau = Flavour(kf_tau).Mass();  in.mb = Flavour(kf_b).Mass();
  in.alpha0 = MODEL::aqed->AqedThomson();
  in.alpha_in = MODEL::s_model->ScalarConstant("alpha_QED");
  in.cms = true;
  for (kf_code q : {kf_d, kf_u, kf_s, kf_c})
    if (Flavour(q).Mass() > 0.)
      msg_Error()<<METHOD<<"(): massive light quark "<<Flavour(q)
                 <<"; the counterterm assumes massless u, d, s, c as in OpenLoops."<<std::endl;
  YFS::photon_ct_scheme sch(YFS::photon_ct_scheme::alpha0);
  switch (ToType<MODEL::ew_scheme::code>(rpa->gen.Variable("EW_REN_SCHEME"))) {
  case MODEL::ew_scheme::alpha0:  sch = YFS::photon_ct_scheme::alpha0;  break;
  case MODEL::ew_scheme::Gmu:     sch = YFS::photon_ct_scheme::Gmu;     break;
  case MODEL::ew_scheme::alphamZ: sch = YFS::photon_ct_scheme::alphamZ; break;
  default:
    msg_Error()<<METHOD<<"(): EW_REN_SCHEME not alpha0, Gmu or alphamZ; "
               <<"no photon counterterm."<<std::endl;
    return 0.;
  }
  std::string log;
  const double c(YFS::OnShellPhotonCounterterm(in, sch, log));
  msg_Info()<<"YFS: RV_PHOTON_CT 2: "<<log<<std::endl;
  return c;
}

void NLO_Base::CalibratePhotonCounterterm() {
  m_rvct_done = true;
  m_rvct = 0.;
  if (PhotonCTMode() == rvphotonct::analytic) {
    m_rvct = AnalyticPhotonCounterterm();
    return;
  }
  if (PhotonCTMode() != rvphotonct::calibrated) return;
  // two soft points, x0 and 2 x0, extrapolated linearly to x = 0
  const double x0(1e-4);
  double dv[2];
  for (int i(0); i < 2; ++i) {
    Vec4D_Vector pt;
    double v1(0.);
    if (!SoftPhotonPoint(m_plab, (i + 1)*x0, pt) || !RealVirtualFactor(pt, v1)) {
      msg_Error()<<METHOD<<"(): calibration point failed, no photon counterterm."
                 <<std::endl;
      return;
    }
    dv[i] = v1 - m_vborn_own;
  }
  m_rvct = -(2.*dv[0] - dv[1]);
  msg_Info()<<"YFS: RV_PHOTON_CT: external-photon counterterm c = "<<m_rvct
            <<" (v_{n+1} - v_B at x = "<<x0<<", "<<2.*x0<<": "<<dv[0]<<", "
            <<dv[1]<<")"<<std::endl;
}

/*
  RV_MODE 1: photon i's non-factorisable real-virtual, B rho_i (v_{n+1} - v),
  on the point and with the ratio CalculateReal() left in m_rvinfo[i].
*/
double NLO_Base::CalculateRealVirtualRemainder(size_t i, const Vec4D &k) {
  if (!m_realvirt || i >= m_rvinfo.size() || !m_rvinfo[i].ok) return 0.;
  if (m_rv_soft_cut > 0. && k.E() < m_rv_soft_cut * sqrt(m_s)) {
    m_softRV++;
    return 0.;
  }
  static const std::string pterm(ATOOLS::Settings::GetMainSettings()["YFS"]
                                 ["RV_PROBE_TERM"].SetDefault("none").Get<std::string>());
  const RVPointInfo &info(m_rvinfo[i]);
  // the other photons' factors of REAL_COMBINE's product: the remainder is
  // photon i's, the others are factorised on it as in the product
  double others(1.);
  for (size_t k(0); k < m_prodfac.size(); ++k) if (k != i) others *= m_prodfac[k];
  if (pterm == "form_factor_only") {
    // diagnostic fast path: no loop call, only the form-factor difference
    const double Keps((m_rrtool && RRMode() == rrmode::exact) ? m_rr_soft_cut*sqrt(m_s) : 0.5*sqrt(m_s));
    Vec4D_Vector legs(info.p.begin(), info.p.end() - 1);
    p_nlodipoles->MakeDipolesII(m_flavs, legs, legs);
    p_nlodipoles->MakeDipoles(m_flavs, legs, legs);
    p_nlodipoles->MakeDipolesIF(m_flavs, legs, legs);
    // the (n+1) loop is never called here, so its provider may not exist
    // yet: the Born loop's IR scale and eps convention for both
    p_nlodipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
    const double a(p_nlodipoles->RealSoftSum(Keps));
    p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
    const double c(m_born*info.rho*(a - p_dipoles->RealSoftSum(Keps))/m_rescale_alpha*others);
    return IsBad(c) ? 0. : c;
  }
  if (!m_vborn_own_ok) { ++m_rvNoVirt; return 0.; }
  double v1(0.);
  if (!RealVirtualFactor(info.p, v1)) {
    if (p_realvirt->FailCut()) m_failcut = true;
    m_zeroRV++;
    if (RVProbeOn()) std::cerr<<"@@@ RVREM failed x="<<2.*k.E()/sqrt(m_s)<<"\n";
    return 0.;
  }
  /*
    One IR convention for the Born virtual and the (n+1) loop (2026-09-29,
    the user's option A). The event carries (1 + v), v = V/B - B(generation
    legs) (m_vborn, the form factor's legs); the exact one-photon O(alpha^2)
    weight (1 + v)(1 + delta) + rho (v_{n+1} - v) needs v_{n+1} with the same
    subtraction on the same legs, so v_{n+1} - v = lt - V/B.
    The pre-2026-09-29 dv was v_{n+1}(own legs) - v_B(own legs): it missed
    rho (v_B - v) (dv_own below; CEEX: REAL_VIRTUAL 2 still reads it).
    Open (notes sec. 11.12): with RR_MODE 1 the double real's soft x hard log
    (RR_SOFT_CUT dependence) is not cancelled by this remainder, nor by the
    coherent real-soft form-factor difference B~(pt) - B~(gen) (bt_pt, bt_gen,
    RV_PROBE_TERM form_factor), which is 10 times the double real's log.
  */
  const double dv_own(v1 + m_rvct - m_vborn_own);
  const double vb_ev(m_vborn + m_virt_subval/m_born);   // V/B of the Born point
  double bt_pt(0.), bt_gen(0.);
  const double Keps((m_rrtool && RRMode() == rrmode::exact) ? m_rr_soft_cut*sqrt(m_s) : 0.5*sqrt(m_s));
  {
    Vec4D_Vector legs(info.p.begin(), info.p.end() - 1);
    p_nlodipoles->MakeDipolesII(m_flavs, legs, legs);
    p_nlodipoles->MakeDipoles(m_flavs, legs, legs);
    p_nlodipoles->MakeDipolesIF(m_flavs, legs, legs);
    p_nlodipoles->p_yfsFormFact->p_virt = p_realvirt->p_loop_me.get();
    bt_pt = p_nlodipoles->RealSoftSum(Keps);
    p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
    bt_gen = p_dipoles->RealSoftSum(Keps);
  }
  const double dv(m_rvLastLt + m_rvct - vb_ev);
  /*
    YFS: RR_LOG_CHECK N (diagnostic): for the first N photons with x > 0.05,
    the soft x hard log the form-factor difference must cancel: d/dlnK of
    B~(pt; K) - B~(gen; K) at K = RR_SOFT_CUT sqrt(s), against the angular
    integral at |k| = K of the soft factors the double real sees, the
    coherent eikonal on this point's legs (exact R_2) minus the product's
    (the crude on the generation legs plus their IF interference), both
    also separately ("@@@ RRLOG"; the three ratios must be one constant).
  */
  { static int nlog(ATOOLS::Settings::GetMainSettings()["YFS"]["RR_LOG_CHECK"]
                    .SetDefault(0).Get<int>());
    if (nlog > 0 && 2.*k.E()/sqrt(m_s) > 0.05) {
      --nlog;
      Vec4D_Vector legs(info.p.begin(), info.p.end() - 1);
      auto bt = [&](double K, bool pt) {
        if (pt) {
          p_nlodipoles->MakeDipolesII(m_flavs, legs, legs);
          p_nlodipoles->MakeDipoles(m_flavs, legs, legs);
          p_nlodipoles->MakeDipolesIF(m_flavs, legs, legs);
          p_nlodipoles->p_yfsFormFact->p_virt = p_realvirt->p_loop_me.get();
          return p_nlodipoles->RealSoftSum(K);
        }
        p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
        return p_dipoles->RealSoftSum(K); };
      const double h(0.5);
      const double dpt((bt(Keps*exp(h), true) - bt(Keps*exp(-h), true))/(2.*h));
      const double dgen((bt(Keps*exp(h), false) - bt(Keps*exp(-h), false))/(2.*h));
      p_nlodipoles->MakeDipolesII(m_flavs, legs, legs);
      p_nlodipoles->MakeDipoles(m_flavs, legs, legs);
      p_nlodipoles->MakeDipolesIF(m_flavs, legs, legs);
      double ipt(0.), igen(0.);
      const int N(20000);
      for (int n(0); n < N; ++n) {
        const double z(1. - (2.*n + 1.)/N), r(sqrt(Max(0., 1. - z*z)));
        const double ph(n*M_PI*(3. - sqrt(5.)));
        const Vec4D q(Keps, Keps*r*cos(ph), Keps*r*sin(ph), Keps*z);
        ipt  += p_nlodipoles->CalculateRealSub(q);
        igen += p_dipoles->CalculateRealSubEEX(q) + p_dipoles->CalculateRealSubIF(q);
      }
      {
        Vec4D_Vector dirs;
        for (int n(0); n < 4000; ++n) {
          const double z(1. - (2.*n + 1.)/4000.), r(sqrt(Max(0., 1. - z*z)));
          const double ph(n*M_PI*(3. - sqrt(5.)));
          dirs.push_back(Vec4D(Keps, Keps*r*cos(ph), Keps*r*sin(ph), Keps*z)*(1./Keps));
        }
        p_nlodipoles->p_yfsFormFact->p_virt = p_realvirt->p_loop_me.get();
        std::cerr<<"@@@ RRLOGPT"<<p_nlodipoles->RealSoftSumReport(Keps, dirs)<<"\n";
        p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
        std::cerr<<"@@@ RRLOGGEN"<<p_dipoles->RealSoftSumReport(Keps, dirs)<<"\n";
      }
      ipt *= 4.*M_PI*Keps*Keps/N;
      igen *= 4.*M_PI*Keps*Keps/N;
      std::ostringstream o;
      o<<std::setprecision(8)<<"@@@ RRLOG x="<<2.*k.E()/sqrt(m_s)<<" fsr="<<(PhotonIsFSR(k)?1:0)
       <<" nph="<<m_photons.size()<<" dBt/dlnK="<<dpt - dgen<<" Idiff="<<ipt - igen
       <<" ratio="<<(dpt - dgen)/(ipt - igen)<<" ratio_pt="<<dpt/ipt<<" ratio_gen="<<dgen/igen
       <<std::setprecision(15)<<" dpt="<<dpt<<" ipt="<<ipt<<" legs=";
      for (const Vec4D &q : legs) o<<q<<" ";
      o<<"\n";
      std::cerr<<o.str();
    } }
  /*
    YFS: RV_PROBE_TERM (diagnostic, default none): the RealVirtual stream
    carries only one piece of the remainder, to measure it on the same
    events: vB_minus_v = rho (v_B - v), the Born-virtual leg term the
    pre-2026-09-29 remainder lacked; form_factor = rho (B~(pt) - B~(gen))/kappa
    (form_factor_only: the same without calling the loop at all);
    old = the pre-2026-09-29 remainder rho (v_{n+1}(own) - v_B).
  */
  double use(dv);
  if (pterm == "vB_minus_v") use = m_vborn_own - m_vborn;
  else if (pterm == "form_factor") use = (bt_pt - bt_gen)/m_rescale_alpha;
  else if (pterm == "old") use = dv_own;
  else if (pterm != "none")
    THROW(fatal_error, "YFS: RV_PROBE_TERM must be none, vB_minus_v, form_factor, "
          "form_factor_only or old");
  const double contrib(m_born*info.rho*use*others);
  if (!IsBad(dv_own)) m_rvdv.push_back(std::make_pair(k, dv_own));
  if (RVProbeOn()) {
    std::ostringstream o;
    // the photon's record, on its own line after the loop call's @@@ RVREM
    o<<std::setprecision(10)<<"@@@ RVPH x="<<2.*k.E()/sqrt(m_s)
     <<" xpt="<<2.*info.p.back().E()/sqrt(m_s)
     <<" fsr="<<(PhotonIsFSR(k)?1:0)<<" nph="<<m_photons.size()
     <<" v="<<m_vborn_own<<" vev="<<m_vborn<<" dv="<<dv<<" dv_own="<<dv_own
     <<" dvA="<<m_rvLastLt + m_rvct - vb_ev<<" dBt="<<(bt_pt - bt_gen)/m_rescale_alpha
     <<" dBt10="<<[&]() {
        Vec4D_Vector legs(info.p.begin(), info.p.end() - 1);
        p_nlodipoles->MakeDipolesII(m_flavs, legs, legs);
        p_nlodipoles->MakeDipoles(m_flavs, legs, legs);
        p_nlodipoles->MakeDipolesIF(m_flavs, legs, legs);
        p_nlodipoles->p_yfsFormFact->p_virt = p_realvirt->p_loop_me.get();
        const double a(p_nlodipoles->RealSoftSum(10.*Keps));
        p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
        return (a - p_dipoles->RealSoftSum(10.*Keps))/m_rescale_alpha; }()
     <<" rho="<<info.rho
     <<" rv/B="<<contrib/m_born<<"\n";
    if (std::abs(dv) > 3e-3 && 2.*info.p.back().E()/sqrt(m_s) < 2e-3) {
      o<<"   p=";
      for (const Vec4D &q : info.p) o<<q<<" ";
      o<<"\n   plab=";
      for (const Vec4D &q : m_plab) o<<q<<" ";
      o<<"\n   photons=";
      for (const YFS::Photon &g : m_photons) o<<g.K()<<(g.IsFSR()?"F ":"I ");
      o<<"\n";
    }
    std::cerr<<o.str();
  }
  if (IsBad(contrib)) { m_zeroRV++; return 0.; }
  m_nonZeroRV++;
  return contrib;
}

double NLO_Base::CalculateRealVirtual() {
  if (!m_realvirt)
    return 0;
  m_rv_hard1 = 0.;
  m_rv_hard2 = 0.;
  if (RVMode() == rvmode::remainder) {
    double rv(0.);
    m_rvdv.clear();
    m_vborn_own_ok = m_vborn_ok && BornVirtualOnOwnLegs(m_vborn_own);
    if (m_vborn_own_ok && !m_rvct_done) CalibratePhotonCounterterm();
    for (size_t i(0); i < m_photons.size(); ++i) {
      YFS::Photon &g(m_photons[i]);
      const double contrib(CalculateRealVirtualRemainder(i, g.K()));
      rv += contrib;
      g.m_beta11 = contrib;
    }
    HardestBetas(m_photons, [](const YFS::Photon &g) { return g.beta11(); },
                 m_rv_hard1, m_rv_hard2);
    if (RVProbeOn() && m_born != 0.) {
      double xmax(0.);
      for (const YFS::Photon &g : m_photons) xmax = Max(xmax, 2.*g.E()/sqrt(m_s));
      std::ostringstream o;
      o<<std::setprecision(10)<<"@@@ RVEV xmax="<<xmax<<" nph="<<m_photons.size()
       <<" v="<<m_vborn<<" rv/B="<<rv/m_born<<"\n";
      std::cerr<<o.str();
    }
    return rv;
  }
  if (m_rv_hard_photon==2) {
    Vec4D k = FixedTestPhoton();
    if (m_check_rv) CheckRealVirtualSub(k);
    m_rv_hard1 = CalculateRealVirtual(k);
    m_rv_hard2 = m_rv_hard1;  // single photon: 2-photon sum == 1-photon
    return m_rv_hard1;
  }
  if (m_rv_hard_photon==1) {
    Vec4D k = MostEnergeticPhoton();
    if (k.E() == 0.) return 0;
    if (m_check_rv) {
      if (k.E() < 0.2 * sqrt(m_s)) return 0;
      CheckRealVirtualSub(k);
    }
    m_rv_hard1 = CalculateRealVirtual(k);
    m_rv_hard2 = m_rv_hard1;  // single photon: 2-photon sum == 1-photon
    return m_rv_hard1;
  }
  double realvirtual(0);
  m_rv_hard2 = 0.;
  for (YFS::Photon &g : m_photons) {
    const Vec4D k(g.K());
    if (m_check_rv) {
      if (k.E() < 0.2 * sqrt(m_s))
        continue;
      CheckRealVirtualSub(k);
    }
    double contrib = CalculateRealVirtual(k);
    realvirtual += contrib;
    g.m_beta11 = contrib;
  }
  HardestBetas(m_photons, [](const YFS::Photon &g) { return g.beta11(); },
               m_rv_hard1, m_rv_hard2);
  return realvirtual;
}

double NLO_Base::CalculateRealVirtual(Vec4D k) {
  if (!m_realvirt) return 0;
  if (m_rv_soft_cut > 0. && k.E() < m_rv_soft_cut * sqrt(m_s)) {
    m_softRV++;
    msg_Debugging() << METHOD << " skipping soft photon: E=" << k.E()
                    << " < RV_SOFT_CUT*sqrt(s)=" << m_rv_soft_cut * sqrt(m_s)
                    << "\n";
    return 0;
  }
  Vec4D_Vector p(m_plab), pi(m_bornMomenta), pf(m_bornMomenta);
  double tot(0), sub(0);
  double norm = 2 * pow(2 * M_PI, 3);
  double flux(1);
  Vec4D kk = k;
  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);
  MapMomenta(p, k);
  double yfspole;
  p.push_back(k);
  CheckMasses(p, 1);
  Vec4D_Vector pp = p;
  pp.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, pp, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, pp, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, pp, m_plab);
  p_nlodipoles->p_yfsFormFact->p_virt = p_realvirt->p_loop_me.get();
  // the dim-reg B~ reads the loop ME's IR scale, which is 0 until the first
  // loop call: without this the first RV of every rank was NaN (and the run's
  // total weight inf); afterwards the scale is the same every call
  if (!(p_realvirt->p_loop_me->IRscale() > 0.)) p_realvirt->Calc(p, m_born);
  double subloc = p_nlodipoles->CalculateRealVirtualSubEps(k);
  yfspole = p_nlodipoles->Get_E1();
  // Eikonal factor S(k) at the reduced kinematics.
  const double eikloc = p_nlodipoles->CalculateRealSub(k);
  const double aB = eikloc * m_oneloop;
  if (m_flux_mode == fluxmode::mapped)
    flux = p_nlodipoles->CalculateFlux(k);
  else if (m_flux_mode == fluxmode::average)
    flux = 0.5 *
           (p_nlodipoles->CalculateFlux(kk) + p_nlodipoles->CalculateFlux(k));
  else
    flux = p_dipoles->CalculateFlux(kk);
  double subb;

  // Both arms of this were the same call; the fsrcount test did nothing.
  subb = p_dipoles->CalculateRealSubEEX(kk);
  if (p.size() != (m_flavs.size() + 1)) {
    msg_Error() << "Mismatch in " << METHOD << std::endl;
  }
  double r = p_realvirt->Calc(p, m_born) / norm * BornPhotonSym(1);
  m_rv = r * flux;
  if (p_realvirt->FailCut())
    m_failcut = true;
  ;
  if (IsBad(r)) {
    m_zeroRV++;
    msg_Error() << "Real-Virtual is " << r << std::endl;
    return 0;
  }
  if (IsZero(r,1e-30)) {
    m_zeroRV++;
    msg_Error() << "Real-Virtual is " << r << std::endl;
    return 0;
  }
  double rtree(0.);
  if (m_realtool)
    rtree = p_real->Calc_R(p) / norm * BornPhotonSym(1);
  else
    msg_Error() << METHOD << ": no real-emission ME available, "
                << "RV YFS subtraction is incomplete.\n";
  m_rvsub = (subloc * rtree * flux + aB) / m_rescale_alpha;
  const double d1 = r * flux - m_rvsub;
  msg_Debugging() << METHOD << " r*flux=" << r * flux
                  << " rtree*flux=" << rtree * flux
                  << " B_fin=" << subloc << " S(k)=" << eikloc
                  << " oneloop=" << m_oneloop << " d1=" << d1 << "\n";
  if (m_submode == submode::local)
    tot = d1 / eikloc;
  else if (m_submode == submode::global)
    tot = d1 / subb;
  else if (m_submode == submode::off)
    tot = (r * flux) / subb;
  if (RVProbeOn()) ProbeRealVirtual(kk, p, r, rtree, subloc, eikloc, flux, subb, tot);
  const double rvmax = ATOOLS::Max(fabs(m_rv), fabs(m_rvsub));
  const double C = (rvmax > 0.) ? fabs(d1) / rvmax : 0.;
  if (m_rv_cancel_hist) {
    if (C > 0.) {
      const double logC = ATOOLS::Max(-16., log10(C));
      m_histograms1d["RV_tot_by_logC_w"]->Insert(logC, tot);
      m_histograms1d["RV_tot_by_logC_n"]->Insert(logC, 1.);
    }
    const double efrac = k.E() / sqrt(m_s);
    m_histograms1d["RV_tot_by_Efrac_w"]->Insert(efrac, tot);
    m_histograms1d["RV_tot_by_Efrac_n"]->Insert(efrac, 1.);
    const double costh_beam = k.CosTheta();
    double maxcos_leg = -1.;
    for (size_t i = 0; i < m_flavs.size(); ++i) {
      if (m_flavs[i].Charge() == 0.) continue;
      maxcos_leg = ATOOLS::Max(maxcos_leg, k.CosTheta(m_plab[i]));
    }
    const bool hardwide = (efrac > 0.1) && (maxcos_leg < 0.9);
    const double subscale =
        ATOOLS::Max(ATOOLS::Max(fabs(m_rvsub), fabs(aB)), 1e-300);
    const double lr = log10(ATOOLS::Max(fabs(m_rv) / subscale, 1e-300));
    m_histograms1d["RV_MEstab_all"]->Insert(lr, 1.);
    m_histograms1d["RV_tot_by_MEstab_w"]->Insert(lr, tot);
    if (hardwide) m_histograms1d["RV_MEstab_hardwide"]->Insert(lr, 1.);
    if (C >= 1.) m_rvHiC++;
    const double bscale = ATOOLS::Max(fabs(m_born), 1e-30);
    if (fabs(tot) > 1.e3 * bscale) {
      m_rvBlowup++;
      if (rtree == 0.) m_rvBlowupRtree0++;
      if (C >= 1.) m_rvBlowupHiC++;
      if (efrac < 0.01) m_rvBlowupSoft++;
      if (hardwide) m_rvBlowupHardWide++;
      int rank = 0;
#ifdef USING__MPI
      if (mpi->Size() > 1) rank = mpi->Rank();
#endif
      std::ofstream bf(std::string(m_debugDIR_NLO) + "/RV_blowups_rank" +
                           std::to_string(rank) + ".txt",
                       std::ios_base::app);
      bf << std::setprecision(10)
         << "E=" << k.E() << " Efrac=" << efrac
         << " costh_beam=" << costh_beam << " maxcos_leg=" << maxcos_leg
         << " pT=" << k.PPerp() << " hardwide=" << (hardwide ? 1 : 0)
         << " tot=" << tot << " d1=" << d1 << " C=" << C
         << " | rv(r*flux)=" << (r * flux) << " rvsub=" << m_rvsub
         << " rtree*flux=" << (rtree * flux) << " subloc(Bfin)=" << subloc
         << " aB(eik*oneloop)=" << aB << " eikloc(S_k)=" << eikloc
         << " subb=" << subb << " flux=" << flux << " oneloop=" << m_oneloop
         << " rescale=" << m_rescale_alpha << " submode=" << (int)m_submode
         << "\n";
    }
  }
  if (m_rv_me_max_ratio > 0.) {
    const double subscale_g =
        ATOOLS::Max(ATOOLS::Max(fabs(m_rvsub), fabs(aB)), 1e-300);
    if (fabs(m_rv) > m_rv_me_max_ratio * subscale_g) {
      m_rvUnstable++;
      msg_Debugging() << METHOD << " skipping RV, unstable loop ME: |rv|="
                      << fabs(m_rv) << " > " << m_rv_me_max_ratio
                      << "*max(|rvsub|,|aB|)=" << (m_rv_me_max_ratio * subscale_g)
                      << " (E=" << k.E() << ", d1=" << d1 << ")\n";
      return 0;
    }
  }
  if (m_rv_cancel_eps > 0. && fabs(d1) < m_rv_cancel_eps * rvmax) {
    m_softRV++;
    msg_Debugging() << METHOD << " skipping RV, cancellation C=" << C
                    << " < RV_CANCEL_EPS=" << m_rv_cancel_eps
                    << " (d1=" << d1 << ", max(|rv|,|rvsub|)=" << rvmax << ")\n";
    return 0;
  }
  if (m_check_poles == 1 && r != 0) {
    double pr1 =
        p_realvirt->p_loop_me->ME_E1() * p_realvirt->m_factor * flux / norm;
    double pr2 = p_realvirt->p_loop_me->ME_E1() * p_realvirt->m_factor;
    // p_nlodipoles->CalculateRealSub(k)*
    const double correctdigit =
        ::countMatchingDigits(pr2, -p_nlodipoles->Get_E1());
    m_histograms1d["RVSinglePoleCD"]->Insert(correctdigit);
    m_histograms1d["RealLoopEpsLP"]->Insert(log10(fabs(pr1)));
    m_histograms1d["RealLoopEpsYFS"]->Insert(log10(fabs(yfspole)));
    if (!IsEqual(pr2, -yfspole, 1e-4)) {
      msg_Out() << "Poles do not cancel in YFS Real-Virtuals" << std::endl
                << "Process =  " << p_realvirt->p_loop_me->Name() << std::endl
                << "Correct Digits =  " << correctdigit << std::endl
                << "One-Loop Provider RV eps^{-1}  = " << pr2 << std::endl
                << "Sherpa RV eps^{-1} = " << -yfspole << std::endl
                << "Sherpa/One-Loop = " << -yfspole / pr2 << std::endl;
      return 0;
    } else {
      msg_Debugging() << std::setprecision(16)
                      << "Poles cancel in YFS Real-Virtuals" << std::endl
                      << "Process =  " << p_realvirt->p_loop_me->Name()
                      << std::endl
                      << "Correct Digits =  " << correctdigit << std::endl
                      << "One-Loop Provider RV eps^{-1}  = " << pr2 << std::endl
                      << "Sherpa RV eps^{-1} = " << p_nlodipoles->Get_E1()
                      << std::endl;
    }
  }
  if (IsZero(tot)) m_zeroRV++;
  else m_nonZeroRV++;
  return tot;
}
