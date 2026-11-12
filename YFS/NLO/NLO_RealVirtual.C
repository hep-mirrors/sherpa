#include <chrono>
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
  YFS: RV_PROBE (diagnostic, default 0): one "@@@ RVPROBE" line per
  real-virtual photon on std::cerr with the loop over tree, the YFS dim-reg
  subtraction on the Born legs and on the (n+1) point's legs, and the pole
  coefficients. dv_ev and dv_pt = v_{n+1} - v must vanish in the soft limit.
*/
static bool RVProbeOn() {
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PROBE"]
                       .SetDefault(0).Get<int>() != 0);
  return on;
}

bool NLO_Base::RVRemainder() const { return m_realvirt && RVIsRemainder(); }

/*
  RV_MODE 1: v_{n+1} at the (n+1)-body point p (photon last), in units of
  v = V_sub/B: the provider's loop over its tree minus the Born virtual's
  subtraction on the point's charged legs. The subtraction is built after
  the loop call so a dim-reg B~ reads this call's IR scale. The 1/eps
  coefficients of loop and subtraction must cancel (counted in m_rvPoleFail).
*/
/*
  YFS: RV_LOOP_FRAME, the frame the (n+1)-point loop is handed to OpenLoops
  in (name or old integer).
  lab (0): the real's point as is.
  canonical (1): CanonicalBeamFrame, the beams' rest frame, generically rotated.
  tilted (2, paper setting): beams exactly on z in their rest frame, tilted
     by 1e-3 rad.
  beam_axis (3): beams exactly on z in their rest frame, no tilt.
  No single frame is stable for photons deep in the electron's dead cone;
  RV_LOOP_CHECK compares two (see RealVirtualFactor).
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
  YFS: RV_LOOP_CHECK. For a photon within RV_LOOP_CHECK_ANGLE (default 1e-4
  rad) of a beam, evaluate V_fin/T in a second frame (tilted if RV_LOOP_FRAME
  is beam_axis, else beam_axis) and drop the photon's real-virtual when the
  two differ by more than RV_LOOP_CHECK_TOL (default 1e-3).
*/
/*
  YFS: RV_LOOP_CHECK_FINAL (default true). Also check photons within
  RV_LOOP_CHECK_ANGLE of an outgoing charged leg (flavours given, p = legs +
  photon), where the loop can be unstable for collinear final-state photons.
  false restores the beam-only check.
*/
static bool RVLoopCheckNeeded(const Vec4D_Vector &p, const Flavour_Vector *fl = NULL)
{
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_LOOP_CHECK"]
                       .SetDefault(true).Get<bool>());
  static const double amax(ATOOLS::Settings::GetMainSettings()["YFS"]
                           ["RV_LOOP_CHECK_ANGLE"].SetDefault(1e-4).Get<double>());
  static const bool fin(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_LOOP_CHECK_FINAL"]
                        .SetDefault(true).Get<bool>());
  if (!on || p.size() < 3) return false;
  const Vec4D &k(p.back());
  size_t n(2);
  if (fin && fl && fl->size() + 1 == p.size()) n = fl->size();
  for (size_t i(0); i < n; ++i) {
    if (i >= 2 && (*fl)[i].Charge() == 0.) continue;
    const double c(Vec3D(p[i])*Vec3D(k)/(Vec3D(p[i]).Abs()*Vec3D(k).Abs()));
    if (acos(Max(-1., Min(1., c))) < amax) return true;
  }
  return false;
}

/*
  YFS: RV_BEAM_CONE (radians; default auto). A real-virtual photon within
  this angle of a beam gets no remainder (dv = 0, counted). auto (-1) is 1e-3
  for the Recola provider, whose loop is unreliable there, and 0 otherwise.
*/
bool NLO_Base::RVInBeamCone(const Vec4D_Vector &p) const
{
  static const double set(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_BEAM_CONE"]
                          .SetDefault(-1.).Get<double>());
  const double cone(set >= 0. ? set
                    : (p_realvirt != NULL && p_realvirt->m_gen == "Recola" ? 1e-3 : 0.));
  if (!(cone > 0.) || p.size() < 3) return false;
  const Vec3D k(p.back());
  for (size_t i(0); i < 2; ++i) {
    const Vec3D b(p[i]);
    const double c(b*k/(b.Abs()*k.Abs()));
    if (acos(Max(-1., Min(1., c))) < cone) return true;
  }
  return false;
}

// ACRAIC: HYBRID_ISR_RMIN r_c: V_2 if min(r_1, r_2) < r_c (always for r_c < 0)
bool NLO_Base::InHybridIsrCone(const Vec4D &k) const
{
  static const double rc(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["HYBRID_ISR_RMIN"]
                         .SetDefault(-1.).Get<double>());
  if (rc < 0.) return true;
  return Min(k[0] - k[3], k[0] + k[3])/sqrt(m_s) < rc;
}

bool NLO_Base::RealVirtualFactor(const Vec4D_Vector &pin, double &v1) {
  if (pin.size() != m_flavs.size() + 1) return false;
  const Vec4D_Vector p(RVLoopPoint(pin, RVLoopFrame()));
  double mur2(sqr(rpa->gen.Ecms()));
  if (m_looptool && p_virt->p_loop_me->IRscale() > 0.)
    mur2 = sqr(p_virt->p_loop_me->IRscale());
  double lt(0.);
  if (RVLoopCheckNeeded(pin, &m_flavs)) {
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
  m_rvLastSub = sub/m_rescale_alpha;
  m_rvLastLoop = p;
  const double e1l(p_realvirt->E1()*p_realvirt->m_factor);
  const double e1s(p_nlodipoles->Get_E1()/m_rescale_alpha);
  // with massification the poles agree only up to O(m^2/s_ij)
  // YFS: RV_POLE_TOL: tolerance of the (n+1) loop's 1/eps against the subtraction's
  static const double poletolset(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_POLE_TOL"]
                                 .SetDefault(0.).Get<double>());
  const double poletol(poletolset > 0. ? poletolset : (p_realvirt->p_massif ? 1e-4 : 1e-6));
  // diagnostic: one line per loop call with the pole mismatch and the photon
  static const bool polelog(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_POLE_LOG"]
                            .SetDefault(false).Get<bool>());
  if (polelog && m_dim_reg && e1l != 0.) {
    const Vec4D &kp(p.back());
    std::cerr<<std::setprecision(8)<<"@@@ RVPOLE x "<<2.*kp[0]/sqrt((p[0] + p[1]).Abs2())
             <<" rel "<<std::abs(e1l + e1s)/std::max(std::abs(e1l), 1e-3)
             <<" v1 "<<lt - sub/m_rescale_alpha<<"\n";
  }
  const bool polefail(m_dim_reg && std::abs(e1l + e1s) > poletol*std::max(std::abs(e1l), 1e-3));
  if (polefail) ++m_rvPoleFail;
  // YFS: RV_POLE_VETO: dv = 0 when the provider's 1/eps_IR fails that check
  static const bool poleveto(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_POLE_VETO"]
                             .SetDefault(false).Get<bool>());
  if (poleveto && polefail && e1l != 0.) { ++m_rvPoleVeto; return false; }
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
    o<<"\n";
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
  YFS: RV_PHOTON_CT (default analytic): the emitted photon's charge
  renormalisation in the real-virtual. With REAL_ALPHA0 1 the real couples
  the photon with alpha(0) but the loop provider renormalises it in the input
  scheme (G_mu), so v_{n+1} - v does not vanish in the soft limit. The
  counterterm c is added to every v_{n+1}.
  none (0): c = 0 (right for an alpha(0) card).
  calibrated (1): c = -(v_{n+1} - v_B) at two soft points (x = 1e-4, 2e-4)
     added to the first event's Born (SoftPhotonPoint), extrapolated to
     x = 0. Forces the soft limit to 0, so it is not a test of it.
  analytic (2): AnalyticPhotonCounterterm (Photon_Counterterm.C), the
     on-shell photon's 2 Re[dZe(alpha(0)) - dZe(scheme)] alpha/(4 pi) from the
     model's masses and widths. Independent of the kinematics, so the soft
     limit stays a test.
*/
static rvphotonct::code PhotonCTMode() {
  static const rvphotonct::code m(ATOOLS::Settings::GetMainSettings()["YFS"]
    ["RV_PHOTON_CT"].SetDefault(rvphotonct::analytic).Get<rvphotonct::code>());
  return m;
}

double NLO_Base::AnalyticPhotonCounterterm(double &cim) const {
  cim = 0.;
  // the provider's own external-photon convention (OpenLoops or Recola)
  const bool recola(p_realvirt != NULL && p_realvirt->m_gen == "Recola");
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
    if (Flavour(q).Mass() > 0. && !recola)
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
  double c(0.);
  if (recola) {
    // Recola sees the light quarks with the masses Sherpa passes and
    // regulates a massless one at its IR scale
    in.provider = YFS::photon_ct_provider::recola;
    in.mu = Flavour(kf_u).Mass(); in.md = Flavour(kf_d).Mass();
    in.ms = Flavour(kf_s).Mass(); in.mc = Flavour(kf_c).Mass();
    in.mu_ir = ATOOLS::Settings::GetMainSettings()["RECOLA_IR_SCALE"]
      .SetDefault(100).Get<double>();
    c = YFS::RecolaPhotonCounterterm(in, sch, log);
  }
  else c = YFS::OnShellPhotonCounterterm(in, sch, log);
  cim = YFS::PhotonCountertermImaginary(in, sch);
  msg_Info()<<"YFS: RV_PHOTON_CT 2 ("<<(recola ? "Recola" : "OpenLoops")<<"): "
            <<log<<"; amplitude-level 2 Im c/2 = "<<cim<<std::endl;
  return c;
}

void NLO_Base::CalibratePhotonCounterterm() {
  m_rvct_done = true;
  m_rvct = 0.;
  m_rvct_im = 0.;
  if (PhotonCTMode() == rvphotonct::analytic) {
    m_rvct = AnalyticPhotonCounterterm(m_rvct_im);
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
/*
  B~(K) and the final-leg form factor (the exponent's final-state and
  initial-final part) on the (n+1)-body point's own legs (pt) and on the
  generation legs (gen), for the terms of the real-virtual remainder.
*/
void NLO_Base::SoftSumsOnBothLegs(const Vec4D_Vector &legs, const double K,
                                  double &bt_pt, double &bt_gen,
                                  double &dY_pt, double &dY_gen)
{
  p_nlodipoles->MakeDipolesII(m_flavs, legs, legs);
  p_nlodipoles->MakeDipoles(m_flavs, legs, legs);
  p_nlodipoles->MakeDipolesIF(m_flavs, legs, legs);
  p_nlodipoles->p_yfsFormFact->p_virt = p_realvirt->p_loop_me.get();
  bt_pt = p_nlodipoles->RealSoftSum(K);
  // remainder_own_legs and RV_PROBE: the final-state and initial-final form
  // factor on the point's legs minus on the generation legs
  const bool ownlegs(RVMode() == rvmode::remainder_own_legs);
  if (ownlegs || RVProbeOn()) dY_pt = p_nlodipoles->FormFactorSumFinalLegs();
  p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
  bt_gen = p_dipoles->RealSoftSum(K);
  if (ownlegs || RVProbeOn()) dY_gen = p_dipoles->FormFactorSumFinalLegs();
}

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
    const double Keps((m_rrtool && RRFromSingles() && m_rr_soft_cut > 0.)
                      ? m_rr_soft_cut*sqrt(m_s) : 0.5*sqrt(m_s));
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
  // hybrid: no (n+1)-point loop for an initial-state photon that gets V_2
  static const bool fsronly(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["REAL_VIRTUAL"]
                            .SetDefault(acraicrv::automatic).Get<acraicrv::code>()
                            == acraicrv::hybrid);
  if (fsronly && i < m_photons.size() && !m_photons[i].IsFSR() && InHybridIsrCone(k)) {
    ++m_rvIsrSkipped;
    return 0.;
  }
  if (RVInBeamCone(info.p)) { ++m_rvBeamCone; return 0.; }
  double v1(0.);
  if (!RealVirtualFactor(info.p, v1)) {
    if (p_realvirt->FailCut()) m_failcut = true;
    m_zeroRV++;
    if (RVProbeOn()) std::cerr<<"@@@ RVREM failed x="<<2.*k.E()/sqrt(m_s)<<"\n";
    return 0.;
  }
  /*
    One IR convention for the Born virtual and the (n+1) loop: the event
    carries (1 + v), v = V/B - B(generation legs), and the one-photon
    O(alpha^2) weight (1 + v)(1 + delta) + rho (v_{n+1} - v) needs v_{n+1}
    with the same subtraction on the same legs, so v_{n+1} - v = lt - V/B.
  */
  const double dv_own(v1 + m_rvct - m_vborn_own);
  const double vb_ev(m_vborn + m_virt_subval/m_born);   // V/B of the Born point
  double bt_pt(0.), bt_gen(0.), dY_pt(0.), dY_gen(0.);
  const double Keps((m_rrtool && RRFromSingles() && m_rr_soft_cut > 0.)
                     ? m_rr_soft_cut*sqrt(m_s) : 0.5*sqrt(m_s));
  SoftSumsOnBothLegs(Vec4D_Vector(info.p.begin(), info.p.end() - 1), Keps,
                     bt_pt, bt_gen, dY_pt, dY_gen);
  const double dv(m_rvLastLt + m_rvct - vb_ev);
  /*
    YFS: RV_PROBE_TERM (diagnostic, default none): the RealVirtual stream
    carries only one piece of the remainder: vB_minus_v = rho (v_B - v),
    form_factor = rho (B~(pt) - B~(gen))/kappa (form_factor_only: without
    the loop call), old = rho (v_{n+1}(own) - v_B).
  */
  /*
    remainder_own_legs: v_{n+1} - v with the (n+1) loop subtracted on the
    point's own legs plus Y_FF+IF(point) - Y_FF+IF(generation) (see RVMode in
    NLO_Base_Internal.H).
  */
  const double dv_ownlegs(dv_own + (m_vborn_own - m_vborn) + (dY_pt - dY_gen));
  double use(RVMode() == rvmode::remainder_own_legs ? dv_ownlegs : dv);
  if (pterm == "vB_minus_v") use = m_vborn_own - m_vborn;
  else if (pterm == "form_factor") use = (bt_pt - bt_gen)/m_rescale_alpha;
  else if (pterm == "old") use = dv_own;
  else if (pterm != "none")
    THROW(fatal_error, "YFS: RV_PROBE_TERM must be none, vB_minus_v, form_factor, "
          "form_factor_only or old");
  const double contrib(m_born*info.rho*use*others);
  if (!IsBad(dv_own)) {
    m_rvdv.push_back(std::make_pair(k, dv_own));
    m_rvdv_spt.push_back((info.p[0] + info.p[1]).Abs2());
  }
  // ACRAIC: REAL_VIRTUAL helicity - the same loop call resolved by helicity,
  // on the same remainder convention as dv_own (own legs on both points)
  if (!IsBad(dv_own) && RealVirtualHelicityRequested()) {
    RVHelicityFactors f;
    BuildRealVirtualHelicityFactors(k, dv_own, info.p, f);
    m_rvhel.push_back(f);
  }
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
     <<" rho="<<info.rho<<" dY="<<dY_pt - dY_gen
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
  if (RVIsRemainder()) {
    double rv(0.);
    m_rvdv.clear();
    m_rvdv_spt.clear();
    m_rvhel.clear();
    m_loopnhel.clear();
    // ACRAIC: RV_FRACTION f: exact loops on a hashed fraction f of the events
    const double rvfrac(RealVirtualFraction());
    m_rvSkipped = false;
    if (rvfrac < 1.) {
      ATOOLS::Vec4D_Vector key(m_plab);
      for (const YFS::Photon &g : m_photons) key.push_back(g.K());
      if (KinematicsHash(key) >= rvfrac) {
        m_rvSkipped = true;
        ++m_rvFracSkipped;
        return 0.;
      }
      ++m_rvFracEvaluated;
    }
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
    BuildLoopNPairs();   // ACRAIC: LOOP_BETA_MAX >= 2 (no effect otherwise)
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
  // YFS: RV_COLLINEAR_CUT: dv = 0 for a photon collinear to a charged leg
  if (m_rv_collinear_cut > 0.) {
    double pkmin(1e300);
    for (size_t i(0); i < m_plab.size() && i < m_flavs.size(); ++i)
      if (m_flavs[i].Charge() != 0.) pkmin = Min(pkmin, (m_plab[i]*k)/m_s);
    if (pkmin < m_rv_collinear_cut) { m_softRV++; return 0; }
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
  // the dim-reg B~ reads the loop ME's IR scale, which is 0 before the first loop call
  if (!(p_realvirt->p_loop_me->IRscale() > 0.)) p_realvirt->Calc(p, m_born);
  double subloc = p_nlodipoles->CalculateRealVirtualSubEps(k);
  // in the loop's coupling, as the remainder subtracts it (USE_MODEL_ALPHA 0)
  yfspole = p_nlodipoles->Get_E1() / m_rescale_alpha;
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
        ::countMatchingDigits(pr2, -yfspole);
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


/*
  ACRAIC: LOOP_BETA_MAX 2. The two-photon one-loop record of every photon
  pair with both x = 2E/sqrt(s) above LOOP_BETA_XCUT (and BETA2_XCUT,
  LOOP_SOFT_FLOOR). Needs REAL_VIRTUAL helicity. The pair point gets the
  one-photon guards per photon: RV_BEAM_CONE and the RV_LOOP_CHECK two-frame
  comparison.
*/
void NLO_Base::BuildLoopNPairs()
{
  if (RequestedMaxLoopPhotons() < 2 || !RealVirtualHelicityRequested()
      || RealVirtualProvider(2) == NULL) return;
  static double xcut(-1.), sfloor(-1.);
  if (xcut < 0.) {
    xcut = ATOOLS::Settings::GetMainSettings()["ACRAIC"]["LOOP_BETA_XCUT"].SetDefault(-1.).Get<double>();
    // ACRAIC uses a pair only above BETA2_XCUT (its tree beta_2)
    xcut = Max(xcut, ATOOLS::Settings::GetMainSettings()["ACRAIC"]["BETA2_XCUT"].SetDefault(0.05).Get<double>());
    sfloor = ATOOLS::Settings::GetMainSettings()["YFS"]["LOOP_SOFT_FLOOR"].SetDefault(1e-4).Get<double>();
  }
  const double rs(sqrt(m_s)), xmin(Max(xcut, sfloor));
  for (size_t i(0); i < m_photons.size(); ++i) {
    if (2.*m_photons[i].E()/rs <= xmin) continue;
    for (size_t j(i + 1); j < m_photons.size(); ++j) {
      if (2.*m_photons[j].E()/rs <= xmin) continue;
      LoopNHelicityFactors f;
      if (BuildLoopNHelicityFactors(i, j, f)) m_loopnhel.push_back(f);
      else ++m_lnh_fail;
    }
  }
  /*
    ACRAIC: LOOP_LADDER_PROBE with LOOP_BETA_MAX >= 3: the three-photon
    record for events with exactly three photons above the cut (diagnostic)
  */
  static const int ladder(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["LOOP_LADDER_PROBE"]
                          .SetDefault(0).Get<int>());
  if (ladder > 0 && RequestedMaxLoopPhotons() >= 3 && RealVirtualProvider(3) != NULL
      && m_photons.size() == 3) {
    bool all(true);
    for (size_t g(0); g < 3; ++g) all = all && 2.*m_photons[g].E()/rs > xmin;
    if (all) {
      LoopNHelicityFactors f;
      if (BuildLoopNSetFactors(std::vector<size_t>{0, 1, 2}, f)) m_loopnhel.push_back(f);
      else ++m_lnh_fail;
    }
  }
  // restore the one-photon loop provider on the form factor
  if (p_realvirt) p_nlodipoles->p_yfsFormFact->p_virt = p_realvirt->p_loop_me.get();
}


bool NLO_Base::BuildLoopNHelicityFactors(size_t i, size_t j, LoopNHelicityFactors &out)
{
  return BuildLoopNSetFactors(std::vector<size_t>{i, j}, out);
}

bool NLO_Base::BuildLoopNSetFactors(const std::vector<size_t> &S, LoopNHelicityFactors &out)
{
  out = LoopNHelicityFactors();
  const size_t np(S.size());
  YFS::RealVirtual *rv(RealVirtualProvider(np));
  if (np < 2 || rv == NULL || !m_lhelok || !m_vborn_own_ok) return false;
  for (size_t g : S) out.ks.push_back(m_photons[g].K());
  out.k1 = out.ks[0]; out.k2 = out.ks[1];
  const Vec4D ki(out.k1), kj(out.k2);
  // the set's (n+|S|)-point in the lab (EventSetPoint); no record if it fails
  Vec4D_Vector P;
  if (!EventSetPoint(S, P) || !CheckMomentumConservation(P)) return false;
  // the one-photon guards, per photon (P: lab, beams first)
  std::vector<Vec4D_Vector> pk;
  for (const Vec4D &k : out.ks) pk.push_back(Vec4D_Vector{P[0], P[1], k});
  for (const Vec4D_Vector &x : pk) if (RVInBeamCone(x)) { ++m_lnh_cone; return false; }
  const Vec4D_Vector q(RVLoopPoint(P, RVLoopFrame()));
  double mur2(sqr(rpa->gen.Ecms()));
  if (m_looptool && p_virt && p_virt->p_loop_me->IRscale() > 0.)
    mur2 = sqr(p_virt->p_loop_me->IRscale());
  double lt(0.);
  bool needcheck(false);
  for (const Vec4D_Vector &x : pk) needcheck = needcheck || RVLoopCheckNeeded(x);
  if (needcheck) {
    static const double tol(ATOOLS::Settings::GetMainSettings()["YFS"]
                            ["RV_LOOP_CHECK_TOL"].SetDefault(1e-3).Get<double>());
    double lt2(0.);
    const Vec4D_Vector q2(RVLoopPoint(P, RVLoopFrame() == rvloopframe::beam_axis
                                         ? rvloopframe::tilted : rvloopframe::beam_axis));
    if (!rv->LoopOverTree(P, q2, mur2, lt2)) return false;
    ++m_lnh_checked;
    // HelicityAmplitudes reads the provider's last point
    if (!rv->LoopOverTree(P, q, mur2, lt)) return false;
    if (!(std::abs(lt - lt2) <= tol)) { ++m_lnh_unstable; return false; }
  }
  else if (!rv->LoopOverTree(P, q, mur2, lt)) return false;
  // the YFS virtual subtraction on the point's own charged legs
  Vec4D_Vector legs(q.begin(), q.end() - np);
  p_nlodipoles->MakeDipolesII(m_flavs, legs, legs);
  p_nlodipoles->MakeDipoles(m_flavs, legs, legs);
  p_nlodipoles->MakeDipolesIF(m_flavs, legs, legs);
  p_nlodipoles->p_yfsFormFact->p_virt = rv->p_loop_me.get();
  const double sub(p_nlodipoles->CalculateVirtualSub()/m_rescale_alpha);
  if (IsBad(sub)) return false;
  const double e1l(rv->E1()*rv->m_factor), e1s(p_nlodipoles->Get_E1()/m_rescale_alpha);
  const double poletol(rv->p_massif ? 1e-4 : 1e-6);
  if (m_dim_reg && std::abs(e1l + e1s) > poletol*std::max(std::abs(e1l), 1e-3))
    ++m_lnh_polefail;
  PHASIC::Virtual_ME2_Base *lme(rv->p_loop_me.get());
  const Flavour_Vector &flavs(rv->m_flavs);
  if (flavs.size() != m_flavs.size() + np) return false;
  std::vector<Complex> a0, a1;
  std::vector<std::vector<int> > hel;
  if (!lme->HelicityAmplitudes(a0, a1, hel) || a0.empty() || a0.size() != a1.size())
    return false;
  const double alpha(rv->m_factor*2.*M_PI);
  const double mu(lme->IRscale()), mu2(mu*mu);
  Complex K;
  double imsh(0.);
  if (!LoopHelicityIR(flavs, q, alpha, mu2, K, imsh)) return false;
  static const double minborn(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["VIRT_HEL_MINBORN"].SetDefault(1e-3).Get<double>());
  const double mshift(lt - rv->m_factor*lme->ME_Finite());
  std::map<std::vector<int>, size_t> bornidx;
  for (size_t c(0); c < m_lhel.size(); ++c) bornidx[m_lhel[c]] = c;
  // the photon counterterm per external photon: |S| of them here
  const double hnp(0.5*np);
  const Complex commonB(-0.5*m_lrun - 0.5*m_vborn_own_sub, -m_limsh);
  const Complex commonR(-0.5*sub + hnp*m_rvct + 0.5*mshift, -imsh + hnp*m_rvct_im);
  const double v2(lt - sub + np*m_rvct);
  double s0(0.), amax(0.), sd(0.);
  for (size_t c(0); c < a0.size(); ++c) { s0 += std::norm(a0[c]); amax = Max(amax, std::abs(a0[c])); }
  if (!(s0 > 0.)) return false;
  std::vector<Complex> rem(a0.size());
  for (size_t c(0); c < a0.size(); ++c) {
    Complex d2;
    if (std::abs(a0[c]) > minborn*amax) d2 = a1[c]/a0[c] + commonR;
    else { d2 = Complex(0.5*v2, 0.); ++m_lnh_nfall; }
    Complex d0(0.5*m_vborn_own, 0.);
    const std::vector<int> hb(hel[c].begin(), hel[c].end() - np);
    const std::map<std::vector<int>, size_t>::const_iterator it(bornidx.find(hb));
    if (it != bornidx.end() && m_lratio_ok[it->second]) d0 = m_lratio[it->second] + commonB;
    rem[c] = d2 - d0;
    if (IsBad(rem[c].real()) || IsBad(rem[c].imag())) return false;
    sd += std::norm(a0[c])*2.*d2.real();
  }
  // the (n+2)-point identity: <2 Re delta_S>_tree = v_{n+2} + 2c
  const double res(sd/s0 - v2);
  m_lnh_res_sum += res;
  m_lnh_res_max = Max(m_lnh_res_max, std::abs(res));
  ++m_lnh_n;
  out.v2 = v2;
  out.hel.swap(hel); out.a0.swap(a0); out.rem.swap(rem);
  out.ok = true;
  return true;
}
