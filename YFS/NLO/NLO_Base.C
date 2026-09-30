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

using namespace YFS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

/*
  YFS: RV_MODE - how the squared-level real-virtual (NLO_Part with E) is built.

  0 (default): the legacy RV, CalculateRealVirtual(Vec4D): per photon
    [r_RV flux - (B~ rtree flux + S~ B v)/kappa]/S~_crude, the whole
    beta_1^(2) - beta_1^(1) = R v_{n+1} - S~ B v, with its own point, crude,
    flux and coupling, and with E present REAL_COMBINE's product and
    VIRTUAL_COMBINE switched off. IR unsafe: the (n+1) loop runs at
    mu = 1 GeV against the Born's mu^2 = s, and the full-EW finite part is
    mu dependent beyond its IR pole (+0.4% G_mu, +0.8% alpha(0) per event
    whose photons are all soft).

  1: the remainder. With v = V_sub/B of the event, v_{n+1} the same object at
    photon j's (n+1)-body point and rho_j = R_j/(crude_j B) the photon's
    exact-over-crude ratio, the exact one-photon O(alpha^2) weight
      beta_0 (1 + v) + [R (1 + v_{n+1}) - S B (1 + v)]/crude
    is identically
      (1 + v)(1 + delta) + rho (v_{n+1} - v),
    delta = (R - S B)/(crude B) the real's bracket.
*/
/*
  YFS: RR_MODE - how the squared-level double real (NLO_Part with W) is built.
  0 (default): the legacy assembly (CalculateRealReal(k1,k2), steered by
    RR_CONVENTIONS), which also switches REAL_COMBINE's product off.
  1: exact at O(alpha^2) for every photon pair. f_k is photon k's factor in
    REAL_COMBINE's product (its bracket plus its IF ratio, = kappa r F/(c B),
    the exact over crude at its own single-real point), and
        rho_ij = kappa^2 R_2(P) F_ij / (c_i c_j B)
    the pair's exact over crude at its point P (EventPairPoint: the event,
    the other photons reduced), c the singles' densities carried to P by the
    ratio of crude eikonals, F_ij the pair's flux (P_beams - K_ISR)^2/P_beams^2.
    The pair adds
        RR_ij = B (rho_ij - f_i f_j) prod_{k != i,j} f_k,
    so that with two photons product + RR = rho_ij, the exact two-photon ME
    over the generator's density, event by event. RR does NOT vanish when
    one photon goes soft next to a hard one: the exact R_2 attaches the soft
    photon to the legs after the hard emission, the product (and the form
    factor) to the generation legs. That soft x hard log is cancelled by the
    real-virtual of the hard photon (RV_MODE 1, CalculateRealVirtualRemainder,
    the form-factor difference at the RR's cut), so RR_MODE 1 needs RV_MODE 1
    and RR_SOFT_CUT > 0, and the test is sigma independent of RR_SOFT_CUT and
    IR_CUTOFF, not an event-level soft limit (notes sec. 11.12).
    Settings: RR_Generator OpenLoops (Comix's R_2 is unreliable for soft
    photons), RR_SOFT_CUT 1e-4 (OpenLoops' R_2 for two ultra-soft collinear
    initial-state photons, x ~ 1e-7, depends on the frame). Born photons'
    channel weights (REAL_BORN_PHOTON_MULTICHANNEL) are not applied to R_2
    (e+e- -> gamma gamma needs them). RR_PROBE 1 prints "@@@ RRPAIR" per pair
    (2 adds P and the Born point) and "@@@ RRID" per two-photon event.
*/
static int RRMode() {
  static const int m(ATOOLS::Settings::GetMainSettings()["YFS"]["RR_MODE"]
                     .SetDefault(0).Get<int>());
  return m;
}

static int RVMode() {
  static const int m(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_MODE"]
                     .SetDefault(0).Get<int>());
  return m;
}

/*
  YFS: RV_PROBE (diagnostic, default 0): one line per real-virtual photon on
  std::cerr, "@@@ RVPROBE", with the loop over the tree lt = V_fin/T of the
  last loop call, the YFS dim-reg subtraction on the event's Born legs (the
  dipoles CalculateRealVirtual builds, Born momenta m_plab) and on the
  (n+1)-body point's own charged legs, the Born's v = V_sub/B, and the
  one-loop pole coefficients of each. dv_ev and dv_pt are v_{n+1} - v with
  the two subtractions: the object that must vanish in the soft limit.
*/
// YFS: RV_LOOP_FRAME, documented at RVLoopPoint
static int RVLoopFrame() {
  static const int m(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_LOOP_FRAME"]
                     .SetDefault(1).Get<int>());
  return m;
}

static bool RVProbeOn() {
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PROBE"]
                       .SetDefault(0).Get<int>() != 0);
  return on;
}

bool NLO_Base::RVRemainder() const { return m_realvirt && RVMode() == 1; }

/*
  The frame the one-loop provider is called in (RV_MODE 1 always, the Born's
  virtual with YFS: VIRTUAL_CANONICAL_FRAME). OpenLoops' one-loop amplitudes
  are not Lorentz invariant to better than ~1% when a leg is almost, but not
  exactly, on a coordinate axis - measured with pyol at Sherpa's own points

  So: boost to the rest frame of the two beams, put the beams exactly back to
  back on the z axis, then turn the whole point by
  a fixed generic rotation, so that no leg lies near a coordinate axis. The
  squared amplitudes are invariant, so nothing else changes.
*/
static Vec4D_Vector BeamRestFrameOnZ(const Vec4D_Vector &p)
{
  Vec4D_Vector q(p);
  if (q.size() < 2) return q;
  Poincare cms(q[0] + q[1]);
  for (Vec4D &v : q) cms.Boost(v);
  if (sqr(q[0][1]) + sqr(q[0][2]) > 0. || q[0][3] < 0.) {
    const Poincare rot(q[0], Vec4D(0., 0., 0., 1.));
    Vec4D a(q[0]), b(q[0]);
    rot.Rotate(a);
    rot.RotateBack(b);
    const bool fwd(sqr(a[1]) + sqr(a[2]) <= sqr(b[1]) + sqr(b[2]));
    for (Vec4D &v : q) { if (fwd) rot.Rotate(v); else rot.RotateBack(v); }
  }
  q[0] = Vec4D(q[0][0], 0., 0.,  Vec3D(q[0]).Abs());
  q[1] = Vec4D(q[1][0], 0., 0., -Vec3D(q[1]).Abs());
  return q;
}

static Vec4D_Vector CanonicalBeamFrame(const Vec4D_Vector &p)
{
  Vec4D_Vector q(BeamRestFrameOnZ(p));
  // fixed generic rotation (z-y-z Euler angles 0.3, 0.7, 1.1)
  static const double ca(cos(0.3)), sa(sin(0.3)), cb(cos(0.7)), sb(sin(0.7)),
                      cc(cos(1.1)), sc(sin(1.1));
  for (Vec4D &v : q) {
    double x(v[1]), y(v[2]), z(v[3]);
    double x1(ca*x - sa*y), y1(sa*x + ca*y);
    double x2(cb*x1 + sb*z), z2(-sb*x1 + cb*z);
    v = Vec4D(v[0], cc*x2 - sc*y1, sc*x2 + cc*y1, z2);
  }
  return q;
}

double massmin = 2220;
double rcount = 1;
double sumw = 0;


// Lambda (Kaellen function) now lives once in YFS/Tools/Dipole.H.

template <typename Get>
static void HardestBetas(const YFS::Photon_Vector &ph, Get get,
                         double &h1, double &h2) {
  std::vector<const YFS::Photon *> s;
  s.reserve(ph.size());
  for (const YFS::Photon &k : ph) s.push_back(&k);
  std::sort(s.begin(), s.end(),
            [](const YFS::Photon *a, const YFS::Photon *b) {
              return a->E() > b->E();
            });
  h1 = s.empty() ? 0. : get(*s[0]);
  h2 = h1 + (s.size() < 2 ? 0. : get(*s[1]));
}

NLO_Base::NLO_Base() {
  p_yfsFormFact = std::make_unique<YFS::YFS_Form_Factor>();
  p_nlodipoles = std::make_unique<YFS::Define_Dipoles>();
  // p_real/p_virt/p_realvirt/p_realreal/p_vv are default-null in the header.
  // p_realreal was the one the list here forgot, so m_rrtool read an
  // uninitialised pointer whenever SetProviders had not run yet.
  m_evts = 0;
  m_recola_evts = 0;
  m_realtool = 0;
  m_realvirt = 0;
  m_looptool = 0;
  m_rrtool = 0;
  m_vvtool = 0;
  m_zeroRV = 0;
  m_zeroRR = 0;
  m_nonZeroRR = 0;
  m_zeroV = 0;
  m_nonZeroRV=0;
  m_real_hard1 = 0.;
  m_rv_hard1 = 0.;
  m_rr_hard2 = 0.;
  m_real_hard2 = 0.;
  m_rv_hard2 = 0.;
  m_zero_real_amp = 0;
  m_ceex_done = false;
  m_softRV = 0;
  m_softRR = 0;
  m_rvUnstable = 0;
  m_rvHiC = 0;
  m_rvBlowup = 0;
  m_rvBlowupRtree0 = 0;
  m_rvBlowupHiC = 0;
  m_rvBlowupSoft = 0;
  m_rvBlowupHardWide = 0;
  BookHistograms();
  if (m_check_poles == 1) {
    if (!ATOOLS::DirectoryExists(m_debugDIR_NLO))
      ATOOLS::MakeDir(m_debugDIR_NLO);
    m_histograms1d["SinglePoleCD"] = std::make_unique<Histogram>(0, 0, 25, 25);
    m_histograms1d["SinglePoleVV"] = std::make_unique<Histogram>(0, 0, 25, 25);
    m_histograms1d["DoublePoleVV"] = std::make_unique<Histogram>(0, 0, 25, 25);
    m_histograms1d["OneLoopEpsLP"] = std::make_unique<Histogram>(0, -1.5, -0.5, 50);
    m_histograms1d["OneLoopEpsYFS"] = std::make_unique<Histogram>(0, -1.5, -0.5, 50);
    m_histograms1d["RealLoopEpsLP"] = std::make_unique<Histogram>(0, -5, 0.0, 50);
    m_histograms1d["RealLoopEpsYFS"] = std::make_unique<Histogram>(0, -5, 0.0, 50);
    m_histograms1d["relativediff"] = std::make_unique<Histogram>(0, -20., -5.0, 50);
    m_histograms1d["RVSinglePoleCD"] = std::make_unique<Histogram>(0, 0, 25, 25);
    m_histograms2d["REAL_SUB"] =
        std::make_unique<Histogram_2D>(0, 0, sqrt(m_s) / 2., 200, 0, 2 * M_PI, 20);
    m_histograms2d["REAL"] =
        std::make_unique<Histogram_2D>(0, 0, sqrt(m_s) / 2., 200, 0, 2 * M_PI, 20);
  }
  if (m_rv_cancel_hist) {
    if (!ATOOLS::DirectoryExists(m_debugDIR_NLO))
      ATOOLS::MakeDir(m_debugDIR_NLO);
    m_histograms1d["RV_tot_by_logC_w"] = std::make_unique<Histogram>(0, -16., 0., 80);
    m_histograms1d["RV_tot_by_logC_n"] = std::make_unique<Histogram>(0, -16., 0., 80);
    m_histograms1d["RV_tot_by_Efrac_w"] = std::make_unique<Histogram>(0, 0., 0.5, 100);
    m_histograms1d["RV_tot_by_Efrac_n"] = std::make_unique<Histogram>(0, 0., 0.5, 100);
    m_histograms1d["RV_MEstab_all"] = std::make_unique<Histogram>(0, -2., 40., 84);
    m_histograms1d["RV_MEstab_hardwide"] = std::make_unique<Histogram>(0, -2., 40., 84);
    m_histograms1d["RV_tot_by_MEstab_w"] = std::make_unique<Histogram>(0, -2., 40., 84);
  }
}

NLO_Base::~NLO_Base() {
  WriteHistograms();
  msg_Out()<<"Total zero V: "<<m_zeroV<<std::endl;
  msg_Out()<<"Total zero RV: "<<m_zeroRV<<std::endl;
  msg_Out()<<"Total zero RR: "<<m_zeroRR<<std::endl;
  msg_Out()<<"Total non-zero RR: "<<m_nonZeroRR<<std::endl;
  msg_Out()<<"Total non-zero RV: "<<m_nonZeroRV<<std::endl;
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    int gbuf[3] = {m_softRV, m_rvUnstable, m_softRR};
    mpi->Allreduce(gbuf, 3, MPI_INT, MPI_SUM);
    m_softRV = gbuf[0];
    m_rvUnstable = gbuf[1];
    m_softRR = gbuf[2];
  }
#endif
  msg_Out()<<"Total soft RV skipped: "<<m_softRV<<std::endl;
  msg_Out()<<"Total unstable-ME RV skipped (RV_ME_MAX_RATIO): "<<m_rvUnstable<<std::endl;
  if (RVMode() == 1 && m_realvirt) {
#ifdef USING__MPI
    if (mpi->Size() > 1) {
      int rb[2] = {m_rvPoleFail, m_rvNoVirt};
      mpi->Allreduce(rb, 2, MPI_INT, MPI_SUM);
      m_rvPoleFail = rb[0];
      m_rvNoVirt = rb[1];
    }
#endif
    msg_Out()<<"RV_MODE 1: pole mismatches > 1e-6: "<<m_rvPoleFail
             <<", photons without a Born virtual: "<<m_rvNoVirt<<std::endl;
    msg_Out()<<"RV_MODE 1 (this rank): loop frame "<<RVLoopFrame()
             <<", two-frame checks "<<m_rvLoopChecked<<", dropped as unstable "
             <<m_rvLoopUnstable<<std::endl;
  }
  msg_Out()<<"Total soft RR pairs skipped: "<<m_softRR<<std::endl;
  if (m_rv_cancel_hist) {
#ifdef USING__MPI
    if (mpi->Size() > 1) {
      int buf[6] = {m_rvHiC,      m_rvBlowup,    m_rvBlowupRtree0,
                    m_rvBlowupHiC, m_rvBlowupSoft, m_rvBlowupHardWide};
      mpi->Allreduce(buf, 6, MPI_INT, MPI_SUM);
      m_rvHiC = buf[0];
      m_rvBlowup = buf[1];
      m_rvBlowupRtree0 = buf[2];
      m_rvBlowupHiC = buf[3];
      m_rvBlowupSoft = buf[4];
      m_rvBlowupHardWide = buf[5];
    }
#endif
    msg_Out()<<"RV photons with C>=1 (subtraction not cancelling): "<<m_rvHiC<<std::endl;
    msg_Out()<<"RV blow-ups |tot|>1e3*|Born|: "<<m_rvBlowup
             <<"  (of these: rtree==0: "<<m_rvBlowupRtree0
             <<", C>=1: "<<m_rvBlowupHiC
             <<", soft E/sqrt(s)<0.01: "<<m_rvBlowupSoft
             <<", HARD WIDE-ANGLE: "<<m_rvBlowupHardWide<<")"<<std::endl;
    if (m_rvBlowup>0)
      msg_Out()<<"  -> hard wide-angle fraction of blow-ups: "
               <<(100.*m_rvBlowupHardWide/m_rvBlowup)<<"%"
               <<" (if >0, instability is NOT confined to soft/collinear)"<<std::endl;
    if (m_rvBlowup>0)
      msg_Out()<<"  -> rtree==0 fraction of blow-ups: "
               <<(100.*m_rvBlowupRtree0/m_rvBlowup)<<"%"
               <<" (mechanism confirmed if ~100%)"<<std::endl;
  }
  msg_Out()<<"Total zero real amplitudes: "<<m_zero_real_amp<<std::endl;
  msg_Out()<<"Total events : "<<m_evts<<std::endl;
  ReportLoopHelicity();
}

void NLO_Base::SetProviders(YFS::Virtual *virt, YFS::Real *real,
                            YFS::RealVirtual *realvirt, YFS::RealReal *realreal,
                            YFS::VirtualVirtual *vv) {
  p_virt     = virt;
  p_real     = real;
  p_realvirt = realvirt;
  p_realreal = realreal;
  p_vv       = vv;
  m_looptool = (p_virt     != NULL);
  m_realtool = (p_real     != NULL);
  m_realvirt = (p_realvirt != NULL);
  m_rrtool   = (p_realreal != NULL);
  m_vvtool   = (p_vv       != NULL);
}

void NLO_Base::Init(Flavour_Vector &flavs, Vec4D_Vector &plab,
                    Vec4D_Vector &born) {
  m_rawbeta.clear();
  m_flavs = flavs;
  m_plab = plab;
  m_bornMomenta = born;
}

double NLO_Base::CalculateVirtual() {
  m_lhelok = false;
  m_vborn_ok = false;
  if (CeexSuppliesVirtual())
    return (m_ceexvirt - 1.) * m_born;
  if (m_eex_virt) {
    // subtract born to avoid double counting
    // already present in eex!!
    return p_dipoles->CalculateEEXVirtual() * m_born - m_born;
  }
  if (!m_looptool)
    return 0;
  double virt;
  double sub;
  p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
  CheckMassReg();
  /*
    YFS: VIRTUAL_CANONICAL_FRAME (default 0): evaluate the Born's one-loop
    in the rest frame of the beams with the beams exactly on the z axis
    (CanonicalBeamFrame, which documents why). The event's Born point after
    ISR has two beams with the same tiny p_T, and there OpenLoops' 2 -> 2
    virtual is off by up to 8% of itself (Z-pole mu mu, full EW, alpha(0):
    0.085 against 0.091 in every rotated frame). Not for CEEX_Virtual:
    helicity, whose helicity labels are the lab's.
  */
  static const bool vcanon(ATOOLS::Settings::GetMainSettings()["YFS"]
                           ["VIRTUAL_CANONICAL_FRAME"].SetDefault(0).Get<int>() != 0);
  if (vcanon && !(m_useceex && m_ceexvirtsrc == ceexvirt::helicity))
    virt = p_virt->CalcInFrame(m_plab, CanonicalBeamFrame(m_plab), m_born);
  else
    virt = p_virt->Calc(m_plab, m_born);
  if (m_check_virt_born) {
    // the provider's Born is pointlike, m_born is dressed with the pion form
    // factor, so compare against the dressed provider Born
    if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born()
                         * ExternalFormFactor(m_plab, m_flavs), 1e-6)) {
      msg_Error() << METHOD
                  << "\n Warning! Loop provider's born is different! YFS "
                     "Subtraction likely fails\n"
                  << "Loop Provider " << ":  " << p_virt->p_loop_me->ME_Born()
                  << "\nSherpa" << ":  " << m_born << std::endl
                  << "PhaseSpace Point = ";
      for (auto _p : m_plab)
        msg_Error() << _p << std::endl;
    }
  }
  if (p_virt->FailCut())
    return 0;
  if (m_virt_sub && p_virt->p_loop_me->Mode() != 1)
    sub = p_dipoles->CalculateVirtualSub();
  else
    sub = 0;
  m_virt_raw = virt; m_virt_subval = sub * m_born / m_rescale_alpha;
  m_oneloop = (virt - sub * m_born / m_rescale_alpha);
  // YFS: CEEX_Virtual: helicity - the same loop call, resolved by helicity
  if (m_useceex && m_ceexvirtsrc == ceexvirt::helicity &&
      p_virt->p_loop_me->Mode() == 0 && !IsZero(virt) && m_born != 0.)
    BuildLoopHelicityFactors(virt, sub);
  if (IsZero(virt)){
    m_zeroV++;
    return 0;
  }
  if (p_virt->p_loop_me->Mode() == 1)
    m_oneloop /= m_rescale_alpha;
  if (IsBad(m_oneloop) || IsBad(sub)) {
    msg_Error() << "YFS Virtual is NaN" << std::endl
                << "Virtual:  " << virt << std::endl
                << "Subtraction: " << sub * m_born << std::endl
                << "PhaseSpace Point: " << std::endl
                << m_plab << std::endl;
  }
  if (m_check_poles == 1) {
    if (m_virt_sub == 0)
      sub = p_dipoles->CalculateVirtualSub();
    double p1 = p_virt->p_loop_me->ME_E1() * p_virt->m_factor;
    double yfspole = p_dipoles->Get_E1();
    int ncorrect = ::countMatchingDigits(p1, -yfspole);
    double reldiff = (p1 + yfspole) / p1;
    if (!IsEqual(p1, -yfspole, 1e-4)) {
      msg_Error() << "Poles do not cancel in YFS Virtuals" << std::endl
                  << "Correct digits =  " << ncorrect << std::endl
                  << "Relative diff =  " << reldiff
                  << std::endl
                  // <<"Process =  "<<p_virt->p_loop_me->Name()<<std::endl
                  << "One-Loop Provider V eps^{-1}  = " << p1 << std::endl
                  << "Sherpa V eps^{-1} = " << yfspole << std::endl
                  << "Sherpa/One-Loop = " << yfspole / p1 << std::endl;
      return 0;
    } else {
      int i = 0;
      msg_Debugging() << std::setprecision(32);
      msg_Debugging() << "Poles cancel in YFS Virtuals to " << ncorrect
                      << " digits" << std::endl
                      << "Relative diff =  " << reldiff << std::endl;
      m_histograms1d["SinglePoleCD"]->Insert(ncorrect);
      m_histograms1d["OneLoopEpsYFS"]->Insert(log10(fabs(yfspole)));
      m_histograms1d["OneLoopEpsLP"]->Insert(log10(fabs(p1)));
      m_histograms1d["relativediff"]->Insert(log10(fabs(reldiff)));
      msg_Debugging() << std::setprecision(32)
                      << "One-Loop Provider V eps^{-1}  = " << p1 << std::endl
                      << "Sherpa V eps^{-1}  = " << yfspole << std::endl;
    }
  }
  // v of this event for the RV_MODE 1 remainder v_{n+1} - v
  if (m_born != 0. && !IsBad(m_oneloop)) {
    m_vborn = m_oneloop/m_born;
    m_vborn_ok = true;
  }
  return m_oneloop;
}

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
    if (RVMode() == 1) m_rvinfo.push_back(RVPointInfo());
    m_lastrvinfo = RVPointInfo();
    if (RRMode() == 1) m_rrphot.push_back(RRPhotonInfo());   // factor 1 if skipped
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
    if (RVMode() == 1) m_rvinfo.back() = m_lastrvinfo;
    if (RRMode() == 1) m_rrphot.back() = m_lastrrinfo;
    real += contrib;
    // YFS: REAL_COMBINE 1: this photon's factor, its bracket plus the IF
    // interference ratio m_ifi_prod picked up for it (1 without IFI_Real)
    if (m_born != 0.) {
      const double rif(m_ifireal && ifi_before != 0. && !IsBad(m_ifi_prod)
                       ? m_ifi_prod/ifi_before : 1.);
      prodw *= contrib/m_born + rif;
      if (RRMode() == 1) m_rrphot.back().factor = contrib/m_born + rif;
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
    YFS: REAL_COMBINE. 0: the O(alpha) sum 1 + sum_j delta_j over
    the event's photons, delta_j = beta_1(k_j)/(S~_j B) plus that photon's IF
    ratio - 1. 1: the product prod_j (1 + delta_j), which agrees with the sum
    at O(alpha) (identical with one photon) and adds the factorised
    beta_2 ~ beta_1 beta_1/beta_0 at O(alpha^2).
  */
  { static const int comb(ATOOLS::Settings::GetMainSettings()["YFS"]
                          ["REAL_COMBINE"].SetDefault(1).Get<int>());
    const bool useprod(comb == 1 || (comb < 0 && BornHasExchangeLine(p_bornproc, m_plab)));
    // RV_MODE 1's real-virtual is a remainder on top of this product (it
    // supplies neither beta_2 nor v x real), so the product stays on with it
    if (useprod && m_born != 0. && !IsBad(prodw)
        && (!m_rrtool || RRMode() == 1) && (!m_realvirt || RVMode() == 1))
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

double NLO_Base::CalculateReal(Vec4D k, bool raw) {
  double norm = 2. * pow(2 * M_PI, 3);
  Vec4D_Vector p(m_plab), pi(m_bornMomenta), pf(m_bornMomenta);
  dipoletype::code fluxtype;
  Vec4D kk = k;
  m_evts += 1;

  msg_Debugging() << METHOD << " raw=" << raw
                  << " k=" << k << " E=" << k.E() << " pt=" << k.PPerp() << "\n";

  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);
  fluxtype = p_nlodipoles->WhichResonant(k);

  msg_Debugging() << METHOD << " fluxtype=" << fluxtype << "\n";
  MapMomenta(p, k);

  p.push_back(k);
  /*
    YFS: ME_PROBE - for one final-state photon, the photon-lepton angle in
    units of that lepton's m/E and the lepton energy, in the event (post-
    emission legs, lab photon) and at the point handed to Calc_R, before and
    after CheckMasses; and the real it returns. Dead-cone referee against
    CEEX, which evaluates Comix on the event legs.
  */
  static const bool meprobe(ATOOLS::Settings::GetMainSettings()["YFS"]
                            ["ME_PROBE"].SetDefault(0).Get<int>() != 0);
  auto dcangle = [&](const Vec4D_Vector &q, const Vec4D &g, double &el) {
    double best(1e99); el = 0.;
    for (size_t i(2); i < q.size() && i < m_flavs.size(); ++i) {
      if (!m_flavs[i].IsChargedLepton()) continue;
      const double ct(Vec3D(q[i])*Vec3D(g)/(Vec3D(q[i]).Abs()*Vec3D(g).Abs()));
      const double th(acos(Max(-1., Min(1., ct)))/(m_flavs[i].Mass()/q[i][0]));
      if (th < best) { best = th; el = q[i][0]; }
    }
    return best; };
  double pe_ev(0.), pe_pre(0.), th_ev(-1.), th_pre(-1.);
  const bool probe_this(meprobe && PhotonIsFSR(kk) && m_photons.size() == 1);
  if (probe_this) {
    th_ev  = dcangle(m_postlab, kk, pe_ev);
    th_pre = dcangle(p, k, pe_pre);
  }
  const Vec4D_Vector p_before(p);
  CheckMasses(p, 1);

  Vec4D_Vector pp = p;
  pp.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, pp, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, pp, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, pp, m_plab);

  double r = p_real->Calc_R(p) / norm;
  if (probe_this) {
    double pe_post(0.);
    const double th_post(dcangle(p, p.back(), pe_post));
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
    YFS: REAL_BORN_PHOTON_MULTICHANNEL 1: the share of the exact real that
    belongs to the generated assignment of this photon, h_gen/sum_b h_b
    (BornPhotonChannelWeight). Only r is partitioned: each channel keeps its
    own beta_0 and its own subtraction, so the soft bracket (where no other
    channel is open and the factor is exactly 1) is unchanged.
  */
  const double r_nomc(r);
  r *= BornPhotonChannelWeight(kk);
  m_real = r;
  if (p_real->FailCut()) m_failcut = true;

  double flux;
  if (m_flux_mode == 1)
    flux = p_nlodipoles->CalculateFlux(k);
  else if (m_flux_mode == 2)
    flux = 0.5 * (p_nlodipoles->CalculateFlux(kk) + p_nlodipoles->CalculateFlux(k));
  else
    flux = p_dipoles->CalculateFlux(kk);

  /*
    Define_Dipoles::CalculateFlux forces fluxtype = initial whenever both ISR
    and FSR are on (the WhichResonant() result above is never used), so a
    FINAL-state photon received the initial-state flux (Q_X - k)^2/Q_X^2 =
    1 - x, as if it had reduced the beam energy. It has not: the Born scale
    s' is untouched by final-state emission, and the two-body phase space at
    the reduced pair mass differs from the crude one only by the muon
    velocity ratio. e+e- -> mu mu at 0.7 GeV (CMD): the real ME of every FSR
    photon was scaled by 0.70-0.78 before the subtraction, the single-FSR-
    photon events had Born+real/CEEX with median 0.75 and a 90th percentile
    of 6.4 where the two must agree event by event, and the photon spectrum
    grew a 2x bump against KKMC at E_gamma/sqrt(s) = 0.12-0.22.
    Setting flux = 1 (YFS: REAL_FSR_FLUX: 1) made the single-FSR-photon
    Born+real/CEEX WORSE (median 0.75 -> 1.67) with the old subtraction.
    With REAL_SUB_EIK 8 (the multichannel crude from the pre-FSR system)
    it is the consistent choice, and 1 is the default since 2026-09-26.
  */
  { static const int fsrflux(ATOOLS::Settings::GetMainSettings()["YFS"]
                             ["REAL_FSR_FLUX"].SetDefault(1).Get<int>());
    if (fsrflux != 0 && PhotonIsFSR(kk)) {
      if (fsrflux == 1) flux = 1.;          // no flux
      else if (fsrflux == 2) flux *= flux;  // (m_ff^2/s')^2
      /*
        4: the flux of the photon's OWN pair, (Q_D - k)^2/Q_D^2 with Q_D the
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
      else if (fsrflux == 4 && m_plab.size() == m_flavs.size()) {
        const YFS::Photon *g(FindPhoton(kk));
        if (g != nullptr && g->Dip() != nullptr) {
          const int l(g->Dip()->Left()), r(g->Dip()->Right());
          if (l >= 2 && r >= 2 && l < (int)m_plab.size() && r < (int)m_plab.size()) {
            const Vec4D Qd(m_plab[l] + m_plab[r]);
            const double q2(Qd.Abs2());
            if (q2 > 0.) flux = (Qd - kk).Abs2()/q2;
          } } }
    } }
  // CalculateRealSub is the plain eikonal; the symmetry factor is on r above
  double subloc = p_nlodipoles->CalculateRealSub(k);
  double subb   = p_dipoles->CalculateRealSubEEX(kk);
  const double flux0(flux), subloc0(subloc), subb0(subb);   // SUB8_TRACE
  /*
    Which eikonal beta_1 subtracts (YFS: REAL_SUB_EIK; default 8 since
    2026-09-26, see mode 8 below). 0: the coherent one at the mapped,
    post-emission point (the old default). 1: the event's crude,
    S~_II + S~_FF on the pre-emission legs - the density the photon was
    generated with, so the weight is r flux/(S~ B) exactly. 2: the coherent
    eikonal on the pre-emission (Born) legs of the event, the legs the form
    factor exponent is built on.
  */
  { static const int subeik(ATOOLS::Settings::GetMainSettings()["YFS"]
                            ["REAL_SUB_EIK"].SetDefault(8).Get<int>());
    if (subeik == 1) subloc = subb;
    else if (subeik == 2) subloc = p_dipoles->CalculateRealSub(kk);
    /*
      3: the event's crude for FINAL-state photons only. With the
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
    else if (subeik == 3 && PhotonIsFSR(kk)) subloc = subb;
    /*
      4: for FINAL-state photons, subtraction AND denominator on the
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
    /*
      5: for FINAL-state photons, the COHERENT eikonal on the pre-emission
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
    else if (subeik == 5 && PhotonIsFSR(kk)) {
      // The same two pieces m_ifi_prod is built from (below): the incoherent
      // crude on the generation legs plus the initial-final interference on
      // those legs. CalculateRealSub(kk) is NOT this - it sits on the
      // post-emission dipole momenta and reproduced mode 4 exactly.
      const double ifg(p_dipoles->CalculateRealSubIF(kk));
      if (!IsBad(ifg) && !IsZero(subb)) subloc = subb + ifg;
    }
    /*
      6: mode 3 with the crude carrying the Born of each ASSIGNMENT. A hard
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
    else if (subeik == 6 && PhotonIsFSR(kk)) {
      double sII(0.), sFF(0.), br(1.);
      if (p_dipoles->HasDipoleII()) {
        YFS::Dipole &D(p_dipoles->GetDipoleII());
        sII = D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1)); }
      for (auto &D : p_dipoles->GetDipoleFF())
        sFF += D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1));
      if (sII + sFF > 0. && ReducedBornISR(kk, br)) {
        const double s6(sII*br + sFF);
        if (s6 > 0. && !IsBad(s6)) { subloc = subb = s6; }
      } else subloc = subb;
    }
    /*
      7: the MULTICHANNEL crude. The generator reaches a given final state
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
      ISR-labelled half is applied where the reduced-point crude is built.
    */
    else if (subeik == 7 && PhotonIsFSR(kk)) {
      double sII(0.), sFF(0.), bI(1.), s2(1.);
      if (p_dipoles->HasDipoleII()) {
        YFS::Dipole &D(p_dipoles->GetDipoleII());
        sII = D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1)); }
      for (auto &D : p_dipoles->GetDipoleFF())
        sFF += D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1));
      const double ifg(p_dipoles->CalculateRealSubIF(kk));
      if (sFF > 0. && ReducedBornAlt(kk, -1, bI, s2) && s2 > 0.) {
        const double s7(sFF + sII*bI/s2);
        if (s7 > 0. && !IsBad(s7)) { subb = s7; subloc = s7 + (IsBad(ifg) ? 0. : ifg); }
      } else { subloc = subb + (IsBad(ifg) ? 0. : ifg); }
      flux = 1.;
    }
    /*
      8: mode 7 with the alternative channel built from the right system.
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
    */
    else if (subeik == 8 && PhotonIsFSR(kk) && p_bornproc != nullptr
             && Sub8Part() != 2) {
      double sII(0.), sFF(0.);
      if (p_dipoles->HasDipoleII()) {
        YFS::Dipole &D(p_dipoles->GetDipoleII());
        sII = D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1)); }
      /*
        YFS: SUB8_FSR_LEGS (diagnostic, 2026-09-27, default 0 = unchanged).
        1: the FSR channel density on the POST-emission legs times the pair
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
      static const int fsrlegs(ATOOLS::Settings::GetMainSettings()["YFS"]
                               ["SUB8_FSR_LEGS"].SetDefault(0).Get<int>());
      /*
        2 and 3 (NOTES-deadcone-crude-2026-09-28.md): the FSR channel density
        the generator really has, post-emission eikonal times pair flux,
        so that the one-photon weight is r/(g B), the density CEEX divides by
        (hard half of SUB8_FSRSUB 4: g (1 + S_IF/crude), on the new g).
        Unlike 1, the SOFT half keeps the generation-leg scale (gsoft below):
        it is the point's coherent current, the soft limit of r, and must not
        be rescaled by a post/pre eikonal ratio, which leaves soft photons
        next to a lepton kicked by a companion with O(1) brackets (1 doubled
        YFS.NLO's error at the Z pole).
        2: the EVENT's density, GeneratorFSRDensity - S~ on the event's
           post-emission legs times the dipole's flux F_D shared among its
           photons. Right at one photon; with companions it describes a
           different emission from the one |M_1|^2 is evaluated for (the
           REAL_FSR_MAP 2 point re-absorbs the companions, which moves the
           lepton by more than its dead cone for electrons): Z-pole e e,
           YFS.NLO fiducial error x4, CEEX/NLO 1.034 at s'/s 0.85-0.93.
        3 (recommended): the POINT's density, PointFSRDensity - S~ on the
           legs of the (n+1)-body point times its own single-emission flux.
           Identical to 2 at one photon; with companions it is the density
           of the emission the point describes, and it reduces to the old
           generation-leg crude where the companions dominate the recoil.
      */
      double sFFgen(0.);
      for (auto &D : p_dipoles->GetDipoleFF())
        sFFgen += D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1));
      if (fsrlegs == 2) {
        if (!GeneratorFSRDensity(kk, *p_dipoles, m_postlab, sFF)) sFF = sFFgen;
      }
      else if (fsrlegs == 3) {
        if (!PointFSRDensity(k, *p_dipoles, pp, sFF)) sFF = sFFgen;
      }
      else if (fsrlegs == 1 && m_postlab.size() == m_plab.size()) {
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
      const double ifg(p_dipoles->CalculateRealSubIF(kk));
      const Vec4D Q(PreFSRSystem()), R(Q - kk);
      const double s2(Q.Abs2() > 0. ? R.Abs2()/Q.Abs2() : 0.);
      double bI(0.);
      if (!(s2 > 0.) || !PreFSRBornRatio(R, bI)) bI = 0.;   // no ISR channel
      const double isrchannel(s2 > 0. ? sII*bI/s2 : 0.);
      const double g(sFF + isrchannel);
      // the scale of the soft (point) subtraction: see SUB8_FSR_LEGS 2, 3
      const double gsoft(fsrlegs >= 2 ? sFFgen + isrchannel : g);
      /*
        Which subtraction (YFS: SUB8_FSRSUB):
        0: mode 5's, on the generation (pre-FSR) legs, scaled with g. The
           one-photon weight is then exactly r/(g_I + g_F) (MODE: FSR cone
           identity 1.00), but r's soft limit sits on the legs of the point
           it is evaluated at, not on the generation legs: with a hard ISR
           photon (250 GeV radiative return) the FSR point's beams are the
           generator's reduced beams, tilted by the ISR transverse momentum,
           and next to a hard FSR companion its leptons are the kicked ones.
           Soft FSR photons there kept O(1) brackets (mean +0.18 per photon
           below x = 1e-3, single photons up to 250), YFS.NLO +6% at 250 GeV
           from this half alone.
        1: the coherent current at that point (the default's subloc),
           scaled with g: soft brackets vanish, the identity for hard photons
           does not hold (Z pole +4.4% on Born+real).
        4 (default): 1 for soft photons, 0 for hard ones,
           w = 1/(1 + (y/y0)^2), y = 2 k.Q_pre/Q_pre^2, y0 = YFS: SUB8_Y0
           (0.01). 250 GeV FSR half 1.062 -> 0.992 of CEEX, Z pole and the
           cone identity unchanged.
      */
      static const int fsub(ATOOLS::Settings::GetMainSettings()["YFS"]
                            ["SUB8_FSRSUB"].SetDefault(4).Get<int>());
      if (g > 0. && !IsBad(g) && subb > 0.) {
        const double cru(subb);
        const double sgen(g*(1. + (IsBad(ifg) ? 0. : ifg/cru)));
        if (fsub == 1) subloc *= gsoft/cru;
        else if (fsub == 4) {
          /*
            4: the point's coherent subtraction for soft photons, mode 5's
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
        flux   = 1.;
      }
      m_sub8[0] = sII; m_sub8[1] = sFF; m_sub8[2] = bI; m_sub8[3] = s2;
    }
    else if (subeik == 4 && PhotonIsFSR(kk) && m_postlab.size() == m_plab.size()) {
      double s(0.);
      if (p_dipoles->HasDipoleII()) {
        YFS::Dipole &D(p_dipoles->GetDipoleII());
        s += D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1)); }
      for (auto &D : p_dipoles->GetDipoleFF()) {
        const int l(D.Left()), r(D.Right());
        if (l >= 2 && r >= 2 && l < (int)m_postlab.size() && r < (int)m_postlab.size())
          s += D.Eikonal(kk, m_postlab[l], m_postlab[r]); }
      if (s > 0. && !IsBad(s)) { subloc = subb = s; }
    } }
  /*
    REAL_FSR_FLUX: 3 - the crude a FINAL-state photon is divided by carries
    F = (q + k)^2/q^2 on its final-state part, q the post-emission pair, and
    r carries no flux. This is the crude CEEX divides by (validated against
    KKMC): on single-FSR-photon events K rho_crude(CEEX)/((S~_II + F S~_FF) B)
    has median 0.96 where K rho_crude/(S~ B) has 1.15 and a tail to 2.1.
  */
  { static const int fsrflux(ATOOLS::Settings::GetMainSettings()["YFS"]
                             ["REAL_FSR_FLUX"].SetDefault(1).Get<int>());
    if (fsrflux == 3 && PhotonIsFSR(kk) && m_postlab.size() == m_plab.size()) {
      Vec4D q; for (size_t i = 2; i < m_postlab.size(); ++i) q += m_postlab[i];
      const double q2(q.Abs2());
      if (q2 > 0.) {
        const double F((q + kk).Abs2()/q2);
        double sII(0.), sFF(0.);
        if (p_dipoles->HasDipoleII()) {
          YFS::Dipole &D(p_dipoles->GetDipoleII());
          sII = D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1)); }
        for (auto &D : p_dipoles->GetDipoleFF())
          sFF += D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1));
        if (sII + F*sFF > 0.) { subb = sII + F*sFF; flux = 1.; m_eikeex = subb; }
      } } }
  m_wifterms.push_back(subb != 0. ? subloc/subb : 1.);
  m_eikeex = subb;
  m_subloc = subloc;

  msg_Debugging() << METHOD << " r=" << r << " flux=" << flux
                  << " (mode=" << m_flux_mode << ")"
                  << " subloc=" << subloc << " subb=" << subb
                  << " born=" << m_born << " alpha=" << m_rescale_alpha << "\n";

  if (!CheckMomentumConservation(p)) {
    msg_Debugging() << METHOD << " momentum conservation failed"
                    << " k.E=" << k.E() << " dip_mass=" << (p[2]+p[3]).Mass() << "\n";
    msg_Error() << "Momentum Conservation fails in " << METHOD << "\n";
    if (m_isr_debug || m_fsr_debug) {
      m_histograms1d["k_E"]->Insert(k.E());
      m_histograms1d["k_pt"]->Insert(k.PPerp());
      m_histograms1d["dip_mass"]->Insert((p[2] + p[3]).Mass());
    }
    return 0;
  }

  if ((p[2] + p[3]).Mass() < massmin)
    massmin = (p[2] + p[3]).Mass();

  if (m_isr_debug || m_fsr_debug) {
    m_histograms1d["k_E_pass"]->Insert(k.E());
    m_histograms1d["k_pt_pass"]->Insert(k.PPerp());
    m_histograms1d["dip_mass_pass"]->Insert((p[2] + p[3]).Mass());
  }

  if (IsZero(r)) {
    msg_Debugging() << METHOD << " r=0, returning 0\n";
    m_zero_real_amp++;
    return 0;
  }
  if (IsBad(r) || IsBad(flux)) {
    msg_Debugging() << METHOD << " bad point: r=" << r << " flux=" << flux << "\n";
    msg_Error() << "Bad point for YFS Real\n"
                << "  Real ME : " << r << "\n"
                << "  Flux    : " << flux << "\n";
    return 0;
  }

  double tot;
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
  double subb_loc(subb);
  if (m_map_reduced) {
    subb_loc = 0.;
    if (p_nlodipoles->HasDipoleII()) {
      YFS::Dipole &D(p_nlodipoles->GetDipoleII());
      subb_loc += D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
    }
    for (auto &D : p_nlodipoles->GetDipoleFF())
      subb_loc += D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
    if (IsZero(subb_loc) || IsBad(subb_loc)) subb_loc = subb;
  }
  // the denominator as it would be without REAL_SUB_EIK's changes to it
  const double denom0(m_map_reduced ? subb_loc : subb0);
  /*
    REAL_SUB_EIK 7, ISR-labelled photon: the final-state channel term of the
    crude carries its own Born and the ISR flux, f = flux B_F/B with B_F the
    Born at (P + k)^2 (the photon returned to beams and pair). At one photon
    r flux/(B (S~_II + f S~_FF)) = r/(g_I + g_F). The same (f - 1) S~_FF is
    added to the subtraction, which keeps the soft bracket at its old value.
  */
  { static const int se7(ATOOLS::Settings::GetMainSettings()["YFS"]
                         ["REAL_SUB_EIK"].SetDefault(8).Get<int>());
    if (se7 == 7 && !PhotonIsFSR(kk)) {
      double bF(1.), s2(1.);
      if (ReducedBornAlt(kk, +1, bF, s2)) {
        const double f(flux*bF);
        double sFFpt(0.);
        if (m_map_reduced) {
          for (auto &D : p_nlodipoles->GetDipoleFF())
            sFFpt += D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
        } else {
          for (auto &D : p_dipoles->GetDipoleFF())
            sFFpt += D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1));
        }
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
    } }
  /*
    REAL_SUB_EIK 8, ISR-labelled photon: the FSR channel term of the crude,
    S~_FF with the Born at (Q_pre + k)^2 in the ISR photon's own units
    (times its flux), the whole subtraction scaled with the denominator.
  */
  { static const int se8(ATOOLS::Settings::GetMainSettings()["YFS"]
                         ["REAL_SUB_EIK"].SetDefault(8).Get<int>());
    static const int t8(ATOOLS::Settings::GetMainSettings()["YFS"]
                        ["SUB8_TRACE"].SetDefault(0).Get<int>());
    static long n8(0);
    const bool fsr(PhotonIsFSR(kk));
    double bF(-1.), f(-1.), scale(1.);
    if (se8 == 8 && !fsr && p_bornproc != nullptr && Sub8Part() != 1) {
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
      if (m_map_reduced && pp.size() == m_flavs.size() && Q.Abs2() > 0.) {
        Vec4D Rp(k);
        for (size_t i(2); i < pp.size(); ++i) Rp += pp[i];
        const double r2(Rp.Abs2());
        if (r2 > 0.) R = sqrt(r2/Q.Abs2())*Q;
      }
      bF = 0.;
      if (!PreFSRBornRatio(R, bF)) bF = 0.;
      f = flux*bF;
      double sFFpt(0.);
      if (m_map_reduced) {
        for (auto &D : p_nlodipoles->GetDipoleFF())
          sFFpt += D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
      } else {
        for (auto &D : p_dipoles->GetDipoleFF())
          sFFpt += D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1));
      }
      static const int isrred(ATOOLS::Settings::GetMainSettings()["YFS"]
                              ["SUB8_ISR_REDUCED"].SetDefault(1).Get<int>());
      const double nb(subb_loc + (f - 1.)*sFFpt);
      if ((isrred || !m_map_reduced) && !IsBad(nb) && nb > 0. && subb_loc > 0.) {
        scale = nb/subb_loc;
        subb_loc = nb;
        subloc  *= scale;
        if (!m_map_reduced) subb *= scale;
      }
    }
    if (se8 == 8 && t8 && n8 < t8) { ++n8;
      double me2ref(-1.);
      { Vec4D_Vector pr; if (PreFSRBornPoint(PreFSRSystem(), pr)) me2ref = BornME2At(pr); }
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
      { double thl(1e9);
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
        o<<" thl="<<thl<<" xisrmax="<<xisr<<" xfsrother="<<xfsr; }
      // which legs the FF eikonal sits on: generation (GetBornMomenta),
      // pre-FSR lab (m_plab), post-FSR lab (m_postlab)
      { double sg(0.), sl(0.), sp(0.), dmax(0.);
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
        o<<" sFFgen="<<sg<<" sFFpre="<<sl<<" sFFpost="<<sp<<" |dBorn|="<<dmax; }
      o<<"\n";
      std::cerr<<o.str(); }
    m_sub8[0] = m_sub8[1] = m_sub8[2] = m_sub8[3] = -1.; }
  m_lastflux = flux; m_lastsubloc = subloc;
  m_lastdenom = (m_submode == submode::local ? subloc
                 : m_submode == submode::global ? subb_loc : subb);
  m_lastsubloc0 = subloc0;
  m_lastdenom0 = (m_submode == submode::local ? subloc0
                  : m_submode == submode::global ? denom0 : subb0);
  if (m_submode == submode::local)
    tot = (r * flux - subloc * m_born / m_rescale_alpha) / subloc;
  else if (m_submode == submode::global)
    tot = (r * flux - subloc * m_born / m_rescale_alpha) / subb_loc;
  else if (m_submode == submode::off)
    tot = (r * flux) / subb;
  else
    msg_Error() << METHOD << " unknown YFS subtraction mode " << m_submode << "\n";

  /*
    YFS: REAL_ALPHA0 (default 1 since 2026-09-27). The real above is Comix's, with
    the model's alpha on every photon (1/131.9 under G_mu), and the
    subtraction is raised to match it (/m_rescale_alpha), while beta_0 and
    the eikonals carry alpha(0) (USE_MODEL_ALPHA 0). beta_1 is then the hard
    remainder at alpha_model: where the generated channel's S~ B exceeds
    the exact real by far (a Born near its t-channel pole), the one-photon
    factor tends to 1 - alpha_model/alpha(0) = -0.0387, not to 0. CEEX
    rescales its Comix reals to alpha(0) (Ceex_Base::ComixPhotonCoupling).
    1 does the same here: beta_1 -> m_rescale_alpha beta_1, the whole real
    correction at alpha(0), as in CEEX. Measured with the same seed, 100k
    events, 0 -> 1: Z-pole mu mu fiducial CEEX/YFS.NLO 1.0056 -> 1.0044,
    m_ff at 60 GeV 1.111 -> 1.079; Bhabha 1.0071 -> 1.0046 and 1.035 ->
    1.007; YFS.NLO rises (its real correction is negative there). 0 restores
    the alpha_model remainder. Exactly 1 with USE_MODEL_ALPHA 1.
  */
  { static const int a0(ATOOLS::Settings::GetMainSettings()["YFS"]
                        ["REAL_ALPHA0"].SetDefault(1).Get<int>());
    if (a0) tot *= m_rescale_alpha; }

  // WEIGHT_PROBE bookkeeping for the hardest photon of the event: the
  // exact-over-crude ratio before the multichannel share, the subtraction
  // over the denominator (the "1 - sub" floor), and G.
  if (kk.E() >= m_bpmc_hardx) {
    const double den(m_submode == submode::local ? subloc : subb_loc);
    m_bpmc_hardx = kk.E();
    m_bpmc_hardG = m_bpmc_lastG;
    m_bpmc_hardR = (m_born != 0. && den != 0.) ? r_nomc*flux/(m_born*den) : 0.;
    m_bpmc_hardsub = (den != 0.) ? subloc/(m_rescale_alpha*den) : 0.;
  }

  { static const double tr(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_TRACE"].Get<double>());
    if (tr > 0.) {
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
      double sII(0.), sFF(0.);
      if (p_dipoles->HasDipoleII()) {
        YFS::Dipole &D(p_dipoles->GetDipoleII());
        sII = D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1)); }
      for (auto &D : p_dipoles->GetDipoleFF())
        sFF += D.Eikonal(kk, D.GetBornMomenta(0), D.GetBornMomenta(1));
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
                 <<" Bpp/B="<<(m_born>0.?BornME2At(pp)/m_born:-1.)
                 <<" born="<<m_born
                 <<" beta1/S~="<<tot<<" beta1/(S~B)="<<(m_born!=0.?tot/m_born:0.)
                 <<" failcut="<<(p_real->FailCut()?1:0)
                 <<"\n     klab="<<kk<<" kmap="<<k
                 <<"\n     Pa="<<m_bornMomenta[0]<<" Pb="<<m_bornMomenta[1]
                 <<" IIborn0="<<p_dipoles->GetDipoleII().GetBornMomenta(0)
                 <<"\n     pa_j="<<p[0]<<" pb_j="<<p[1]<<" Qlab="<<Q
                 <<"\n";
    } }
  { static const bool ds(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>()!=0);
    if (ds) {
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
    } }

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
      {
        const int ib = std::min(4, (int)(10.*kk.E()/sqrt(m_s)));
        if (ib >= 0) {
          ++m_ifi_x_n[ib];
          m_ifi_x_r[ib] += ratio;
          const double ex = subloc/(m_rescale_alpha*subb);
          m_ifi_x_e[ib] += IsBad(ex) ? ratio : ex;
        }
      }
      if (ratio < m_ifi_min) m_ifi_min = ratio;
      if (ratio > m_ifi_max) m_ifi_max = ratio;
    }
    msg_Debugging() << METHOD << " IFI_Real ratio=" << ratio << "\n";
  }

  msg_Debugging() << METHOD << " submode=" << m_submode
                  << " r*flux=" << r*flux
                  << " sub=" << subloc * m_born / m_rescale_alpha
                  << " tot=" << tot << "\n";

  if (m_isr_debug || m_fsr_debug) {
    if (m_isr_debug)
      m_histograms2d["Real_Flux"]->Insert(
          flux, sqrt(p_dipoles->GetDipoleII().Sprime()));
  }

  if (m_no_subtraction) {
    msg_Debugging() << METHOD << " no_subtraction: returning r/subloc=" << r/subloc << "\n";
    return r / subloc;
  }

  if (IsBad(tot)) {
    msg_Debugging() << METHOD << " tot is NaN/Inf"
                    << " r*flux=" << r*flux
                    << " subloc*born=" << subloc*m_born
                    << " subb=" << subb << "\n";
    msg_Error() << "NLO real is NaN\n"
                << "  R        : " << r << "\n"
                << "  Local  S : " << subloc * m_born << "\n"
                << "  Global S : " << subb << "\n";
  }

  if (m_isr_debug || m_fsr_debug) {
    m_histograms2d["IFI_EIKONAL"]->Insert(k.Y(), k.PPerp(),
                                          p_nlodipoles->CalculateRealSubIF(k));
    m_histograms2d["REAL_SUB"]->Insert((p[0] + p[1]).Mass(), k.E(), tot / m_born);
    m_histograms2d["REAL"]->Insert(k.E(), k.Theta(), r);
    m_histograms2d["REAL_SUB"]->Insert(k.E(), k.Theta(), tot);
  }

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

  if (raw) {
    double rawval = r * flux - subloc * m_born / m_rescale_alpha;
    msg_Debugging() << METHOD << " raw: returning " << rawval << "\n";
    return rawval;
  }

  msg_Debugging() << METHOD << " returning tot=" << tot << "\n";
  if (RVMode() == 1 && m_born != 0. && m_lastdenom != 0.) {
    // RV_MODE 1: the point and the photon's exact-over-crude ratio, in the
    // same units and coupling as tot (REAL_ALPHA0's kappa included)
    static const int a0(ATOOLS::Settings::GetMainSettings()["YFS"]
                        ["REAL_ALPHA0"].SetDefault(1).Get<int>());
    m_lastrvinfo.ok = true;
    m_lastrvinfo.p = p;
    m_lastrvinfo.rho = (a0 ? m_rescale_alpha : 1.)*r*flux/(m_lastdenom*m_born);
  }
  if (RRMode() == 1 && m_born != 0. && m_lastdenom != 0. && !IsBad(tot)) {
    m_lastrrinfo.ok = true;
    m_lastrrinfo.delta = tot/m_born;
    m_lastrrinfo.flux = flux;
    m_lastrrinfo.denom = m_lastdenom;
    m_lastrrinfo.subloc = m_lastsubloc/m_lastdenom;
    { static const int a0rr(ATOOLS::Settings::GetMainSettings()["YFS"]
                            ["REAL_ALPHA0"].SetDefault(1).Get<int>());
      m_lastrrinfo.rho = (a0rr ? m_rescale_alpha : 1.)*r*flux/(m_lastdenom*m_born); }
    m_lastrrinfo.pt = p;
  }
  return tot;
}

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
  YFS: RV_LOOP_FRAME - the frame RV_MODE 1 hands the (n+1)-point loop to
  OpenLoops in.
  0: the lab, the real's point as it is (the legacy RV's frame). Beams on z
     up to the rounding of the rotate-back (|p_T|/E ~ 1e-33): the
     hp_mode 1 offset of 0.54% on e+e- -> mu mu gamma, and photons within
     ~3e-6 rad of a beam lose precision on the axis.
  1 (default): CanonicalBeamFrame, the beams' rest frame turned by a fixed
     generic rotation. 
  2: the beams' rest frame with the beams exactly on z, then tilted by
     1e-3 ra
  3: the beams' rest frame with the beams exactly on z (no tilt): stable for
     every nu nu gamma point seen
  No single frame is right everywhere for photons deep in the electron's dead
  cone; YFS: RV_LOOP_CHECK 1 compares two (see RealVirtualFactor).
*/

static Vec4D_Vector TiltedBeamFrame(const Vec4D_Vector &p)
{
  Vec4D_Vector q(BeamRestFrameOnZ(p));
  static const double c(cos(1e-3)), s(sin(1e-3));
  for (Vec4D &v : q) v = Vec4D(v[0], c*v[1] + s*v[3], v[2], -s*v[1] + c*v[3]);
  return q;
}

static Vec4D_Vector RVLoopPoint(const Vec4D_Vector &p, int frame)
{
  switch (frame) {
  case 0: return p;
  case 2: return TiltedBeamFrame(p);
  case 3: return BeamRestFrameOnZ(p);
  default: return CanonicalBeamFrame(p);
  }
}

/*
  YFS: RV_LOOP_CHECK (default 0). 1: for a photon within RV_LOOP_CHECK_ANGLE
  (default 1e-4 rad, ~10 m_e/E at the Z pole) of a beam, evaluate V_fin/T a
  second time, in frame 2 if RV_LOOP_FRAME is 3 and in frame 3 otherwise,
  and drop the photon's real-virtual (counted, reported at the end) when the
  two differ by more than RV_LOOP_CHECK_TOL (default 1e-3). The remainder of
  such photons is O(1e-3 .. 1e-2) of |M_1|^2 at most; OpenLoops' value there
  can be O(1e4) in one frame.
*/
static bool RVLoopCheckNeeded(const Vec4D_Vector &p)
{
  static const int on(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_LOOP_CHECK"]
                      .SetDefault(0).Get<int>());
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
    const Vec4D_Vector q(RVLoopPoint(pin, RVLoopFrame() == 3 ? 2 : 3));
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
  RV_MODE 1: the Born's v with the subtraction on the Born point's OWN legs,
  v_B = V_fin(m_plab)/B - B~(m_plab legs)/kappa, from the loop CalculateVirtual
  already evaluated. The event's v (m_vborn) subtracts B~
  on the legs the YFS form factor is built on: the II dipole and the initial
  legs of the IF dipoles at the FULL beams (YFS_Handler::MakeYFS), while the
  loop is at m_plab, whose beams are reduced by the ISR photons. For a
  one-photon event the (n+1)-body point's legs are those event legs, and
  v_{n+1} - v would be exact given that convention; but with a hard ISR
  companion the point of a soft photon sits at the reduced beams, and
  v_{n+1} - v -> B~(reduced legs) - B~(full-beam legs) != 0 as the photon
  goes soft: 1e-4 .. 6e-3 per soft photon, i.e. RV_SOFT_CUT dependence.
  v_{n+1} - v_B -> 0 for every photon (both IR subtracted on their own legs:
  the object the standalone referee olref.py computes), and differs from the
  exact one-photon identity only by rho (v - v_B), the Born virtual's own
  leg convention.
*/
bool NLO_Base::BornVirtualOnOwnLegs(double &vb) {
  if (!m_looptool || m_born == 0.) return false;
  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
  const double sub(p_nlodipoles->CalculateVirtualSub());
  vb = m_virt_raw/m_born - sub/m_rescale_alpha;
  return !IsBad(vb);
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
  YFS: RV_PHOTON_CT (default 0) - the emitted photon's charge renormalisation
  in the real-virtual. With REAL_ALPHA0 1 the real couples the photon with
  alpha(0), but the loop provider renormalises it in the input scheme; in a
  G_mu card (OpenLoops, Z pole) v_{n+1} - v -> -0.0075 as the photon goes
  soft. The counterterm c is added to every v_{n+1} (one external photon).
  0: none (right for an alpha(0) card, c = 0 there).
  1: calibrated, once per run: c = -(v_{n+1} - v_B) at points with a photon
     of x = 1e-4 and 2e-4 added to the first event's Born (SoftPhotonPoint),
     extrapolated linearly to x = 0. It forces the soft limit to 0 by
     construction, so a soft-limit test is then no test, and it is taken at
     one point of one process.
  2: analytic (AnalyticPhotonCounterterm, Photon_Counterterm.C): the on-shell
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
  Real_Generator's real included (-8% on sigma_fid); 1 and 2 keep the 22
  registration.
*/
static int PhotonCTMode() {
  static const int m(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PHOTON_CT"]
                     .SetDefault(0).Get<int>());
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
  if (PhotonCTMode() == 2) { m_rvct = AnalyticPhotonCounterterm(); return; }
  if (PhotonCTMode() != 1) return;
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
    const double Keps((m_rrtool && RRMode() == 1) ? m_rr_soft_cut*sqrt(m_s) : 0.5*sqrt(m_s));
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
  const double Keps((m_rrtool && RRMode() == 1) ? m_rr_soft_cut*sqrt(m_s) : 0.5*sqrt(m_s));
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
  if (RVMode() == 1) {
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
  if (m_flux_mode == 1)
    flux = p_nlodipoles->CalculateFlux(k);
  else if (m_flux_mode == 2)
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

void NLO_Base::ProbeRealVirtual(const Vec4D &kk, const Vec4D_Vector &p,
                                double r, double rtree, double subloc,
                                double eikloc, double flux, double subb,
                                double tot) {
  PHASIC::Virtual_ME2_Base *lme(p_realvirt->p_loop_me.get());
  const double lt(lme->ME_Finite()*p_realvirt->m_factor);
  const double e1lt(lme->ME_E1()*p_realvirt->m_factor);
  const double kap(m_rescale_alpha);
  const double e1ev(p_nlodipoles->Get_E1());
  // the same subtraction on the (n+1)-body point's own legs
  Vec4D_Vector pp(p);
  pp.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, pp, pp);
  p_nlodipoles->MakeDipoles(m_flavs, pp, pp);
  p_nlodipoles->MakeDipolesIF(m_flavs, pp, pp);
  const double subpp(p_nlodipoles->CalculateRealVirtualSubEps(p.back()));
  const double e1pp(p_nlodipoles->Get_E1());
  const double v(m_born != 0. ? m_oneloop/m_born : 0.);
  std::ostringstream o;
  o<<std::setprecision(10)<<"@@@ RVPROBE x="<<2.*kk.E()/sqrt(m_s)
   <<" fsr="<<(PhotonIsFSR(kk)?1:0)<<" nph="<<m_photons.size()
   <<" lt="<<lt<<" Trv/rtree="<<(rtree != 0. ? r/(lt*rtree) : 0.)
   <<" sub_ev="<<subloc/kap<<" sub_pt="<<subpp/kap
   <<" v="<<v<<" subB="<<(m_born != 0. ? m_virt_subval/m_born : 0.)
   <<" vraw="<<(m_born != 0. ? m_virt_raw/m_born : 0.)
   <<" dv_ev="<<lt - subloc/kap - v<<" dv_pt="<<lt - subpp/kap - v
   <<" eikB/rtree="<<(rtree != 0. ? eikloc*m_born/(kap*rtree) : 0.)
   <<" flux="<<flux<<" subb="<<subb<<" tot/B="<<(m_born != 0. ? tot/m_born : 0.)
   <<" E1lt="<<e1lt<<" E1ev="<<e1ev/kap<<" E1pt="<<e1pp/kap
   <<" irs="<<lme->IRscale()<<" kap="<<kap
   <<" Tloop="<<lme->ME_Born()<<" rtree="<<rtree*2.*pow(2.*M_PI,3)
   <<"\n   p=";
  for (const Vec4D &q : p) o<<q<<" ";
  o<<"\n   plab=";
  for (const Vec4D &q : m_plab) o<<q<<" ";
  o<<"\n";
  std::cerr<<o.str();
}

/*
  RR_MODE 1: the point pt (legs, then photons) with photon `which` (index
  into pt) taken out and its momentum absorbed by the final state: the final
  legs, in their rest frame, rescaled (masses kept) to the invariant mass
  (Q_f + k)^2 and boosted onto Q_f + k. Beams and the other photons are
  untouched. A soft reduction: as k -> 0 the result tends to pt without k.
*/
static bool AbsorbPhoton(const Vec4D_Vector &pt, size_t nlegs, size_t which,
                         Vec4D_Vector &out)
{
  if (which < nlegs || which >= pt.size() || nlegs < 3) return false;
  const Vec4D k(pt[which]);
  Vec4D Qf;
  for (size_t i(2); i < nlegs; ++i) Qf += pt[i];
  const Vec4D T(Qf + k);
  const double MT2(T.Abs2());
  if (!(MT2 > 0.) || !(Qf.Abs2() > 0.)) return false;
  const double MT(sqrt(MT2));
  Poincare rest(Qf);
  Vec4D_Vector f;
  std::vector<double> m2;
  for (size_t i(2); i < nlegs; ++i) {
    Vec4D v(pt[i]);
    rest.Boost(v);
    f.push_back(v);
    m2.push_back(Max(0., v.Abs2()));
  }
  double msum(0.);
  for (double m : m2) msum += sqrt(m);
  if (!(MT > msum)) return false;
  double xi(1.);
  for (int it(0); it < 60; ++it) {
    double E(0.), dE(0.);
    for (size_t i(0); i < f.size(); ++i) {
      const double p2(Vec3D(f[i]).Sqr()), e(sqrt(xi*xi*p2 + m2[i]));
      E += e;
      if (e > 0.) dE += xi*p2/e;
    }
    if (!(dE > 0.)) return false;
    const double step((E - MT)/dE);
    xi -= step;
    if (std::abs(step) < 1e-15*Max(1., xi)) break;
  }
  if (!(xi > 0.)) return false;
  out.assign(pt.begin(), pt.begin() + 2);
  Poincare onto(T);
  for (size_t i(0); i < f.size(); ++i) {
    const Vec3D p3(xi*Vec3D(f[i]));
    Vec4D v(sqrt(p3.Sqr() + m2[i]), p3);
    onto.BoostBack(v);
    out.push_back(v);
  }
  for (size_t i(nlegs); i < pt.size(); ++i) if (i != which) out.push_back(pt[i]);
  return true;
}

/*
  RR_MODE 1: the point pt with final-state photon `which` recombined with the
  leg of its radiating pair it is closer to in angle (the emitter; the other
  leg is the spectator): in the rest frame of emitter + spectator + photon
  the two legs are rebuilt back to back along emitter + photon, on shell.
  Only the pair changes, and the emitter keeps the direction of
  emitter + photon, so the collinear structure of the OTHER photon is kept
  (absorbing into the whole final state moved the legs by the hard photon's
  recoil, and the eikonals of a collinear companion by factors of 50).
*/
static bool RecombineWithEmitter(const Vec4D_Vector &pt, size_t which, int l, int r,
                                 const ATOOLS::Flavour_Vector &fl, Vec4D_Vector &out)
{
  if (which >= pt.size() || l < 2 || r < 2 || l == r
      || l >= (int)fl.size() || r >= (int)fl.size()) return false;
  const Vec4D k(pt[which]);
  auto angle = [](const Vec4D &a, const Vec4D &b) {
    return Vec3D(a)*Vec3D(b)/(Vec3D(a).Abs()*Vec3D(b).Abs()); };
  const int e(angle(pt[l], k) >= angle(pt[r], k) ? l : r), sp(e == l ? r : l);
  const Vec4D Q(pt[e] + pt[sp] + k);
  const double Q2(Q.Abs2()), me(fl[e].Mass()), ms(fl[sp].Mass());
  if (!(Q2 > sqr(me + ms)) || !(Q[0] > 0.)) return false;
  Poincare rest(Q);
  Vec4D q(pt[e] + k);
  rest.Boost(q);
  Vec3D n(q);
  if (!(n.Abs() > 0.)) return false;
  n = n/n.Abs();
  const double pcm(0.5*sqrt(Lambda(Q2, me*me, ms*ms)/Q2));
  Vec4D pe(sqrt(me*me + pcm*pcm), pcm*n), ps(sqrt(ms*ms + pcm*pcm), -pcm*n);
  rest.BoostBack(pe);
  rest.BoostBack(ps);
  out.clear();
  for (size_t i(0); i < pt.size(); ++i) {
    if (i == which) continue;
    if ((int)i == e) out.push_back(pe);
    else if ((int)i == sp) out.push_back(ps);
    else out.push_back(pt[i]);
  }
  return true;
}

/*
  RR_MODE 1: the beams of pt rebuilt (on shell, back to back along their own
  axis in the rest frame of the outgoing momenta) so that pt balances
  exactly. The event's momenta balance to ~1e-9 GeV, which is the energy
  scale of its softest photons (x ~ 1e-7); an imbalance of that size moved
  R_2 of such a pair by factors of 2 against S S B.
*/
static bool BalanceOnBeams(Vec4D_Vector &pt, const ATOOLS::Flavour_Vector &fl)
{
  if (pt.size() < 3) return false;
  Vec4D T;
  for (size_t i(2); i < pt.size(); ++i) T += pt[i];
  const double T2(T.Abs2()), m1(fl[0].Mass()), m2(fl[1].Mass());
  if (!(T2 > sqr(m1 + m2)) || !(T[0] > 0.)) return false;
  Poincare rest(T);
  Vec4D a(pt[0]);
  rest.Boost(a);
  Vec3D n(a);
  if (!(n.Abs() > 0.)) return false;
  n = n/n.Abs();
  const double pcm(0.5*sqrt(Lambda(T2, m1*m1, m2*m2)/T2));
  Vec4D pa(sqrt(m1*m1 + pcm*pcm),  pcm*n), pb(sqrt(m2*m2 + pcm*pcm), -pcm*n);
  rest.BoostBack(pa);
  rest.BoostBack(pb);
  pt[0] = pa;
  pt[1] = pb;
  return true;
}

/*
  RR_MODE 1: the point pt with the initial-state photons in `drop` (indices
  into pt) taken out the way the single real's scaled point does it
  (REAL_MAP 2, MapMomentaScaled): the final system F (final legs and the
  final-state photons, isfsr[slot] for pt[nlegs + slot]) keeps its mass s',
  the kept initial-state photons K and the beams P are scaled by
  x = sqrt(s'/(P - K)^2) and F is boosted rigidly onto x (P - K). Angles to
  the beams and energy fractions of the kept photons, and every invariant
  inside F, are unchanged; with two photons left it is the single real's
  point of the other. Used for the photons other than the pair
  (EventPairPoint): taking a hard one out of its beam instead left that beam
  at a fraction of its energy, and a hard collinear companion inside the
  widened dead cone (4.8e-6 rad, m/E from 1.1e-5 to 5.5e-5) had 400 times
  less crude than in its single real: one pair +171 B in 50k Z-pole nu nu
  events.
*/
static bool ScaleOutISR(const Vec4D_Vector &pt, size_t nlegs, const std::vector<bool> &isfsr,
                        const std::vector<size_t> &drop, const ATOOLS::Flavour_Vector &fl,
                        Vec4D_Vector &out)
{
  if (nlegs < 3 || pt.size() < nlegs || isfsr.size() != pt.size() - nlegs) return false;
  auto dropped = [&](size_t l) {
    return std::find(drop.begin(), drop.end(), l) != drop.end(); };
  for (size_t l : drop) if (l < nlegs || l >= pt.size() || isfsr[l - nlegs]) return false;
  const Vec4D P(pt[0] + pt[1]);
  Vec4D F, K;
  for (size_t l(2); l < pt.size(); ++l) {
    if (dropped(l)) continue;
    if (l < nlegs || isfsr[l - nlegs]) F += pt[l];
    else K += pt[l];
  }
  const Vec4D R(P - K);
  const double sp(F.Abs2()), R2(R.Abs2());
  if (!(sp > 0.) || !(R2 > 0.) || !(R[0] > 0.) || !(F[0] > 0.)) return false;
  const double x(Min(1., sqrt(sp/R2)));
  const double sj(x*x*P.Abs2()), m1(fl[0].Mass()), m2(fl[1].Mass());
  if (!(x > 0.) || !(sj > sqr(m1 + m2))) return false;
  Poincare toP(P);
  Vec4D a(pt[0]);
  toP.Boost(a);
  Vec3D n(a);
  if (!(n.Abs() > 0.)) return false;
  n = n/n.Abs();
  const double pcm(0.5*sqrt(Lambda(sj, m1*m1, m2*m2)/sj));
  Vec4D pa(sqrt(pcm*pcm + m1*m1), pcm*n), pb(sqrt(pcm*pcm + m2*m2), -pcm*n);
  toP.BoostBack(pa);
  toP.BoostBack(pb);
  const Vec4D Fp(pa + pb - x*K);
  if (!(Fp[0] > 0.) || !(Fp.Abs2() > 0.)) return false;
  Poincare fromF(F), toFp(Fp);
  out.assign({pa, pb});
  for (size_t l(2); l < pt.size(); ++l) {
    if (dropped(l)) continue;
    Vec4D q(pt[l]);
    if (l < nlegs || isfsr[l - nlegs]) { fromF.Boost(q); toFp.BoostBack(q); }
    else q = x*q;
    out.push_back(q);
  }
  return BalanceOnBeams(out, fl);
}

/*
  RR_MODE 1: the (n+2)-body point of the pair (i, j) taken from the EVENT:
  the full beams, the post-emission final legs and every photon, with the
  photons other than i and j reduced away: each final-state one recombined
  with its emitter, then the initial-state ones scaled out together
  (ScaleOutISR, the single real's REAL_MAP 2 convention). For a two-photon
  event it is the event itself, so the photons sit at the angles, relative
  to the legs, the generator produced them at (the two-photon map rebuilt
  the pair and put hard collinear photons deep into the dead cone of legs
  they were not generated on). False when the event does not balance.
*/
bool NLO_Base::EventPairPoint(size_t i, size_t j, Vec4D_Vector &P)
{
  const size_t nl(m_flavs.size());
  if (m_postlab.size() != nl || m_bornMomenta.size() < 2) return false;
  Vec4D_Vector ev(m_bornMomenta.begin(), m_bornMomenta.begin() + 2);
  for (size_t l(2); l < nl; ++l) ev.push_back(m_postlab[l]);
  std::vector<size_t> ids;
  for (size_t g(0); g < m_photons.size(); ++g) { ev.push_back(m_photons[g].K()); ids.push_back(g); }
  /*
    The event balances to ~1e-6 GeV (11% of Z-pole mu mu pairs fail
    CheckMomentumConservation's tolerance: p_T up to 4e-6 GeV); the beams are
    rebuilt on their axis in the rest frame of the outgoing momenta, a tilt
    of ~1e-7 rad, far inside every dead cone (m_e/E ~ 1e-5). Larger
    imbalances are not the event.
  */
  {
    Vec4D bal(ev[0] + ev[1]);
    for (size_t l(2); l < ev.size(); ++l) bal -= ev[l];
    if (!(Max(Max(std::abs(bal[0]), std::abs(bal[1])), Max(std::abs(bal[2]), std::abs(bal[3])))
          < 1e-6*sqrt(m_s)) || !BalanceOnBeams(ev, m_flavs)) return false;
  }
  // the other final-state photons recombined with their emitters (last
  // first so that indices stay valid), then the other initial-state photons
  // scaled out together
  for (size_t a(ids.size()); a-- > 0; ) {
    const YFS::Photon &g(m_photons[ids[a]]);
    if (ids[a] == i || ids[a] == j || !g.IsFSR()) continue;
    Vec4D_Vector out;
    bool ok(g.Dip() != nullptr
            && RecombineWithEmitter(ev, nl + a, g.Dip()->Left(), g.Dip()->Right(), m_flavs, out));
    if (!ok) ok = AbsorbPhoton(ev, nl, nl + a, out);
    if (!ok) return false;
    ev = out;
    ids.erase(ids.begin() + a);
  }
  std::vector<bool> isfsr;
  std::vector<size_t> drop;
  for (size_t a(0); a < ids.size(); ++a) {
    isfsr.push_back(m_photons[ids[a]].IsFSR());
    if (ids[a] != i && ids[a] != j) drop.push_back(nl + a);
  }
  if (!drop.empty()) {
    Vec4D_Vector out;
    if (!ScaleOutISR(ev, nl, isfsr, drop, m_flavs, out)) return false;
    ev = out;
    std::vector<size_t> keep;
    for (size_t id : ids) if (id == i || id == j) keep.push_back(id);
    ids = keep;
  }
  if (ev.size() != nl + 2 || ids.size() != 2) return false;
  P = ev;
  if (ids[0] != i) std::swap(P[nl], P[nl + 1]);   // photon i first, j last
  return true;
}

double NLO_Base::CrudeOnLegs(const Vec4D_Vector &pt, const Vec4D &k)
{
  Vec4D_Vector legs(pt.begin(), pt.begin() + m_flavs.size());
  p_nlodipoles->MakeDipolesII(m_flavs, legs, legs);
  p_nlodipoles->MakeDipoles(m_flavs, legs, legs);
  double s(0.);
  if (p_nlodipoles->HasDipoleII()) {
    YFS::Dipole &D(p_nlodipoles->GetDipoleII());
    s += D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
  }
  for (auto &D : p_nlodipoles->GetDipoleFF())
    s += D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
  return s;
}

/*
  RR_MODE 1: the pair (i, j)'s exact remainder, see RRMode. P is the event's
  pair point (EventPairPoint; the two-photon MapMomenta point when the event
  does not balance); R_2 is an RR_Generator tree call.
*/
double NLO_Base::RealRealRemainder(size_t i, size_t j)
{
  if (i >= m_rrphot.size() || j >= m_rrphot.size()) return 0.;
  const RRPhotonInfo &gi(m_rrphot[i]), &gj(m_rrphot[j]);
  if (!gi.ok || !gj.ok || !(gi.denom > 0.) || !(gj.denom > 0.) || m_born == 0.) return 0.;
  Vec4D ki(m_photons[i].K()), kj(m_photons[j].K());
  const double phemin(PhotonEminNLO());
  if (phemin > 0. && (ki.E() < phemin || kj.E() < phemin)) return 0.;
  if (Min(ki.E(), kj.E()) < m_rr_soft_cut*sqrt(m_s)) {
    m_softRR++;
    return 0.;
  }
  const size_t nl(m_flavs.size());
  Vec4D_Vector P;
  const bool evpt(EventPairPoint(i, j, P));
  if (!evpt) {
    // fall back to the two-photon map
    P = m_plab;
    MapMomenta(P, ki, kj);
    P.push_back(ki);
    P.push_back(kj);
    CheckMasses(P, 2);
  }
  if (!CheckMomentumConservation(P) || !BalanceOnBeams(P, m_flavs)) { m_zeroRR++; return 0.; }
  const double r2(p_realreal->Calc_R(P)/(2.*pow(2.*M_PI, 6))*BornPhotonSym(2));
  if (p_realreal->FailCut()) { m_failcut = true; return 0.; }
  if (IsBad(r2)) { m_zeroRR++; return 0.; }
  /*
    Each photon's density c is that of its own single-real point: for an
    initial-state photon of a multi-photon event the scaled point (REAL_MAP
    2) carries the photon scaled up, the density ~ 1/x^2 with it. At P it is
    carried by the ratio of the incoherent crude eikonals, photon and legs of
    P over photon and legs of the point (1 when P is the point). Carrying it
    on the generation legs (full beams, m_plab) instead made the double real
    +350% of sigma_fid (a collinear FSR photon against legs it was not
    emitted from); the density on the physical legs is the one the one-photon
    identity uses (SUB8_FSR_LEGS).
  */
  auto carried = [&](const RRPhotonInfo &g, const Vec4D &kP) {
    if (g.pt.size() != m_flavs.size() + 1) return g.denom;
    const double a(CrudeOnLegs(P, kP)), b(CrudeOnLegs(g.pt, g.pt.back()));
    return (a > 0. && b > 0.) ? g.denom*a/b : g.denom; };
  const double ci(carried(gi, P[nl])), cj(carried(gj, P[nl + 1]));
  /*
    The pair's flux: the singles' (Q - k)^2/Q^2 for an initial-state photon
    (full beams Q, FSR photons 1 with REAL_FSR_FLUX 1) taken with both
    initial-state photons of the pair, (Q - K_ISR)^2/Q^2, the event's s'/s
    for two photons. Other flux conventions keep F_i F_j.
  */
  double Fij(gi.flux*gj.flux);
  {
    static const int fsrflux(ATOOLS::Settings::GetMainSettings()["YFS"]
                             ["REAL_FSR_FLUX"].Get<int>());
    if (m_flux_mode == 0 && fsrflux == 1) {
      const Vec4D Q(P[0] + P[1]);
      Vec4D K;
      if (!m_photons[i].IsFSR()) K += P[nl];
      if (!m_photons[j].IsFSR()) K += P[nl + 1];
      if (Q.Abs2() > 0.) Fij = (Q - K).Abs2()/Q.Abs2();
    }
  }
  static const int a0(ATOOLS::Settings::GetMainSettings()["YFS"]
                      ["REAL_ALPHA0"].SetDefault(1).Get<int>());
  const double kap(m_rescale_alpha);
  const double rhoij((a0 ? kap*kap : 1.)*r2*Fij/(ci*cj*m_born));
  double others(1.);
  for (size_t k(0); k < m_rrphot.size(); ++k)
    if (k != i && k != j) others *= m_rrphot[k].factor;
  const double rem(m_born*(rhoij - gi.factor*gj.factor)*others);
  if (m_photons.size() == 2) m_rrExactPair = rhoij;
  static const int probe(ATOOLS::Settings::GetMainSettings()["YFS"]["RR_PROBE"]
                         .SetDefault(0).Get<int>());
  if (probe) {
    std::ostringstream o;
    o<<std::setprecision(10)<<"@@@ RRPAIR xi="<<2.*m_photons[i].E()/sqrt(m_s)
     <<" xj="<<2.*m_photons[j].E()/sqrt(m_s)<<" fi="<<(m_photons[i].IsFSR()?1:0)
     <<" fj="<<(m_photons[j].IsFSR()?1:0)<<" nph="<<m_photons.size()
     <<" rem/B="<<rem/m_born<<" rhoij="<<rhoij<<" facti="<<gi.factor
     <<" factj="<<gj.factor<<" others="<<others<<" ci/c="<<ci/gi.denom
     <<" cj/c="<<cj/gj.denom<<" Fij/FiFj="<<Fij/(gi.flux*gj.flux)
     <<" evpt="<<(evpt?1:0);
    if (probe >= 2) {
      o<<std::setprecision(17)<<" P=";
      for (const Vec4D &q : P) o<<q<<" ";
      o<<" PB=";
      for (const Vec4D &q : m_plab) o<<q<<" ";
    }
    o<<"\n";
    std::cerr<<o.str();
  }
  if (IsBad(rem)) { m_zeroRR++; return 0.; }
  if (!IsZero(rem)) m_nonZeroRR++;
  return rem;
}

double NLO_Base::CalculateRealReal() {
  m_rrExactPair = std::numeric_limits<double>::quiet_NaN();
  if (!m_rrtool)
    return 0;
  if (RRMode() == 1) {
    static bool checked(false);
    if (!checked) {
      checked = true;
      if (!(m_rr_soft_cut > 0.))
        THROW(fatal_error, "YFS: RR_MODE 1 needs RR_SOFT_CUT > 0 (the lower edge of "
              "the soft x hard form-factor difference in the real-virtual).");
      if (!m_realvirt || RVMode() != 1)
        msg_Error()<<"YFS: RR_MODE 1 without RV_MODE 1: the double real's soft x "
                   <<"hard log is not cancelled; the result depends on RR_SOFT_CUT."
                   <<std::endl;
    }
  }
  double rr(0);
  m_rr_hard2 = 0.;
  const YFS::Photon_Vector &photons(m_photons);
  if (photons.size() == 0)
    return 0;
  const YFS::Photon *h1 = nullptr, *h2 = nullptr;
  for (const YFS::Photon &g : photons) {
    if (!h1 || g.E() > h1->E())      { h2 = h1; h1 = &g; }
    else if (!h2 || g.E() > h2->E()) { h2 = &g; }
  }
  if (RRMode() == 1) {
    for (size_t i(0); i < photons.size(); ++i)
      for (size_t j(i + 1); j < photons.size(); ++j) {
        const double contrib(RealRealRemainder(i, j));
        rr += contrib;
        const YFS::Photon *pi(&photons[i]), *pj(&photons[j]);
        if ((pi == h1 && pj == h2) || (pi == h2 && pj == h1)) m_rr_hard2 = contrib;
      }
    return rr;
  }
  for (int i = 0; i < photons.size(); ++i) {
    for (int j = i + 1; j < photons.size(); ++j) {
      Vec4D k = photons[i].K();
      Vec4D kk = photons[j].K();
      const double phemin(PhotonEminNLO());
      if (phemin>0.0 && (k.E()<phemin || kk.E()<phemin)) continue;
      double contrib = CalculateRealReal(k, kk);
      static const bool betacheck(ATOOLS::Settings::GetMainSettings()["YFS"]["BETA_RECURSION"].Get<int>()!=0);
      if (betacheck) {
        const double gen(CalculateRealN((1u<<i) | (1u<<j)));
        std::cerr<<"@@@ BETA2 hand="<<contrib<<" rec="<<gen
                 <<" dev="<<(contrib!=0. ? std::abs(gen/contrib-1.) : -1.)
                 <<std::endl;
      }
      rr += contrib;
      const YFS::Photon *pi = &photons[i], *pj = &photons[j];
      if ((pi == h1 && pj == h2) || (pi == h2 && pj == h1))
        m_rr_hard2 = contrib;
      if (m_check_rr_sub == 2) {
        // accumulating scatter: record each photon of the pair with the pair
        // residual, so energetic collinear photons in large residuals show up
        RecordSubScatter(k,  contrib, "rr", m_rr_eik);
        RecordSubScatter(kk, contrib, "rr", m_rr_eik);
      }
      if (m_check_rr_sub == 1) {
        // k*=2;
        // kk*=2;
        if (k.E() < 0.2 * sqrt(m_s))
          continue;
        if (kk.E() < 0.2 * sqrt(m_s))
          continue;
        if (!m_failcut)
          CheckRealRealSub(k, kk);
      }
    }
  }
  return rr;
}

double NLO_Base::CalculateRealReal(Vec4D k1, Vec4D k2) {
  if (m_rr_soft_cut > 0.) {
    const double emin = Min(k1.E(), k2.E());
    if (emin < m_rr_soft_cut * sqrt(m_s)) {
      m_softRR++;
      msg_Debugging() << METHOD << " skipping soft photon pair: min(E1,E2)="
                      << emin << " < RR_SOFT_CUT*sqrt(s)="
                      << m_rr_soft_cut * sqrt(m_s) << "\n";
      return 0;
    }
  }
  const double norm = 2. * pow(2 * M_PI, 6);
  m_rr_eik = 0.;  // reset; set once the crude eikonals are computed below
  Vec4D_Vector p(m_plab);
  Vec4D_Vector pp = p;
  Vec4D kk1 = k1, kk2 = k2;

  msg_Debugging() << METHOD
                  << " k1=" << k1 << " E1=" << k1.E() << " pt1=" << k1.PPerp()
                  << " k2=" << k2 << " E2=" << k2.E() << " pt2=" << k2.PPerp() << "\n";

  MapMomenta(p, k1, k2);

  p.push_back(k1);
  p.push_back(k2);

  Vec4D_Vector _p = p;
  _p.pop_back();
  _p.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, _p, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, _p, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, _p, m_plab);

  const double subloc1 = p_nlodipoles->CalculateRealSub(k1);
  const double subloc2 = p_nlodipoles->CalculateRealSub(k2);

  double flux;
  if (m_flux_mode == 1)
    flux = p_nlodipoles->CalculateFlux(k1 + k2);
  else
    flux = p_dipoles->CalculateFlux(k1) * p_dipoles->CalculateFlux(k2);

  msg_Debugging() << METHOD << " subloc1=" << subloc1 << " subloc2=" << subloc2
                  << " flux=" << flux << " (mode=" << m_flux_mode << ")\n";

  if (!CheckMomentumConservation(p)) {
    msg_Debugging() << METHOD << " momentum conservation failed, returning 0\n";
    m_zeroRR++;
    return 0;
  }

  double r = p_realreal->Calc_R(p) / norm * BornPhotonSym(2);
  if (p_realreal->FailCut()) {
    msg_Debugging() << METHOD << " FailCut triggered, returning 0\n";
    m_failcut = true;
    m_zeroRR++;
    return 0;
  }
  if (IsBad(r) || IsBad(flux)) {
    msg_Debugging() << METHOD << " bad point: r=" << r << " flux=" << flux << "\n";
    m_zeroRR++;
    return 0;
  }

  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);


  const double sub1  = p_dipoles->CalculateRealSubEEX(kk1);
  const double sub2  = p_dipoles->CalculateRealSubEEX(kk2);
  m_rr_eik = sub1 * sub2;
  const double real1 = CalculateReal(kk1, /*raw*/true);
  const double F1(m_lastflux), SL1(m_lastsubloc), SB1(m_lastdenom);
  const double SL01(m_lastsubloc0), SB01(m_lastdenom0);
  // the eikonal of k2 on the legs of photon 1's one-photon point (the
  // dipoles CalculateReal just built there), and vice versa below
  const double A2raw(p_nlodipoles->CalculateRealSub(kk2));
  const double real2 = CalculateReal(kk2, /*raw*/true);
  const double F2(m_lastflux), SL2(m_lastsubloc), SB2(m_lastdenom);
  const double SL02(m_lastsubloc0), SB02(m_lastdenom0);
  const double A1raw(p_nlodipoles->CalculateRealSub(kk1));
  // what REAL_SUB_EIK / REAL_FSR_FLUX did to each photon's eikonal and
  // denominator in the single real, as a factor (exactly 1 in mode 0)
  auto ratio = [](double a, double b) { return (b != 0. && !IsBad(a/b)) ? a/b : 1.; };
  const double rs1(ratio(SL1, SL01)), rs2(ratio(SL2, SL02));
  const double rd1(ratio(SB1, SB01)), rd2(ratio(SB2, SB02));
  m_recola_evts += 1;
  /*
    YFS: RR_CONVENTIONS (1). beta_2 = R_2 F_1 F_2 - S_2 beta_1(k_1)
    - S_1 beta_1(k_2) - S_1 S_2 B over the crude of each photon, with the flux,
    subtraction eikonal and crude of each photon taken from the single real's
    own calculation of that photon (m_last*). The old assembly (0) built its
    own: the plain eikonal at the double-mapped point, an initial-state flux
    product for every photon and the EEX crude, while beta_1 inside it came
    from the single real. With REAL_SUB_EIK 8 and REAL_FSR_FLUX 1 the two no
    longer cancelled for hard photons: Z-pole mu mu YFS.NLO+RR -187 +- 175%
    and 404 +- 57% pb against 1045 +- 5.7% with the old single-real defaults
    (2026-09-26). The coupling powers are also per term: an eikonal is
    alpha(0), the reals the model's, so S beta_1 carries one m_rescale_alpha
    and S S B two (the old code divided all of it by one).
  */
  /*
    Status 2026-09-26 (Z-pole mu mu, 20k events, Check_RR_Sub 1 for the soft
    limits): 0 is the old assembly. It is small (+0.5%) and passes both
    single-soft limits only with the old single-real settings AND
    REAL_FSR_MAP 0; with REAL_FSR_MAP 2 or 4 (the two-photon FSR map does not
    reduce to the one-photon one) |beta_2|/B stalls at 0.7 as one photon goes
    soft, and the correction is -11%. 3 (default) passes every soft limit with
    REAL_FSR_MAP 0 under REAL_SUB_EIK 8 / REAL_FSR_FLUX 1, but is still -2.8%;
    7 adds the denominator factor (+3.6 +- 3.5%); 8 (eikonals on the other
    photon's legs) made it worse. Open: the FSR maps, and REAL_SUB_EIK 8's
    multichannel crude inside beta_2.
  */
  static const int rrconv(ATOOLS::Settings::GetMainSettings()["YFS"]
                          ["RR_CONVENTIONS"].SetDefault(3).Get<int>());
  if (rrconv) {
    // bitmask while the right combination is established: 1 = the single
    // real's flux, 2 = its subtraction eikonal, 4 = its denominator; a clear
    // bit keeps the old double-real piece (flux product, eikonal at the
    // double-mapped point, EEX crude)
    // 1: the single real's flux; 2: the double real's own eikonal (at its
    // double-mapped point) times the single real's modification factor;
    // 4: the same for the denominator (the EEX crude of the event)
    const double fF(rrconv & 1 ? F1 * F2 : flux);
    const double s1(rrconv & 2 ? subloc1 * rs1 : subloc1), s2(rrconv & 2 ? subloc2 * rs2 : subloc2);
    const double d1(rrconv & 4 ? sub1 * rd1 : sub1), d2(rrconv & 4 ? sub2 * rd2 : sub2);
    if (IsZero(real1) || IsZero(real2) || !(d1 > 0.) || !(d2 > 0.)) {
      m_zeroRR++;
      return 0;
    }
    const double ra(m_rescale_alpha);
    double num(r * fF - (s2 * real1 + s1 * real2) / ra
               - s1 * s2 * m_born / (ra * ra));
    /*
      8: the soft limits. As k2 -> soft, R_2 -> S(k2; legs of photon 1's
      point) R_1(k1), so the eikonal multiplying beta_1(k1) must be k2's on
      THOSE legs (a2), not at the double-mapped point where k1 is removed;
      for FSR the final legs recoil against a hard k1 and the two differ at
      O(1), leaving a finite residual for every soft photon (Z-pole mu mu
      YFS.NLO+RR -11%, nu nu with initial-state radiation only +0.005%).
      beta_2 vanishes in both single-soft limits with
      c = a2 S_1 + a1 S_2 - S_1 S_2, S_i the single real's subtraction.
    */
    if (rrconv & 8) {
      const double a1(A1raw * rs1), a2(A2raw * rs2);
      num = r * fF - (a2 * real1 + a1 * real2) / ra
            - (a2 * SL1 + a1 * SL2 - SL1 * SL2) * m_born / (ra * ra);
    }
    const double tot1(num / (d1 * d2));
    m_rr_eik = d1 * d2;
    if (IsBad(tot1)) {
      msg_Error() << METHOD << " NNLO RR is NaN: r=" << r << " F1=" << F1 << " F2=" << F2
                  << " SL1=" << SL1 << " SL2=" << SL2 << " SB1=" << SB1 << " SB2=" << SB2 << "\n";
      return 0;
    }
    if (!IsZero(tot1)) m_nonZeroRR++;
    return tot1;
  }

  msg_Debugging() << METHOD << " r=" << r
                  << " sub1=" << sub1 << " sub2=" << sub2
                  << " real1=" << real1 << " real2=" << real2
                  << " born=" << m_born << "\n";

  if (IsZero(real1) || IsZero(real2)) {
    msg_Debugging() << METHOD << " real1 or real2 is zero, returning 0\n";
    m_zeroRR++;
    return 0;
  }

  const double fullsub = -subloc2 * real1 - subloc1 * real2 - subloc1 * subloc2 * m_born;
  const double tot     = (r * flux + fullsub / m_rescale_alpha) / sub1 / sub2;

  msg_Debugging() << METHOD << " fullsub=" << fullsub
                  << " r*flux=" << r * flux
                  << " tot=" << tot << "\n";

  if (IsBad(tot))
    msg_Error() << METHOD << " NNLO RR is NaN: r=" << r << " flux=" << flux
                << " fullsub=" << fullsub << " sub1=" << sub1 << " sub2=" << sub2 << "\n";

  if (!IsZero(tot)) m_nonZeroRR++;
  return tot;
}

// ======================================================================
//  General-n real corrections: one subset recursion for every multiplicity
// ======================================================================

YFS::Real_Correction *NLO_Base::RealProvider(size_t nphotons) const
{
  if (nphotons == 1) return p_real;
  if (nphotons == 2) return p_realreal;
  if (nphotons < m_realprov.size()) return m_realprov[nphotons];
  return NULL;
}

void NLO_Base::SetRealProvider(size_t nphotons, YFS::Real_Correction *prov)
{
  if (m_realprov.size() <= nphotons) m_realprov.resize(nphotons+1, NULL);
  m_realprov[nphotons] = prov;
}

size_t NLO_Base::MaxRealPhotons() const
{
  size_t n(0);
  if (p_real)     n = 1;
  if (p_realreal) n = 2;
  for (size_t i(m_realprov.size()); i-- > 3; )
    if (m_realprov[i] != NULL) { n = Max(n, i); break; }
  return n;
}

size_t NLO_Base::RequestedMaxRealPhotons()
{
  static const int n
    (ATOOLS::Settings::GetMainSettings()["YFS"]["NLO_MAX_PHOTONS"]
     .SetDefault(2).Get<int>());
  return n < 1 ? 1 : (size_t)n;
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
  if (m_flux_mode == 1) {
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

double NLO_Base::CalculateVV() {
  if (!m_vvtool)
    return 0;
  if (m_eex_virt) {
    return p_dipoles->CalculateEEXVirtual() * m_born - m_born;
  }
  double virt;
  double sub;
  // CheckMassReg();
  if (!HasISR())
    virt = p_vv->Calc(m_bornMomenta, m_born);
  else
    virt = p_vv->Calc(m_plab, m_born);
  if (m_check_virt_born) {
    // the provider's Born is pointlike, m_born is dressed with the pion form
    // factor, so compare against the dressed provider Born
    if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born()
                         * ExternalFormFactor(m_plab, m_flavs), 1e-6)) {
      msg_Error() << METHOD
                  << "\n Warning! Loop provider's born is different! YFS "
                     "Subtraction likely fails\n"
                  << "Loop Provider " << ":  " << p_virt->p_loop_me->ME_Born()
                  << "\nSherpa" << ":  " << m_born << std::endl
                  << "PhaseSpace Point = ";
      for (auto _p : m_plab)
        msg_Error() << _p << std::endl;
    }
  }
  if (p_vv->FailCut())
    return 0;
  if (m_virt_sub && p_virt->p_loop_me->Mode() != 1)
    sub = p_dipoles->CalculateVirtualSub();
  else
    sub = 0;
  double sub2 = p_dipoles->CalculateVVSubEps();
  // m_oneloop = (virt- sub * m_born/m_rescale_alpha );
  m_oneloop = (virt - sub * CalculateVirtual() / m_rescale_alpha -
               0.5 * sub * sub * m_born / m_rescale_alpha);
  if (p_virt->p_loop_me->Mode() == 1) {
    m_oneloop /= m_rescale_alpha;
  }
  if (IsBad(m_oneloop) || IsBad(sub)) {
    msg_Error() << "YFS Virtual is NaN" << std::endl
                << "Virtual:  " << virt << std::endl
                << "Subtraction: " << sub * m_born << std::endl
                << "PhaseSpace Point: " << std::endl
                << m_plab << std::endl;
  }
  double loope1 =
      p_vv->p_loop_me->ME_E1() *
      p_vv->m_factor; //*p_vv->m_factor;//+p_virt->p_loop_me->ME_E1()*p_virt->m_factor;;
  double loope2 =
      2. * p_vv->p_loop_me->ME_E2() * p_vv->m_factor * p_vv->m_factor;
  double yfse1 = p_dipoles->Get_E1();
  double yfse2 = p_dipoles->GetVV_E2();
  if (m_check_poles == 1) {
    if (m_virt_sub == 0)
      sub = p_dipoles->CalculateVirtualSub();
    const double p1 = p_vv->p_loop_me->ME_E1() * p_vv->m_factor;
    const double p2 =
        2. * p_vv->p_loop_me->ME_E2() * p_vv->m_factor * p_vv->m_factor;
    const double yfspole1 = (p_dipoles->Get_E1());
    const double yfspole2 = p_dipoles->GetVV_E2();
    PRINT_VAR(p1 / yfspole1);
    int ncorrect1 = ::countMatchingDigits(p1, yfspole1, 32);
    int ncorrect2 = ::countMatchingDigits(p2, -yfspole2, 32);
    if (!IsEqual(p2, -yfspole2, 1e-6) || ncorrect1 < 10) {
      msg_Error() << "Poles do not cancel in YFS Double Virtuals" << std::endl
                  << "Correct digits \epsion^{-1} =  " << ncorrect1 << std::endl
                  << "Correct digits \epsion^{-2} =  " << ncorrect2
                  << std::endl;
      return 0;
    } else {
      int i = 0;
      msg_Debugging() << std::setprecision(32);
      msg_Out() << "Poles cancel in YFS double Virtuals to " << ncorrect2
                << " digits" << std::endl;
      m_histograms1d["SinglePoleVV"]->Insert(ncorrect1);
      m_histograms1d["DoublePoleVV"]->Insert(ncorrect2);
    }
  }
  return 0;
}

void NLO_Base::RandomRotate(Vec4D &p) {
  Vec4D t1 = p;
  // rotate around x
  p[2] = cos(m_ranTheta) * t1[2] - sin(m_ranTheta) * t1[3];
  p[3] = sin(m_ranTheta) * t1[2] + cos(m_ranTheta) * t1[3];
  t1 = p;
  // rotate around z
  p[1] = cos(m_ranPhi) * t1[1] - sin(m_ranPhi) * t1[2];
  p[2] = sin(m_ranPhi) * t1[1] + cos(m_ranPhi) * t1[2];
}

void NLO_Base::CheckMappingRecoil(const Vec4D_Vector &p, const Vec4D &ksum) {
  Vec4D q;
  for (size_t i(2); i < p.size(); ++i) q += p[i];
  const Vec3D imbalance(Vec3D(ksum) + Vec3D(q));
  const double scale(Max(Vec3D(ksum).Abs(), Vec3D(q).Abs()));
  if (scale <= 0.) return;
  const double rel(imbalance.Abs() / scale);
  if (rel > 1e-5)
    msg_Error() << METHOD << "(): YFS mapping recoil imbalance " << rel
                << " (photons " << Vec3D(ksum) << ", hard system " << Vec3D(q)
                << ")" << std::endl;
}

/*
  The (n+1)-body point for beta_1(k_j) of FINAL-state photons (REAL_FSR_MAP: 1).

  p[2..] are the final legs BEFORE final-state emission (m_reallab), so their
  sum is the ISR-reduced total momentum Q, of mass sqrt(s'). The rest-frame
  map below adds k on top and rebuilds the beams at sqrt((Q+k)^2), i.e. ABOVE
  sqrt(s'): e+e- -> mu mu at 0.7 GeV gave sqrt(s_j) = 0.83-0.89 GeV for hard
  wide-angle FSR photons, the real ME and flux were taken at that unphysical
  point while the subtraction used S~ B at s', and the fixed-order photon
  spectrum grew a heavy-weight bump above E_gamma ~ 0.12 GeV (independent of
  the ME provider and of REAL_MAP, which touches ISR photons only).

  Here Q is kept: in the Q rest frame the final legs keep their directions and
  are rescaled by one factor xi so that they carry the mass sqrt((Q-K)^2),
  then boosted to move with Q-K; the beams are rebuilt at sqrt(Q^2) = sqrt(s').
  With a single FSR photon this is the YFS final-state recoil of the event.
*/
bool NLO_Base::MapMomentaFSR(Vec4D_Vector &p, Vec4D_Vector &k) {
  if (k.empty() || p.size() < 4) return false;
  for (const Vec4D &kj : k) if (!PhotonIsFSR(kj)) return false;
  m_map_reduced = false;
  /*
    Start from the legs AFTER final-state emission, not before. With the
    pre-emission legs the photon-lepton angle of the point is not the event's,
    and the collinear structure of |M_1|^2 is off by large factors: on e+e- ->
    mu mu at 0.7 GeV, single-FSR-photon events had r/rho_1(CEEX) spread over
    0.003-0.5 where single-ISR-photon events give 0.01613 on every event. The
    other FSR photons are absorbed into the pair below (their sum is added to
    the target mass); with one FSR photon the point IS the event, as for CEEX.
  */
  Vec4D Kall;
  for (const Vec4D &kj : m_FSRPhotons) Kall += kj;
  /*
    REAL_FSR_MAP: 3 - the SINGLE-EMISSION point: the pre-emission legs with
    photon j's own recoil applied by the rescaling below and nothing else.
    The post-emission construction above re-absorbs the other photons by
    rescaling, which keeps the post-emission DIRECTIONS; a soft photon whose
    parent lepton was kicked by a hard companion then sits at a point whose
    collinear structure is that of the post-emission legs while the crude it
    is subtracted against (REAL_SUB_EIK 3, the generation density) is the
    pre-emission one, and its bracket (r_j - S~ B)/S~ is O(1) instead of
    vanishing: Z-pole mu mu, nfsr >= 3 all soft, 10th percentile of the
    fixed-order factor -5. Starting from the pre-emission legs the point
    tends to the crude configuration as k_j -> 0, and with one photon it is
    the generator's own post-emission point (same rescaling recipe).
  */
  { static const int fsrmap(ATOOLS::Settings::GetMainSettings()["YFS"]
                            ["REAL_FSR_MAP"].Get<int>());
    if (fsrmap == 3 && m_plab.size() == p.size()) {
      for (size_t i = 2; i < p.size(); ++i) p[i] = m_plab[i];
      Kall = Vec4D();
    } else if (m_postlab.size() == p.size()) {
      for (size_t i = 2; i < p.size(); ++i) p[i] = m_postlab[i];
    } else Kall = Vec4D(); }
  Vec4D Q;
  for (size_t i = 2; i < p.size(); ++i) Q += p[i];
  Q += Kall;
  const double sq(Q.Abs2());
  if (!(sq > 0.)) return false;
  Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
  Poincare boostQ(Q);
  for (size_t i = 0; i < p.size(); ++i) { pRot.RotateBack(p[i]); boostQ.Boost(p[i]); }
  Vec4D K;
  for (Vec4D &kj : k) { pRot.RotateBack(kj); boostQ.Boost(kj); K += kj; }
  const double M(sqrt(sq));
  const Vec4D Qp(Vec4D(M, 0., 0., 0.) - K);
  const double Mp2(Qp.Abs2());
  if (!(Mp2 > 0.)) return false;
  const double Mp(sqrt(Mp2));
  // Into the rest frame of the legs' own sum L, where their 3-momenta cancel
  // and a common rescaling keeps them cancelling; the rescaled legs are then
  // boosted onto Q - K. (Boosting into the Q - K frame instead left the
  // absorbed photons' 3-momentum unbalanced whenever K != all FSR photons.)
  Vec4D L;
  for (size_t i = 2; i < p.size(); ++i) L += p[i];
  if (!(L.Abs2() > 0.)) return false;
  Poincare boostL(L), boostQp(Qp);
  for (size_t i = 2; i < p.size(); ++i) boostL.Boost(p[i]);
  // xi: sum_i sqrt(m_i^2 + xi^2 |p_i|^2) = Mp, by bisection (monotone in xi)
  std::vector<double> m2, pp2;
  double msum(0.);
  for (size_t i = 2; i < p.size(); ++i) {
    const double mi(m_flavs[i].Mass());
    m2.push_back(mi*mi); pp2.push_back(Vec3D(p[i]).Sqr()); msum += mi;
  }
  if (msum >= Mp) return false;
  auto etot = [&](double xi) { double e(0.);
    for (size_t j = 0; j < m2.size(); ++j) e += sqrt(m2[j] + xi*xi*pp2[j]);
    return e; };
  double lo(0.), hi(1.);
  while (etot(hi) < Mp) hi *= 2.;
  for (int it = 0; it < 200; ++it) {
    const double mid(0.5*(lo+hi));
    (etot(mid) < Mp ? lo : hi) = mid;
  }
  const double xi(0.5*(lo+hi));
  for (size_t i = 2; i < p.size(); ++i) {
    const Vec3D v(xi*Vec3D(p[i]));
    Vec4D f(sqrt(m2[i-2] + v.Sqr()), v);
    boostQp.BoostBack(f);
    p[i] = f;
  }
  const double sign_z = (m_bornMomenta[0][3] < 0 ? -1 : 1);
  const double m1 = m_flavs[0].Mass(), m2b = m_flavs[1].Mass();
  const double lamCM = 0.5 * sqrt(Lambda(sq, m1*m1, m2b*m2b) / sq);
  p[0] = {sqrt(lamCM*lamCM + m1*m1), 0, 0,  sign_z * lamCM};
  p[1] = {sqrt(lamCM*lamCM + m2b*m2b), 0, 0, -sign_z * lamCM};
  Poincare pRot2(m_bornMomenta[0], Vec4D(0., 0., 0, 1.));
  for (size_t i = 0; i < p.size(); ++i) { pRot2.Rotate(p[i]); boostLab.BoostBack(p[i]); }
  for (Vec4D &kj : k) { pRot2.Rotate(kj); boostLab.BoostBack(kj); }
  { Vec4D res(p[0] + p[1]);
    for (size_t i = 2; i < p.size(); ++i) res -= p[i];
    for (const Vec4D &kj : k) res -= kj;
    if (Vec3D(res).Abs() > 1e-9 || std::abs(res[0]) > 1e-9) {
      Vec4D post; for (size_t i = 2; i < m_postlab.size(); ++i) post += m_postlab[i];
      Vec4D pre;  for (size_t i = 2; i < m_plab.size(); ++i)    pre  += m_plab[i];
      Vec4D isr;  for (const Vec4D &kj : m_ISRPhotons) isr += kj;
      std::cerr<<std::setprecision(10)<<"@@@ FSRMAP residual="<<res
               <<" xi="<<xi<<" Mp="<<Mp<<" msum="<<msum
               <<" nFSR="<<m_FSRPhotons.size()<<" nk="<<k.size()
               <<"\n   Kall(lab)="<<Kall<<" post(lab)="<<post<<" pre(lab)="<<pre
               <<" isr(lab)="<<isr<<" P="<<(m_bornMomenta[0]+m_bornMomenta[1])
               <<"\n   post+Kall+isr-P="<<(post+Kall+isr-m_bornMomenta[0]-m_bornMomenta[1])
               <<" pre+isr-P="<<(pre+isr-m_bornMomenta[0]-m_bornMomenta[1])<<std::endl;
    } }
  return true;
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
  static const int on(ATOOLS::Settings::GetMainSettings()["YFS"]
                      ["REAL_BORN_PHOTON_SYM"].SetDefault(1).Get<int>());
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
  Born-photon multichannel (YFS: REAL_BORN_PHOTON_MULTICHANNEL).

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
  static const int mode(ATOOLS::Settings::GetMainSettings()["YFS"]
                        ["REAL_BORN_PHOTON_MULTICHANNEL"].SetDefault(1).Get<int>());
  if (mode == 0 || !p_dipoles || !p_dipoles->HasDipoleII()) return 1.;
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
  return mode == 1 ? 1./G : 1.;   // mode 2: compute and report only
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

/*
  The (n+1)-body point for beta_1(k_j) of FINAL-state photons when the final
  state has MORE THAN ONE radiating dipole (REAL_FSR_MAP: 2, the default).

  Each final-state photon is radiated by one resonant pair (its dipole,
  YFS::Photon::Dip()): Dipole::GenerateEmissions samples it in that pair's
  rest frame and YFS_Handler::CalculateFSR lets only that pair recoil, so in
  the event every pair separately satisfies  Q_D(pre) = q_1' + q_2' + K_D.
  MapMomentaFSR above absorbs the OTHER photons by rescaling ALL final legs
  together: with two pairs the pair that did not radiate k_j is moved off its
  pre-emission momentum and k_j's pair does not land on Q_D - k_j, so both
  resonance propagators of |M_1|^2 are evaluated off the event's invariants
  while the subtraction S~ B sits on the event's ones. On the Z pole a shift
  of a fraction of Gamma_Z is a large factor, which is what the numbers show.
  e+e- -> mu mu tau tau at 250 GeV (doc/examples/YFS/zpole/br_table/
  250_eemumu, 8k events, mu-pair mass within 2 GeV of M_Z):
      nfsr = 1: <BR>/<CEEX> = 0.987, median BR/CEEX = 1.06   (the event)
      nfsr = 2: <BR>/<CEEX> = 0.21,  median 1.03 - the mean is a tail;
      nfsr = 3: 0.17;   nfsr >= 4: 0.12;
      a hard FSR photon (x > 0.05) with other FSR photons: median BR = -0.66,
      <BR> = -3.5 against <CEEX> = 0.75.
  In the Rivet Z1_mass overlay the Born+real dipped 25% at the mu-pair peak
  where LO, EEX and CEEX agree in shape.

  Here every pair other than k_j's is put back at its pre-emission momentum
  (m_plab, its photons re-absorbed), and k_j's pair is rebuilt at
  Q_D(pre) - k_j: its two legs keep the direction they have in the event in
  their own rest frame (the post-emission legs of m_postlab boosted to rest,
  which is the direction the generator drew), are given the momentum of a
  two-body decay of mass sqrt((Q_D - k_j)^2), and are boosted onto Q_D - k_j.
  Several selected photons (CalculateRealReal) are grouped by dipole and
  absorbed by their own pair each. Legs in no radiating dipole are untouched.
  When k_j is the only photon of its pair the pair is the event's, and with
  one radiating pair in the process the whole construction is MapMomentaFSR.
  The beams are rebuilt at sqrt(Q^2) = sqrt(s') exactly as there.
*/
/*
  The pole scheme's (n+1)-body point for photons the W-PAIR dipole radiated.
  Such a photon's dipole names the two charged leptons (Left/Right are the
  lepton positions, DipoleSet::BuildPole), but the pair that recoiled is the
  W's, so the flat map above this would rebuild the lepton pair, the wrong
  pre-emission point by the W recoil. Here: the W pair before the selected
  photons, Qp = W- + W+ + K, the two W's back to back in Qp's frame along
  their post-emission direction with their (preserved) masses, and each W's
  daughters carried from the post-emission W to the rebuilt one - the
  inverse of Define_Dipoles::ApplyPoleRecoil. Photons from any other dipole
  (none exist in the pole scheme) fall back to the flat map.
*/
bool NLO_Base::MapMomentaFSRPole(Vec4D_Vector &p, Vec4D_Vector &k) {
  if (p_dipoles == nullptr || !p_dipoles->PoleActive()) return false;
  const YFS::DipoleSet::WWLegs &w(p_dipoles->WW());
  if (!w.ok || w.lm >= p.size() || w.lp >= p.size() || w.nm >= p.size()
      || w.np >= p.size()) return false;
  Vec4D K;
  for (const Vec4D &kj : k) {
    const YFS::Photon *g(FindPhoton(kj));
    if (g == nullptr || !g->IsFSR() || g->Dip() == nullptr) return false;
    if (!g->Dip()->GetFlav(0).IsVector()) return false;   // not the W pair
    K += kj;
  }
  m_map_reduced = false;
  for (size_t i = 2; i < p.size(); ++i) p[i] = m_plab[i];
  const std::size_t dau[2][2] = {{w.lm, w.nm}, {w.lp, w.np}};
  Vec4D wpost[2], wpre[2];
  for (int i(0); i < 2; ++i) wpost[i] = m_postlab[dau[i][0]] + m_postlab[dau[i][1]];
  const Vec4D Qp(wpost[0] + wpost[1] + K);
  const double Mp2(Qp.Abs2());
  if (!(Mp2 > 0.) || !(Qp[0] > 0.)) return false;
  const double m1(wpost[0].Mass()), m2(wpost[1].Mass());
  if (!(m1 > 0.) || !(m2 > 0.) || m1 + m2 >= sqrt(Mp2)) return false;
  // the W direction: the post-emission W's in their own rest frame
  Vec4D q1(wpost[0]), q2(wpost[1]);
  const Vec4D L(q1 + q2);
  if (!(L.Abs2() > 0.) || !(L[0] > 0.)) return false;
  Poincare boostL(L);
  boostL.Boost(q1); boostL.Boost(q2);
  Vec3D n(Vec3D(q1) - Vec3D(q2));
  if (!(n.Abs() > 0.)) return false;
  n = n/n.Abs();
  const double pcm(0.5*sqrt(Lambda(Mp2, m1*m1, m2*m2)/Mp2));
  wpre[0] = Vec4D(sqrt(m1*m1 + pcm*pcm),  pcm*n);
  wpre[1] = Vec4D(sqrt(m2*m2 + pcm*pcm), -pcm*n);
  Poincare boostQp(Qp);
  boostQp.BoostBack(wpre[0]); boostQp.BoostBack(wpre[1]);
  for (int i(0); i < 2; ++i) {
    Poincare toRest(wpost[i]), fromPre(wpre[i]);
    for (int j(0); j < 2; ++j) {
      Vec4D q(m_postlab[dau[i][j]]);
      toRest.Boost(q);
      fromPre.BoostBack(q);
      p[dau[i][j]] = q;
    }
  }
  return RebuildBeamsAtPreFSR(p, k, "FSRMAPP");
}

bool NLO_Base::MapMomentaFSRDipole(Vec4D_Vector &p, Vec4D_Vector &k) {
  if (k.empty() || p.size() < 4) return false;
  if (m_postlab.size() != p.size() || m_plab.size() != p.size()) return false;
  if (p_dipoles && p_dipoles->PoleActive()) return MapMomentaFSRPole(p, k);
  // The selected photons, grouped by the dipole that radiated them.
  std::map<std::pair<int,int>, Vec4D> ksel;
  for (const Vec4D &kj : k) {
    const YFS::Photon *g(FindPhoton(kj));
    if (g == nullptr || !g->IsFSR() || g->Dip() == nullptr) return false;
    const int l(g->Dip()->Left()), r(g->Dip()->Right());
    if (l < 2 || r < 2 || l >= (int)p.size() || r >= (int)p.size() || l == r)
      return false;
    ksel[std::make_pair(std::min(l,r), std::max(l,r))] += kj;
  }
  m_map_reduced = false;
  // Every final leg at its pre-emission momentum: the pairs that did not
  // radiate a selected photon have their own photons re-absorbed by this.
  for (size_t i = 2; i < p.size(); ++i) p[i] = m_plab[i];
  for (const auto &sel : ksel) {
    const int l(sel.first.first), r(sel.first.second);
    const Vec4D K(sel.second);
    const Vec4D Qd(m_plab[l] + m_plab[r]);     // the pair before it radiated
    const Vec4D Qp(Qd - K);                    // the pair with only K taken out
    const double Mp2(Qp.Abs2());
    if (!(Mp2 > 0.) || !(Qp[0] > 0.)) return false;
    const double Mp(sqrt(Mp2));
    const double m1(m_flavs[l].Mass()), m2(m_flavs[r].Mass());
    if (m1 + m2 >= Mp) return false;
    // The decay direction the event has: the post-emission legs in their
    // own rest frame are back to back along it.
    Vec4D q1(m_postlab[l]), q2(m_postlab[r]);
    const Vec4D L(q1 + q2);
    if (!(L.Abs2() > 0.) || !(L[0] > 0.)) return false;
    Poincare boostL(L);
    boostL.Boost(q1); boostL.Boost(q2);
    Vec3D n(Vec3D(q1) - Vec3D(q2));
    if (!(n.Abs() > 0.)) return false;
    n = n/n.Abs();
    const double pcm(0.5*sqrt(Lambda(Mp2, m1*m1, m2*m2)/Mp2));
    Vec4D f1(sqrt(m1*m1 + pcm*pcm),  pcm*n);
    Vec4D f2(sqrt(m2*m2 + pcm*pcm), -pcm*n);
    Poincare boostQp(Qp);
    boostQp.BoostBack(f1); boostQp.BoostBack(f2);
    p[l] = f1; p[r] = f2;
  }
  return RebuildBeamsAtPreFSR(p, k, "FSRMAPD");
}

// The beams at sqrt(Q^2), Q the pre-emission final state (= s'), in the
// same frame convention as MapMomentaFSR: into the Q rest frame, the beams
// along the Born axis there, and back. Shared by the flat and the pole map.
bool NLO_Base::RebuildBeamsAtPreFSR(Vec4D_Vector &p, Vec4D_Vector &k,
                                    const char *tag) {
  Vec4D Q;
  for (size_t i = 2; i < p.size(); ++i) Q += m_plab[i];
  const double sq(Q.Abs2());
  if (!(sq > 0.)) return false;
  Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
  Poincare boostQ(Q);
  for (size_t i = 2; i < p.size(); ++i) { pRot.RotateBack(p[i]); boostQ.Boost(p[i]); }
  for (Vec4D &kj : k) { pRot.RotateBack(kj); boostQ.Boost(kj); }
  const double sign_z = (m_bornMomenta[0][3] < 0 ? -1 : 1);
  const double mb1 = m_flavs[0].Mass(), mb2 = m_flavs[1].Mass();
  const double lamCM = 0.5 * sqrt(Lambda(sq, mb1*mb1, mb2*mb2) / sq);
  p[0] = {sqrt(lamCM*lamCM + mb1*mb1), 0, 0,  sign_z * lamCM};
  p[1] = {sqrt(lamCM*lamCM + mb2*mb2), 0, 0, -sign_z * lamCM};
  Poincare pRot2(m_bornMomenta[0], Vec4D(0., 0., 0, 1.));
  for (size_t i = 0; i < p.size(); ++i) { pRot2.Rotate(p[i]); boostLab.BoostBack(p[i]); }
  for (Vec4D &kj : k) { pRot2.Rotate(kj); boostLab.BoostBack(kj); }
  { Vec4D res(p[0] + p[1]);
    for (size_t i = 2; i < p.size(); ++i) res -= p[i];
    for (const Vec4D &kj : k) res -= kj;
    const double scale(Max(1., p[0][0]));
    if (Vec3D(res).Abs() > 1e-9*scale || std::abs(res[0]) > 1e-9*scale) {
      static long nprint(0);
      if (++nprint <= 20)
        std::cerr<<std::setprecision(10)<<"@@@ "<<tag<<" residual="<<res
                 <<" nk="<<k.size()
                 <<" nFSR="<<m_FSRPhotons.size()<<" P="<<(m_bornMomenta[0]+m_bornMomenta[1])
                 <<" Q="<<Q<<std::endl;
    } }
  return true;
}

/*
  The (n+1)-body point at which beta_1(k_j) is evaluated when the event has
  OTHER photons, for a photon radiated from the initial state (REAL_MAP: 1).

  What the point has to reproduce is the soft-photon radiation pattern of
  the event: the weight is beta_1(k_j)/S~(k_j) with S~ the eikonal of the
  EVENT, so |M_1|^2 at X_j must carry the same collinear structure for k_j
  or the ratio is unbounded. The massless eikonal sum_i Q_i p_i/(p_i.k) is
  invariant under a rescaling p_i -> x_i p_i of every leg and, for the beams,
  under boosts along the beam axis (it is 4/k_perp^2), but NOT under a
  rotation or a transverse boost of k relative to the axis.

  The rest-frame construction below (REAL_MAP: 0) boosts final state + k_j
  into their common rest frame and rebuilds the beams along the lab axis
  THERE. With a hard spectator photon that frame moves at beta ~ 0.7 and the
  aberrated k_j can land on the rebuilt beam axis: measured on a 250 GeV
  radiative-return event, a 7.4 GeV photon at 34 degrees to the beams (lab
  eikonal 4.3e-5) came out at 20 mrad (eikonal 1.1e-2), and a -21% real
  correction became -52 times the Born. The tail of that ratio is 1/f, so
  the cross-section error never shrank with statistics (2.5-3.5% at 100k for
  250 GeV nu nu, mu mu, tau tau; per-event factors down to -2600).

  Here instead:
    - k_j and the beam DIRECTIONS stay as in the lab;
    - the other photons are taken out of the beams through their
      longitudinal projection, P~ = P - (E_K, (K.n) n), which is the
      light-cone removal that leaves x_a P_a + x_b P_b on the axis;
    - one common factor lambda fixes (lambda P~ - k_j)^2 = s', so the pair
      invariant - the Z propagator - is exactly the event's and matches the
      Born of the subtraction;
    - the final state is boosted rigidly from Q_lab to Q' = P' - k_j, so its
      internal invariants and masses are untouched.
  The II eikonal is then preserved up to the electron-mass terms (whose
  ratio lies in [x^4, 1]: the dead cone only shrinks the term), and at n = 1
  there is nothing to remove, lambda = 1, and the point is the event itself,
  so the one-photon stream is unchanged. What is given up is the angle of
  k_j to the final-state fermions, which changes by the bounded Doppler
  factor of the Q_lab -> Q' boost; for a photon radiated from the FINAL
  state that is the singular structure, so those keep the rest-frame path.
*/
bool NLO_Base::MapMomentaBeamAxis(Vec4D_Vector &p, const Vec4D_Vector &k) {
  if (k.empty() || m_bornMomenta.size() < 2 || p.size() < 3) return false;
  for (const Vec4D &kj : k) if (PhotonIsFSR(kj)) return false;
  const Vec4D P(m_bornMomenta[0] + m_bornMomenta[1]);
  const Vec3D pa3(m_bornMomenta[0]);
  if (!(pa3.Abs() > 0.)) return false;
  const Vec3D n(pa3/pa3.Abs());
  Vec4D Q, ksel;
  for (size_t i(2); i < p.size(); ++i) Q += p[i];
  for (const Vec4D &kj : k) ksel += kj;
  const double sp(Q.Abs2());
  if (!(sp > 0.)) return false;
  // Everything the event radiated that is not among the selected photons:
  // by conservation, so hidden or unlisted photons are removed as well.
  const Vec4D K(P - Q - ksel);
  const double Kn(Vec3D(K)*n);
  const Vec4D Pt(P - Vec4D(K[0], Kn*n));
  const double Pt2(Pt.Abs2());
  if (!(Pt2 > 0.) || !(Pt[0] > 0.)) return false;
  // lambda^2 Pt^2 - 2 lambda Pt.ksel - (sp - ksel^2) = 0, the positive root.
  const double b(Pt*ksel), c(sp - ksel.Abs2());
  const double disc(b*b + Pt2*c);
  if (!(disc >= 0.)) return false;
  const double lambda((b + sqrt(disc))/Pt2);
  if (!(lambda > 0.) || IsBad(lambda)) return false;
  const Vec4D Pp(lambda*Pt);
  const Vec4D Qp(Pp - ksel);
  if (!(Qp[0] > 0.) || !IsEqual(Qp.Abs2(), sp, 1e-8)) return false;
  // On-shell beams back-to-back along +-n in the P' rest frame; P' has no
  // transverse momentum, so that frame is a longitudinal boost away and the
  // axis survives the boost back.
  const double sj(Pp.Abs2());
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (sj <= sqr(m1 + m2)) return false;
  const double lamCM(0.5*sqrt(Lambda(sj, m1*m1, m2*m2)/sj));
  Vec4D pa(sqrt(lamCM*lamCM + m1*m1),  lamCM*n);
  Vec4D pb(sqrt(lamCM*lamCM + m2*m2), -lamCM*n);
  Poincare toPp(Pp);
  toPp.BoostBack(pa); toPp.BoostBack(pb);
  // The final state, rigidly, from its lab total Q to Q'.
  Poincare fromQ(Q), toQp(Qp);
  Vec4D_Vector q(p.begin() + 2, p.end());
  Vec4D qsum;
  for (Vec4D &qi : q) { fromQ.Boost(qi); toQp.BoostBack(qi); qsum += qi; }
  Vec4D bal(pa + pb - qsum - ksel);
  const double scale(Max(Pp[0], 1.));
  for (int mu(0); mu < 4; ++mu)
    if (dabs(bal[mu]) > 1e-7*scale) {
      msg_Error() << METHOD << "(): momentum imbalance " << bal
                  << " in the beam-axis reduction, using the rest-frame map."
                  << std::endl;
      return false;
    }
  p[0] = pa; p[1] = pb;
  for (size_t i(2); i < p.size(); ++i) p[i] = q[i-2];
  m_map_reduced = (K[0] > 1e-12*P[0]);
  { static const bool dm(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>()!=0);
    if (dm) std::cerr<<"@@@ MAPQ sqrt_sqq="<<sqrt(sj)
                     <<" sqrt_s="<<sqrt(m_s)
                     <<" sqrt_sborn="<<(m_bornMomenta[2]+m_bornMomenta[3]).Mass()
                     <<" Eksum="<<ksel[0]<<" lambda="<<lambda<<" beamaxis=1"<<std::endl; }
  return true;
}

/*
  The KKMC convention for beta_1(k_j) in an n-photon event (REAL_MAP: 2).

  KKMC's EEX beta_1 for an initial-state photon is S~(k) B(s') times
  [(1-a)^2 + (1-b)^2]/2 - 1 with a = k.p_a/(p_a.p_b), b = k.p_b/(p_a.p_b)
  taken with the FULL beams and B at the event's s': the other photons enter
  only through s'. The beam-axis reduction above judges the photon against
  beams already reduced by the others, a' = a/x_b, which for two hard photons
  counts the energy loss twice. Leading-log check, two collinear ISR photons
  with fractions x_1, x_2 of the beam, exact factor
  (1+(1-x_1)^2)/2 * (1+(1-x_2/(1-x_1))^2)/2 symmetrised:
      x_1 = x_2 = 0.44 : exact 0.34, KKMC 0.32, reduced beams 0.04
      x_1 = 0.8, x_2 = 0.05 : exact 0.45, KKMC 0.47, reduced beams 0.29
  so the KKMC convention is the better O(alpha) truncation (the remainder is
  the genuine beta_2), and it is what YFS.EEX and CEEX effectively carry.
  Measured, 250 GeV nu nu, 50k events, same sample: the beam-axis reduction
  gives YFS.BR 2.838 +- 0.3%, this construction 3.166 +- 0.2%, CEEX (no
  virtual) 3.272 +- 0.4%, EEX 3.513 (with its +7% virtual); the rest-frame
  map gave 2.925 +- 10%. The remaining 3% to CEEX sits in the events where
  two hard photons share the radiative return, where CEEX's coherent
  amplitude-level sum carries an effective beta_2 that O(alpha) cannot.

  A momentum-conserving point with the full-beam a, b and the pair at s' is
  the n = 1 configuration {P_a, P_b, k_j, Q} scaled by one factor x:
      x^2 (P - k_j)^2 = s'   ->   x = sqrt(s'/(P - K_sel)^2) <= 1,
  a and b are scale invariant, the pair sits at s' by construction, and the
  final state is boosted rigidly from Q_lab to Q' = x(P - k_j). The photon in
  the point is x k_j, so the eikonal of the point is the event's over x^2 and
  the SAME 1/x^2 sits in |M_1|^2: CalculateReal therefore takes the ratio
  |M_1|^2 flux/(S~ B) at the point and multiplies the event's own eikonal
  (m_map_reduced). At n = 1, x = 1 and nothing changes. What is not preserved
  is (Q + k_j)^2, the invariant of the final-state emission diagrams of the
  full real ME; for photons radiated from the final state that IS the
  resonance, so those keep the rest-frame construction.
*/
bool NLO_Base::MapMomentaScaled(Vec4D_Vector &p, Vec4D_Vector &k, bool anylabel) {
  if (k.empty() || m_bornMomenta.size() < 2 || p.size() < 3) return false;
  if (!anylabel) for (const Vec4D &kj : k) if (PhotonIsFSR(kj)) return false;
  const Vec4D P(m_bornMomenta[0] + m_bornMomenta[1]);
  const Vec3D pa3(m_bornMomenta[0]);
  if (!(pa3.Abs() > 0.)) return false;
  const Vec3D n(pa3/pa3.Abs());
  Vec4D Q, ksel;
  for (size_t i(2); i < p.size(); ++i) Q += p[i];
  for (const Vec4D &kj : k) ksel += kj;
  const double sp(Q.Abs2());
  const Vec4D R(P - ksel);
  const double R2(R.Abs2());
  if (!(sp > 0.) || !(R2 > 0.) || !(R[0] > 0.)) return false;
  const Vec4D K(P - Q - ksel);
  // x <= 1 up to rounding: (P - K_sel)^2 = (Q + K_others)^2 >= Q^2.
  const double x(Min(1., sqrt(sp/R2)));
  if (!(x > 0.) || IsBad(x)) return false;
  const double sj(x*x*P.Abs2());
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (sj <= sqr(m1 + m2)) return false;
  // Beams x P_a, x P_b rebuilt on shell along +-n in the P rest frame, then
  // returned to the lab (the identity for balanced beams).
  const double lamCM(0.5*sqrt(Lambda(sj, m1*m1, m2*m2)/sj));
  Vec4D pa(sqrt(lamCM*lamCM + m1*m1),  lamCM*n);
  Vec4D pb(sqrt(lamCM*lamCM + m2*m2), -lamCM*n);
  Poincare toP(P);
  toP.BoostBack(pa); toP.BoostBack(pb);
  const Vec4D Qp(pa + pb - x*ksel);
  if (!(Qp[0] > 0.) || !(Qp.Abs2() > 0.)) return false;
  Poincare fromQ(Q), toQp(Qp);
  Vec4D_Vector q(p.begin() + 2, p.end());
  Vec4D qsum;
  for (Vec4D &qi : q) { fromQ.Boost(qi); toQp.BoostBack(qi); qsum += qi; }
  Vec4D bal(pa + pb - qsum - x*ksel);
  const double scale(Max(P[0], 1.));
  for (int mu(0); mu < 4; ++mu)
    if (dabs(bal[mu]) > 1e-7*scale) {
      msg_Error() << METHOD << "(): momentum imbalance " << bal
                  << " in the scaled reduction, using the rest-frame map."
                  << std::endl;
      return false;
    }
  p[0] = pa; p[1] = pb;
  for (size_t i(2); i < p.size(); ++i) p[i] = q[i-2];
  for (Vec4D &kj : k) kj = x*kj;
  m_map_reduced = (K[0] > 1e-12*P[0]);
  { static const bool dm(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>()!=0);
    if (dm) std::cerr<<"@@@ MAPQ sqrt_sqq="<<sqrt(sj)
                     <<" sqrt_s="<<sqrt(m_s)
                     <<" sqrt_sborn="<<(m_bornMomenta[2]+m_bornMomenta[3]).Mass()
                     <<" Eksum="<<ksel[0]<<" x="<<x<<" scaled=1"<<std::endl; }
  return true;
}

/*
  REAL_MAP: 3. The scaled construction keeps the photon's fractions a, b
  relative to the full beams but moves (Q + k_j)^2, the invariant mass of
  final state plus photon. For a charged final state the full real ME also
  has the diagrams with the photon radiated from the final state, whose Z
  propagator sits at exactly that invariant: with the Z-pole Born, moving it
  by a few GeV changes |M_1|^2 by orders of magnitude. Measured: 250 GeV
  mu mu with REAL_MAP: 2 kept the central value (BR 3.21, CEEX 3.22) but a
  5% error from 466 events with |BR| > 10, and the Z-pole mu mu Born+real
  moved by 0.5%.

  Here both are kept - and the measurement is the reason this is NOT the
  default: 250 GeV mu mu, 100k, REAL_MAP: 3 gives BR 4.98 +- 26% against
  3.21 +- 4.8% and 2.82 +- 5.4% (two versions of the denominator, same
  seeds) for REAL_MAP: 2 (CEEX 3.22) - the 5% is a heavy tail, not a
  Gaussian error, and the central value moves with it. Keeping the invariant exactly
  puts more ISR-labelled events onto the Z pole of their final-state
  emission diagrams, where the full |M_1|^2 is shared between the I and F
  labels in the ratio S~_II/S~_tot but divided by the ISR label's own
  S~_tot B((mu mu)^2), which is orders of magnitude below the physics. The
  label sum is exact (both labels together give the exact real once), so
  this is variance, not bias; the cure is a multi-channel denominator
  S~_II B_I + S~_FF B_F, which needs the Born at the other label's
  kinematics and is not attempted here. At the Z pole modes 2 and 3 agree
  per event to four digits. Unknowns: the beam energies x_a, x_b (directions fixed)
  and one scale kappa for the photon. Conditions:
      P'^2 = (Q + k_j)^2            (final-state-emission invariant),
      (P' - kappa k_j)^2 = s'       (the pair at the event's s'),
      log kappa = w log x_a + (1-w) log x_b,  w = (1 + cos theta_k)/2,
  the last being the scaled convention for the beam the photon is collinear
  with (kappa = x_a exactly along beam a: its fraction of that beam is the
  event's) and a smooth interpolation in between. The photon's energy in the
  (Q + k_j) rest frame is then the event's too, so the final-state emission
  is as hard as it really was; only the angle between photon and fermions
  changes, by the bounded Doppler factor of the Q_lab -> Q' boost. Solved by
  iteration on kappa; the 2x2 system for P' along the axis is quadratic,
  the root with the smaller longitudinal momentum is the one continuous
  with P' = P at n = 1. When it has no real solution (a photon with more
  transverse momentum than energy in the (Q + k_j) frame) the scaled
  construction is used instead.
*/
bool NLO_Base::MapMomentaInvariant(Vec4D_Vector &p, Vec4D_Vector &k) {
  if (k.empty() || m_bornMomenta.size() < 2 || p.size() < 3) return false;
  for (const Vec4D &kj : k) if (PhotonIsFSR(kj)) return false;
  const Vec4D P(m_bornMomenta[0] + m_bornMomenta[1]);
  const Vec3D pa3(m_bornMomenta[0]);
  if (!(pa3.Abs() > 0.)) return false;
  const Vec3D n(pa3/pa3.Abs());
  Vec4D Q, ksel;
  for (size_t i(2); i < p.size(); ++i) Q += p[i];
  for (const Vec4D &kj : k) ksel += kj;
  const double sp(Q.Abs2());
  const Vec4D R(Q + ksel);
  const double sj(R.Abs2());
  if (!(sp > 0.) || !(sj > sp) || !(R[0] > 0.)) return false;
  const Vec4D K(P - R);
  if (!(K[0] > 1e-12*P[0])) return MapMomentaScaled(p, k);   // n = 1: the event
  const double Es(ksel[0]), ks(Vec3D(ksel)*n), D(Es*Es - ks*ks), ksel2(ksel.Abs2());
  if (!(Es > 0.) || !(D > 0.)) return false;
  const double ck(Vec3D(ksel).Abs() > 0. ? ks/Vec3D(ksel).Abs() : 0.);
  const double w(0.5*(1. + ck));
  const double Ea(m_bornMomenta[0][0]), Eb(m_bornMomenta[1][0]);
  const double pza(Vec3D(m_bornMomenta[0])*n), pzb(Vec3D(m_bornMomenta[1])*n);
  // E' = xa Ea + xb Eb, pz' = xa pza + xb pzb (pzb < 0 for balanced beams).
  const double det(Ea*pzb - Eb*pza);
  if (IsZero(det)) return false;
  double kappa(sqrt(sp/(P - ksel).Abs2())), xa(1.), xb(1.), Ep(0.), pzp(0.);
  bool ok(false);
  for (int it(0); it < 50; ++it) {
    const double C((sj + kappa*kappa*ksel2 - sp)/(2.*kappa));   // P'.ksel
    const double disc(C*C - D*sj);
    if (!(disc >= 0.)) break;
    const double r1((C*ks + Es*sqrt(disc))/D), r2((C*ks - Es*sqrt(disc))/D);
    pzp = dabs(r1) < dabs(r2) ? r1 : r2;
    Ep  = (C + pzp*ks)/Es;
    if (!(Ep > 0.) || !(Ep*Ep - pzp*pzp > 0.)) break;
    // P' = xa P_a + xb P_b along the axis.
    xa = (Ep*pzb - Eb*pzp)/det;
    xb = (Ea*pzp - pza*Ep)/det;
    if (!(xa > 0.) || !(xb > 0.)) break;
    const double kn(exp(w*log(xa) + (1. - w)*log(xb)));
    if (dabs(kn - kappa) < 1e-12*kn) { kappa = kn; ok = true; break; }
    kappa = kn;
  }
  if (!ok || IsBad(kappa) || !(kappa > 0.)) return MapMomentaScaled(p, k);
  const Vec4D Pp(Ep, pzp*n);
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (sj <= sqr(m1 + m2)) return MapMomentaScaled(p, k);
  const double lamCM(0.5*sqrt(Lambda(sj, m1*m1, m2*m2)/sj));
  Vec4D pa(sqrt(lamCM*lamCM + m1*m1),  lamCM*n);
  Vec4D pb(sqrt(lamCM*lamCM + m2*m2), -lamCM*n);
  Poincare toPp(Pp);
  toPp.BoostBack(pa); toPp.BoostBack(pb);
  const Vec4D Qp(pa + pb - kappa*ksel);
  if (!(Qp[0] > 0.) || !IsEqual(Qp.Abs2(), sp, 1e-6)) return MapMomentaScaled(p, k);
  Poincare fromQ(Q), toQp(Qp);
  Vec4D_Vector q(p.begin() + 2, p.end());
  Vec4D qsum;
  for (Vec4D &qi : q) { fromQ.Boost(qi); toQp.BoostBack(qi); qsum += qi; }
  Vec4D bal(pa + pb - qsum - kappa*ksel);
  const double scale(Max(P[0], 1.));
  for (int mu(0); mu < 4; ++mu)
    if (dabs(bal[mu]) > 1e-7*scale) {
      msg_Error() << METHOD << "(): momentum imbalance " << bal
                  << ", using the scaled reduction." << std::endl;
      return MapMomentaScaled(p, k);
    }
  p[0] = pa; p[1] = pb;
  for (size_t i(2); i < p.size(); ++i) p[i] = q[i-2];
  for (Vec4D &kj : k) kj = kappa*kj;
  m_map_reduced = true;
  { static const bool dm(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>()!=0);
    if (dm) std::cerr<<"@@@ MAPQ sqrt_sqq="<<sqrt(sj)
                     <<" sqrt_s="<<sqrt(m_s)
                     <<" sqrt_sborn="<<(m_bornMomenta[2]+m_bornMomenta[3]).Mass()
                     <<" Eksum="<<ksel[0]<<" kappa="<<kappa<<" xa="<<xa<<" xb="<<xb
                     <<" invariant=1"<<std::endl; }
  return true;
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D_Vector &k) {
  static const int mapmode(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_MAP"].Get<int>());
  m_map_reduced = false;
  if (mapmode == 1 && MapMomentaBeamAxis(p, k)) return;
  if (mapmode == 2 && MapMomentaScaled(p, k)) return;
  if (mapmode == 3 && MapMomentaInvariant(p, k)) return;
  { static const int fsrmap(ATOOLS::Settings::GetMainSettings()["YFS"]
                            ["REAL_FSR_MAP"].Get<int>());
    if (fsrmap == 2 && MapMomentaFSRDipole(p, k)) return;
    // 2 falls back to the common rescaling when the per-dipole construction
    // does not apply (a photon without a dipole, a pair below threshold).
    if (fsrmap >= 1 && MapMomentaFSR(p, k)) return; }
  if (mapmode != 0) {
    // Only initial-state photons are reduced by the new constructions; a
    // final-state photon (or an early exit) comes here by design. Counted so
    // that a solver that never converges cannot masquerade as the legacy map.
    static long nfall(0);
    if (++nfall == 1000 && msg_LevelIsDebugging())
      msg_Debugging() << METHOD << "(): 1000 photons on the rest-frame map with REAL_MAP "
                      << mapmode << std::endl;
  }
  Vec4D Q;
  Vec4D QQ;
  Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
  for (size_t i = 2; i < p.size(); ++i) Q += p[i];
  Vec4D ksum;
  for (const Vec4D &ki : k) ksum += ki;
  Q += ksum;
  const double sq = Q.Abs2();
  Poincare boostQ(Q);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
  for (size_t i = 0; i < p.size(); ++i) {
    pRot.RotateBack(p[i]);
    boostQ.Boost(p[i]);
  }
  for (Vec4D &ki : k) {
    pRot.RotateBack(ki);
    boostQ.Boost(ki);
  }
  // ksum must be recomputed: the photons have been boosted since.
  ksum = Vec4D();
  for (const Vec4D &ki : k) ksum += ki;
  // CheckMappingRecoil(p, ksum);
  for (size_t i = 2; i < p.size(); ++i) QQ += p[i];
  QQ += ksum;
  const double sqq = QQ.Abs2();
  if (!IsEqual(sqq, sq, 1e-6)) {
    msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD
                << std::endl;
  }
  { static const bool dm(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>()!=0);
    if (dm) std::cerr<<"@@@ MAPQ sqrt_sqq="<<(sqq>0?sqrt(sqq):-1.)
                     <<" sqrt_s="<<sqrt(m_s)
                     <<" sqrt_sborn="<<(m_bornMomenta[2]+m_bornMomenta[3]).Mass()
                     <<" Eksum="<<ksum[0]<<std::endl; }
  const double sign_z = (m_bornMomenta[0][3] < 0 ? -1 : 1);
  const double m1 = m_flavs[0].Mass();
  const double m2 = m_flavs[1].Mass();
  const double lamRaw = sqq*sqq + sqr(m1*m1) + sqr(m2*m2)
                        - 2.*sqq*m1*m1 - 2.*sqq*m2*m2 - 2.*m1*m1*m2*m2;
  if (lamRaw < 0.)
    msg_Error()<<METHOD<<"(): below-threshold Kaellen argument = "<<lamRaw
               <<" (sqq = "<<sqq<<", threshold = "<<sqr(m1+m2)<<")"<<std::endl;
  const double lamCM = 0.5 * sqrt(Lambda(sqq, m1 * m1, m2 * m2) / sqq);
  const double E1 = lamCM * sqrt(1 + m1 * m1 / sqr(lamCM));
  const double E2 = lamCM * sqrt(1 + m2 * m2 / sqr(lamCM));
  p[0] = {E1, 0, 0, sign_z * lamCM};
  p[1] = {E2, 0, 0, -sign_z * lamCM};
  Poincare pRot2(m_bornMomenta[0], Vec4D(0., 0., 0, 1.));
  for (size_t i = 0; i < p.size(); ++i) {
    pRot2.Rotate(p[i]);
    boostLab.BoostBack(p[i]);
  }
  for (Vec4D &ki : k) {
    pRot2.Rotate(ki);
    boostLab.BoostBack(ki);
  }
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k) {
  Vec4D_Vector ks{k};
  MapMomenta(p, ks);
  k = ks[0];
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k,
                          const Vec4D_Vector &ref) {
  const Vec4D_Vector save(m_bornMomenta);
  m_bornMomenta = ref;
  Vec4D_Vector ks{k};
  MapMomenta(p, ks);
  k = ks[0];
  m_bornMomenta = save;
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k1, Vec4D &k2) {
  Vec4D_Vector ks{k1, k2};
  MapMomenta(p, ks);
  k1 = ks[0];
  k2 = ks[1];
}

void NLO_Base::CheckMasses(Vec4D_Vector &p, int realmode) {
  bool allonshell = true;
  std::vector<double> masses;
  Flavour_Vector flavs = m_flavs;
  if (realmode >= 1)
    flavs.push_back(Flavour(kf_photon));
  if (realmode >= 2)
    flavs.push_back(Flavour(kf_photon));
  if (p.size() != flavs.size())
    msg_Error() << "Mismatch between mass and flavour vectors in " << METHOD
                << std::endl;
  for (int i = 0; i < p.size(); ++i) {
    masses.push_back(flavs[i].Mass());
    if (!IsEqual(p[i].Mass(), flavs[i].Mass()) && flavs[i].Mass() != 0) {
      allonshell = false;
    }
  }
  if (!allonshell) {
    m_stretcher.StretchMomenta(p, masses);
    // for (int i = 0; i < p.size(); ++i) {
    // }
  }
}

bool NLO_Base::CheckPhotonForReal(const Vec4D &k) {
  for (int i = 0; i < m_plab.size(); ++i) {
    if (m_flavs[i].IsChargedLepton()) {
      double sik = (k + m_plab[i]).Abs2();
      if (sik  < m_hardmin*m_plab[i].Abs2()) {
        msg_Out() << "Rejecting photon k = " << k << std::endl
                  << "sik = " << sik << std::endl;
        return false;
      }
    }
  }
  return true;
}

bool NLO_Base::CheckPhotonForReal(const Vec4D &k, const Vec4D_Vector &p) {
  for (int i = 0; i < p.size(); ++i) {
    if (m_flavs[i].IsChargedLepton()) {
      double sik = (k + p[i]).Abs2();
      if (sik < m_hardmin * p[i].Abs2()) {
        msg_Out() << "Rejecting photon k = " << k << std::endl
                  << "sik = " << sik << std::endl;
        return false;
      }
      // if(p[i].PPerp() < m_hardmin) return false;
    }
  }
  // if(k.PPerp() < m_hardmin) return false;
  return true;
}

bool NLO_Base::CheckMomentumConservation(Vec4D_Vector p) {
  Vec4D incoming = p[0] + p[1];
  Vec4D outgoing;
  for (int i = 2; i < p.size(); ++i) {
    if (p[i].E() < 0 || IsBad(p[i].E())) {
      msg_Error() << "Energy less than zero!: " << p[i] << std::endl;
      return false;
    }
    outgoing += p[i];
  }
  Vec4D diff = incoming - outgoing;
  if (!IsEqual(incoming, outgoing, 1e-8)) {
    msg_Error() << METHOD << std::endl
                << "Momentum not conserverd in YFS NLO" << std::endl
                << "Incoming momentum = " << incoming << std::endl
                << "Outgoing momentum = " << outgoing << std::endl
                << "Difference = " << diff << std::endl
                << "Vetoing Event " << std::endl;
    return false;
  }
  return true;
}


Vec4D NLO_Base::MostEnergeticPhoton() const {
  Vec4D hardest;
  for (const auto &k : m_ISRPhotons)
    if (k.E() > hardest.E()) hardest = k;
  for (const auto &k : m_FSRPhotons)
    if (k.E() > hardest.E()) hardest = k;
  return hardest;
}


Vec4D NLO_Base::FixedTestPhoton() const {
  double E = m_rv_test_x * sqrt(m_s) / 2.;
  double st = sin(m_rv_test_theta), ct = cos(m_rv_test_theta);
  Vec4D k(E, E * st * cos(m_rv_test_phi), E * st * sin(m_rv_test_phi),
         E * ct);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
  Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
  pRot.Rotate(k);
  boostLab.BoostBack(k);
  return k;
}

