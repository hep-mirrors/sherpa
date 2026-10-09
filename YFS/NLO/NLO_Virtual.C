#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
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
#include "EXTRA_XS/NNLO/ee_aa_QED_VV.H"
#include <iomanip>
#include <sstream>
#include "MODEL/Main/Model_Base.H"
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "YFS/NLO/NLO_Base_Internal.H"

using namespace YFS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

double NLO_Base::CalculateVirtual() {
  m_lhelok = false;
  m_vborn_ok = false;
  m_vewho = 0.;
  m_b02 = YFS::Photonic_Factors();
  m_vp = YFS::VP_Factors();
  m_ghelok = false;
  if (AcraicSuppliesVirtual())
    return (m_acraicvirt - 1.) * m_born;
  if (m_eex_virt) {
    // subtract born to avoid double counting
    // already present in eex!!
    return p_dipoles->CalculateEEXVirtual() * m_born - m_born;
  }
  if (!m_looptool)
    return 0;
  double virt;
  double sub;
  /*
    YFS: VIRTUAL_SUB_LEGS (default own). The Born virtual's IR subtraction
    has to cancel the loop's poles, so it is built on the loop's own point
    m_plab, the post-ISR Born.
  */
  static const virtsublegs::code vsublegs(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["VIRTUAL_SUB_LEGS"].SetDefault(virtsublegs::own).Get<virtsublegs::code>());
  Define_Dipoles *vsub(p_dipoles);
  if (vsublegs == virtsublegs::own) {
    p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
    p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);
    p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
    vsub = p_nlodipoles.get();
  }
  p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
  vsub->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
  CheckMassReg();
  /*
    YFS: VIRTUAL_CANONICAL_FRAME (default true since 2026-10-06): evaluate the Born's one-loop
    in the rest frame of the beams with the beams exactly on the z axis
    (CanonicalBeamFrame, which documents why). The event's Born point after
    ISR has two beams with the same tiny p_T, and there OpenLoops' 2 -> 2
    virtual is off by up to 8% of itself (Z-pole mu mu, full EW, alpha(0):
    0.085 against 0.091 in every rotated frame). Not for ACRAIC_Virtual:
    helicity, whose helicity labels are the lab's.
  */
  static const bool vcanon(ATOOLS::Settings::GetMainSettings()["YFS"]
                           ["VIRTUAL_CANONICAL_FRAME"].SetDefault(true).Get<bool>());
  if (vcanon && !(m_useacraic && m_acraicvirtsrc == acraicvirt::helicity))
    virt = p_virt->CalcInFrame(m_plab, CanonicalBeamFrame(m_plab), m_born);
  else
    virt = p_virt->Calc(m_plab, m_born);
  // (not with YFS: VIRTUAL_MASSIFICATION: the provider's Born is then the
  // massless one at the massless image of the point)
  if (m_check_virt_born && !p_virt->p_massif) {
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
    sub = vsub->CalculateVirtualSub();
  else
    sub = 0;
  m_virt_raw = virt; m_virt_subval = sub * m_born / m_rescale_alpha;
  if (m_acraic_compare && m_virt_sub && p_virt->p_loop_me->Mode() != 1) {
    double ii, ff, ifi;
    vsub->VirtualSubParts(ii, ff, ifi);
    m_virt_sub_ii = ii * m_born / m_rescale_alpha;
    m_virt_sub_ff = ff * m_born / m_rescale_alpha;
    m_virt_sub_if = ifi * m_born / m_rescale_alpha;
  }
  m_oneloop = (virt - sub * m_born / m_rescale_alpha);
  // YFS: VIRTUAL_POLE_PROBE (diagnostic, NLO_Checks.C)
  ProbeVirtualPoles(vsub);
  // YFS: VIRTUAL_CHECK: drop the virtual of a numerically unstable loop
  CheckVirtualStability(virt, sub, vcanon);
  // YFS: VIRTUAL_PROBE_ABOVE (diagnostic, NLO_Checks.C)
  ProbeLargeVirtual(virt, sub);
  // YFS: ACRAIC_Virtual: helicity - the same loop call, resolved by helicity;
  // ACRAIC: REAL_VIRTUAL helicity needs the Born's per-helicity ratios too
  if (m_useacraic && (m_acraicvirtsrc == acraicvirt::helicity ||
                      RealVirtualHelicityRequested()) &&
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
    if (!m_virt_sub)
      sub = vsub->CalculateVirtualSub();
    double p1 = p_virt->E1() * p_virt->m_factor;
    // in the loop's coupling, as m_oneloop subtracts it (USE_MODEL_ALPHA 0:
    // the dipoles are at alpha(0), the loop at the model's alpha)
    double yfspole = vsub->Get_E1() / m_rescale_alpha;
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
  AddEWHigherOrder();
  BuildBeta02();
  BuildDysonVP();
  BuildEWHigherOrderHelicity();
  return m_oneloop;
}

/*
  YFS: EW_HIGHER_ORDER: the electroweak correction beyond the provider's NLO
  (GRIFFIN NNLO - NLO), on the Born-level virtual, so both columns carry it
  (ACRAIC reads this virtual as its external v). It is added AFTER m_vborn, the
  v of the RV_MODE remainder v_{n+1} - v: GRIFFIN is 2->2 only, so there is no
  (n+1)-point counterpart, and putting it into v would make a spurious
  remainder. GRIFFIN depends on s and cos(theta) only, so the frame does not
  matter; its value is normalised like the GRIFFIN-as-provider path (mode 1).
*/
/*
    YFS: VIRTUAL_CHECK (default false) - the Born-level loop a second time,
    in another generic frame (the beams' rest frame on z under a different
    fixed rotation, VirtualCheckFrame), whenever |v| = |V_sub/B| exceeds
    VIRTUAL_CHECK_ABOVE (default 0: every event). If the two loops differ by
    more than VIRTUAL_CHECK_TOL (default 1e-2) in units of the Born, the
    point is numerically unstable and its virtual is dropped (v = 0; counted
    and reported at the end with the dropped share of sum m_born).
*/
void NLO_Base::CheckVirtualStability(double &virt, const double sub, const bool vcanon)
{
  static const bool vchk(ATOOLS::Settings::GetMainSettings()["YFS"]
                         ["VIRTUAL_CHECK"].SetDefault(false).Get<bool>());
  static const double vabove(ATOOLS::Settings::GetMainSettings()["YFS"]
                             ["VIRTUAL_CHECK_ABOVE"].SetDefault(0.).Get<double>());
  static const double vtol(ATOOLS::Settings::GetMainSettings()["YFS"]
                           ["VIRTUAL_CHECK_TOL"].SetDefault(1e-2).Get<double>());
  if (vchk && m_born != 0. && p_virt->p_loop_me->Mode() == 0) {
    m_vLoopBornSum += std::abs(m_born);
    if (std::abs(m_oneloop/m_born) > vabove) {
      ++m_vLoopChecked;
      const double virt2(p_virt->CalcInFrame(m_plab, VirtualCheckFrame(m_plab), m_born));
      if (IsBad(virt2) || IsBad(virt)
          || std::abs(virt2 - virt) > vtol*std::abs(m_born)) {
        ++m_vLoopUnstable;
        m_vLoopUnstableBorn += std::abs(m_born);
        m_virt_raw = sub * m_born / m_rescale_alpha;
        m_oneloop = 0.;
        virt = m_virt_raw;
      }
      // leave the provider at the frame the event used (later calls read it)
      if (vcanon && !(m_useacraic && m_acraicvirtsrc == acraicvirt::helicity))
        p_virt->CalcInFrame(m_plab, CanonicalBeamFrame(m_plab), m_born);
      else
        p_virt->Calc(m_plab, m_born);
    }
  }
}

void NLO_Base::AddEWHigherOrder() {
  if (!p_virt || !p_virt->HasEWHigherOrder() || IsBad(m_oneloop)) return;
  const double ho(p_virt->EWHigherOrder(m_plab)/m_rescale_alpha);
  if (!IsBad(ho)) {
    m_oneloop += ho;
    if (m_born != 0.) m_vewho = ho/m_born;
  }
}

/*
  RV_MODE 1: the Born's v with the subtraction on the Born point's OWN legs,
  v_B = V_fin(m_plab)/B - B~(m_plab legs)/kappa, from the loop CalculateVirtual
  already evaluated..
*/
bool NLO_Base::BornVirtualOnOwnLegs(double &vb) {
  if (!m_looptool || m_born == 0.) return false;
  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
  const double sub(p_nlodipoles->CalculateVirtualSub());
  vb = m_virt_raw/m_born - sub/m_rescale_alpha;
  m_vborn_own_sub = sub/m_rescale_alpha;
  return !IsBad(vb);
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
  // (not with YFS: VIRTUAL_MASSIFICATION: the provider's Born is then the
  // massless one at the massless image of the point)
  if (m_check_virt_born && !p_virt->p_massif) {
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
  // the two-loop pole check (the old one compared the 1/eps of the double
  // virtual with the one loop's): NLO_Base::CheckVVPoles
  if (m_check_poles == 1) CheckVVPoles();
  return 0;
}


void NLO_Base::CheckVVPoles()
{
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]["VV_POLE_CHECK"]
                       .SetDefault(true).Get<bool>());
  static const int probe(ATOOLS::Settings::GetMainSettings()["YFS"]["VV_POLE_PROBE"]
                         .SetDefault(0).Get<int>());
  if (!on || p_vv == NULL || p_virt == NULL) return;
  EXTRAXS::ee_aa_QED_VV *aa(dynamic_cast<EXTRAXS::ee_aa_QED_VV*>(p_vv->p_loop_me.get()));
  if (aa == NULL) return;
  aa->Calc(m_plab);
  if (!aa->Valid()) { ++m_vvp_fail; return; }
  // YFS's B on this event's dipoles, in the Born loop's IR convention
  p_dipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_dipoles->MakeDipoles(m_flavs, m_plab, m_plab);
  p_dipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me.get();
  if (!m_dim_reg) {
    static bool warned(false);
    if (!warned) { warned = true;
      msg_Error()<<METHOD<<"(): YFS: VV_POLE_CHECK needs YFS: Dim_Reg true;"
                 <<" no check.\n"; }
    return;
  }
  p_dipoles->CalculateVirtualSub();
  const double mu(p_virt->p_loop_me->IRscale()), s(aa->S());
  // sigma = 2 B/(alpha/pi) in D = 4 - 2 eps: Sherpa's IR slot is minus the
  // 1/eps coefficient; the finite part moved from mu^2 to s
  const double sm1(-2.*p_dipoles->Get_E1()/m_alpi);
  double s0(2.*p_dipoles->Get_Finite()/m_alpi);
  if (mu > 0.) s0 += sm1*std::log(s/(mu*mu));
  const double me(m_flavs[0].Mass());
  if (!(me > 0.)) { ++m_vvp_fail; return; }
  double m1m[3], m2m[3];
  aa->MassiveSeries(me, m1m, m2m);
  const double want[3] = {sm1, 0.5*sm1*sm1, sm1*m1m[1]};
  const double have[3] = {m1m[0], m2m[0], m2m[1]};
  ++m_vvp_n;
  for (int k(0); k < 3; ++k) {
    const double rel(std::abs(have[k] - want[k])/Max(std::abs(want[k]), 1e-300));
    m_vvp_max[k] = Max(m_vvp_max[k], rel);
    // matching digits, 0 .. 16
    const int d(rel > 0. ? Max(0, Min(16, (int)std::floor(-std::log10(rel)))) : 16);
    ++m_vvp_digits[k][d];
  }
  static int nprobe(0);
  if (nprobe < probe) {
    ++nprobe;
    std::ostringstream o;
    o<<std::setprecision(13)<<"@@@ VVPOLE x="<<aa->X()<<" mu="<<mu
     <<" sigma_m1="<<sm1<<" sigma_0="<<s0
     <<" m1m_m1="<<m1m[0]<<" m1m_0="<<m1m[1]<<" v1="<<m1m[1] - s0
     <<" m2m_m2="<<m2m[0]<<" half_sig2="<<want[1]
     <<" m2m_m1="<<m2m[1]<<" sig_m1m="<<want[2]
     <<" loopprov_fin="<<p_virt->p_loop_me->ME_Finite()<<"\n";
    std::cerr<<o.str();
  }
}

void NLO_Base::ReportVVPoles()
{
  // same decision on every rank, so the reductions stay matched
  if (p_vv == NULL || dynamic_cast<EXTRAXS::ee_aa_QED_VV*>(p_vv->p_loop_me.get()) == NULL)
    return;
  double buf[2] = {(double)m_vvp_n, (double)m_vvp_fail};
  double mx[3] = {m_vvp_max[0], m_vvp_max[1], m_vvp_max[2]};
  double hist[3*17];
  for (int k(0); k < 3; ++k) for (int d(0); d < 17; ++d) hist[17*k + d] = (double)m_vvp_digits[k][d];
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    mpi->Allreduce(buf, 2, MPI_DOUBLE, MPI_SUM);
    mpi->Allreduce(mx, 3, MPI_DOUBLE, MPI_MAX);
    mpi->Allreduce(hist, 3*17, MPI_DOUBLE, MPI_SUM);
  }
#endif
  {
    double aa[2] = {(double)m_b02aa_n, (double)m_b02aa_skip};
#ifdef USING__MPI
    if (mpi->Size() > 1) mpi->Allreduce(aa, 2, MPI_DOUBLE, MPI_SUM);
#endif
    if (aa[0] > 0.)
      msg_Out()<<"YFS: BETA0_2_DIPHOTON - two-loop hard function on "<<(long)aa[0]
               <<" Born points, skipped on "<<(long)aa[1]<<" (m_e^2/min(|t|,|u|) >"
               <<" BETA0_2_DIPHOTON_POWER_CUT)\n";
  }
  if (buf[0] + buf[1] <= 0.) return;
  const char *name[3] = {"one loop  eps^-1 = sigma_-1        ",
                         "two loop  eps^-2 = sigma_-1^2/2    ",
                         "two loop  eps^-1 = sigma_-1 m1m_0  "};
  msg_Out()<<"YFS: VV_POLE_CHECK - two-loop soft poles of the double virtual"
           <<" against YFS's subtraction, "<<(long)buf[0]<<" events ("
           <<(long)buf[1]<<" not checked):\n";
  for (int k(0); k < 3; ++k) {
    int dmin(16);
    for (int d(0); d < 17; ++d) if (hist[17*k + d] > 0.) { dmin = d; break; }
    msg_Out()<<"  "<<name[k]<<" max rel. deviation "<<std::setprecision(3)<<mx[k]
             <<", fewest matching digits "<<dmin<<"\n";
  }
  msg_Out()<<std::flush;
}
