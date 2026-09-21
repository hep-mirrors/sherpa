// Self-checks, counters and debug histograms for NLO_Base.
//
// Split out of NLO_Base.C, which was 1935 lines with a quarter of them given
// to machinery that never runs in production: the *Sub cross-checks, the
// blow-up counters and the debug histograms. Same class, second translation
// unit, so nothing about the physics path changed - only where it is written
// down.
//
// Note this is NOT everything named Check*: CheckMasses,
// CheckMomentumConservation and CheckPhotonForReal are on the live path and
// stay in NLO_Base.C, as does m_ifi_prod, which is a physics reweight rather
// than the m_ifi_* profile counters kept here.

#include "YFS/NLO/NLO_Base.H"
#include "MODEL/Main/Model_Base.H"

#include "ATOOLS/Math/Histogram_2D.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"

#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace ATOOLS;
using namespace YFS;
using namespace std;

// Diagnostic output streams and the residual accumulator, used only by the
// checks below.
std::ofstream out_sub, out_real, out_finite;
std::ofstream out_recola;

struct SubCheckAccumulator {
  size_t n = 0;
  double e_first = 0, e_last = 0, r_first = 0, r_last = 0, r_min = 0, r_max = 0;
  // Optional RV cancellation tracking. For each candidate relative threshold,
  // e_below[t] is the largest photon energy at which the cancellation ratio
  // C = |d1|/max(|rv|,|rvsub|) has already fallen below thr[t]. Since C shrinks
  // as the photon softens, a cut placed at e_below[t] removes exactly the
  // region where fewer than ~-log10(thr[t]) digits survive - the safe home for
  // RV_CANCEL_EPS (=thr[t]) or an equivalent RV_SOFT_CUT.
  bool track_cancel = false;
  static const int NTHR = 4;
  double thr[NTHR] = {1e-8, 1e-10, 1e-12, 1e-14};
  double e_below[NTHR] = {0, 0, 0, 0};
  double c_first = 0, c_last = 0;
  void Add(double e, double r) {
    if (n == 0) {
      e_first = e;
      r_first = r;
      r_min = r_max = r;
    }
    e_last = e;
    r_last = r;
    r_min = ATOOLS::Min(r_min, r);
    r_max = ATOOLS::Max(r_max, r);
    ++n;
  }
  void AddCancel(double e, double C) {
    if (!track_cancel) { c_first = C; track_cancel = true; }
    c_last = C;
    for (int t = 0; t < NTHR; ++t)
      if (C < thr[t] && e > e_below[t]) e_below[t] = e;
  }
  void Print(const std::string &label, double roots = 0.) const {
    if (n == 0) {
      msg_Info() << om::brown << label << ": no points evaluated." << om::reset << "\n";
      return;
    }
    const bool converged = r_last < 1e-3 * ATOOLS::Max(r_first, 1e-300);
    msg_Info() << om::bold << "=== " << label << " ===" << om::reset << "\n"
               << "  points             : " << n << "\n"
               << "  photon energy      : " << e_first << " -> " << e_last << "\n"
               << "  |residual|/Born    : " << (converged ? om::green : om::red)
               << r_first << " -> " << r_last << om::reset
               << "  (range [" << r_min << ", " << r_max << "])\n";
    if (r_first > 0. && !IsZero(r_last))
      msg_Info() << "  suppression factor : " << (r_first / r_last) << "x\n";
    msg_Info() << "  " << (converged ? om::green : om::red) << om::bold
               << (converged
                       ? "OK: residual -> 0 in the soft limit"
                       : "WARNING: residual does not appear to vanish in the "
                         "soft limit - check the subtraction!")
               << om::reset << "\n";
    if (track_cancel) {
      msg_Info() << "  cancellation C     : " << c_first << " -> " << c_last
                 << "  (C = |d1|/max(|rv|,|rvsub|))\n"
                 << "  suggested RV cut placement (keep ~-log10(eps) digits):\n";
      for (int t = 0; t < NTHR; ++t) {
        msg_Info() << "    RV_CANCEL_EPS=" << thr[t] << "  -> cut at E >= "
                   << e_below[t];
        if (roots > 0.)
          msg_Info() << "  (RV_SOFT_CUT ~ " << e_below[t] / roots << ")";
        if (e_below[t] == 0.)
          msg_Info() << "  [C never fell below this in the scanned range]";
        msg_Info() << "\n";
      }
    }
  }
};


void NLO_Base::BookHistograms() {
  if (m_isr_debug || m_fsr_debug) {
    m_histograms2d["IFI_EIKONAL"] = std::make_unique<Histogram_2D>(0, -1., 1., 20, 0, 5., 20);
    m_histograms2d["REAL_SUB"] =
        std::make_unique<Histogram_2D>(0, 0, sqrt(m_s), 200, 0, sqrt(m_s) / 2., 20);
    m_histograms2d["REAL_COLL_RATIO"] =
        std::make_unique<Histogram_2D>(0, 0, 2. * M_PI, 20, 0, sqrt(m_s) / 2., 125);
    m_histograms2d["REAL_COLL_RATIO"] =
        std::make_unique<Histogram_2D>(0, 0, sqrt(m_s) / 2., 125, 0, 10, 20);
    m_histograms2d["REAL_RATIO"] =
        std::make_unique<Histogram_2D>(0, 0, 15, 16, 0, sqrt(m_s) / 2., 200);
    m_histograms2d["Real_Flux"] = std::make_unique<Histogram_2D>(0, 0, 1.1, 50, 80, 100, 100);
    m_histograms1d["k_E"] = std::make_unique<Histogram>(0, 0, sqrt(m_s) / 2, sqrt(m_s) / 2);
    m_histograms1d["k_pt"] = std::make_unique<Histogram>(0, 0, sqrt(m_s) / 2, sqrt(m_s) / 2);
    m_histograms1d["dip_mass"] = std::make_unique<Histogram>(0, 0, sqrt(m_s), sqrt(m_s));
    m_histograms1d["k_E_pass"] =
        std::make_unique<Histogram>(0, 0, sqrt(m_s) / 2, sqrt(m_s) / 2);
    m_histograms1d["k_pt_pass"] =
        std::make_unique<Histogram>(0, 0, sqrt(m_s) / 2, sqrt(m_s) / 2);
    m_histograms1d["dip_mass_pass"] = std::make_unique<Histogram>(0, 0, sqrt(m_s), sqrt(m_s));
    if (!ATOOLS::DirectoryExists(m_debugDIR_NLO))
      ATOOLS::MakeDir(m_debugDIR_NLO);
  }
}

void NLO_Base::WriteHistograms() {
  if (!(m_isr_debug || m_fsr_debug || m_check_real_sub || m_check_poles ||
        m_rv_cancel_hist))
    return;
    string name;
    for (auto &hit : m_histograms2d) {
      name = string(m_debugDIR_NLO) + "/" + hit.first + string(".dat");
      // hit.second->MPISync();
      hit.second->Finalize();
      hit.second->Output(name);
    }
    for (auto &hit : m_histograms1d) {
      name = string(m_debugDIR_NLO) + "/" + hit.first + string(".dat");
      hit.second->MPISync();
      hit.second->Finalize();
      hit.second->Output(name);
    }
    // Owned by the maps; nothing to delete.
}


void NLO_Base::CheckRealSub(Vec4D k, int mode) {
  // if(k.E() < 20) return;
  // k*=100;
  double real;
  std::string filename = "Real_subtracted_";
  std::string filename1 = "Sub_term_";
  std::string filename2 = "Real_ME_";
  for (auto f : m_flavs) {
    filename += f.IDName();
    filename += "_";
    filename1 += f.IDName();
    filename1 += "_";
    filename2 += f.IDName();
    filename2 += "_";
  }
  filename += ".txt";
  filename1 += ".txt";
  filename2 += ".txt";
  if (ATOOLS::FileExists(filename))
    ATOOLS::Remove(filename);
  if (ATOOLS::FileExists(filename1))
    ATOOLS::Remove(filename1);
  if (ATOOLS::FileExists(filename2))
    ATOOLS::Remove(filename2);
  out_finite.open(filename, std::ios_base::app);
  out_sub.open(filename1, std::ios_base::app);
  out_real.open(filename2, std::ios_base::app);
  SubCheckAccumulator acc;
  for (double i = 1; i < 20; i += 0.1) {
    k = k / i;
    real = CalculateReal(k, /*raw*/mode >= 3);
    if (k.E() <= 1e-16)
      break;
    out_finite << k.E() << "," << fabs(real) / m_born << std::endl;
    out_real << k.E() << "," << (m_real)*p_nlodipoles->CalculateFlux(k)
             << std::endl;
    out_sub << k.E() << "," << m_subloc * m_born / m_rescale_alpha << std::endl;
    acc.Add(k.E(), fabs(real) / m_born);
  }
  acc.Print("Real subtraction check (" + filename + ")");
  out_finite.close();
  out_real.close();
  out_sub.close();
  exit(0);
}

void NLO_Base::CheckRealVirtualSub(Vec4D k) {
  // if(k.E() < 20) return;
  // k*=100;
  double real;
  std::string filename = "RealVirtual_subtracted";
  std::string filename1 = "SubRV_term";
  std::string filename2 = "RV_ME";
  for (auto f : m_flavs) {
    filename += "_";
    filename += f.IDName();
    filename1 += "_";
    filename1 += f.IDName();
    filename2 += "_";
    filename2 += f.IDName();
  }
  filename += ".txt";
  filename1 += ".txt";
  filename2 += ".txt";
  if (ATOOLS::FileExists(filename))
    ATOOLS::Remove(filename);
  if (ATOOLS::FileExists(filename1))
    ATOOLS::Remove(filename1);
  if (ATOOLS::FileExists(filename2))
    ATOOLS::Remove(filename2);
  out_finite.open(filename, std::ios_base::app);
  out_sub.open(filename1, std::ios_base::app);
  out_real.open(filename2, std::ios_base::app);
  SubCheckAccumulator acc;
  for (double i = 1; i < 20; i += 0.01) {
    k = k / i;
    if (k.E() <= m_isrcut*sqrt(m_s))
      break;
    // Run this tuning scan with RV_CANCEL_EPS/RV_SOFT_CUT disabled so soft
    // points are not skipped - CalculateRealVirtual still fills m_rv/m_rvsub
    // before any guard, but a guard would zero the returned residual here.
    real = CalculateRealVirtual(k);
    if( IsBad(real) || IsBad(m_rv) || IsBad(m_rvsub) ) continue;
    // Relative cancellation of the two diverging RV pieces at this energy.
    const double rvmax = ATOOLS::Max(fabs(m_rv), fabs(m_rvsub));
    const double C = (rvmax > 0.) ? fabs(m_rv - m_rvsub) / rvmax : 0.;
    out_finite << std::setprecision(16) << k.E() << "," << fabs(real) / m_born
               << "," << C << std::endl;
    out_real << std::setprecision(16) << k.E() << "," << m_rv << std::endl;
    out_sub << std::setprecision(16) << k.E() << "," << m_rvsub << std::endl;
    acc.Add(k.E(), fabs(real) / m_born);
    if (C > 0.) acc.AddCancel(k.E(), C);
  }
  acc.Print("Real-Virtual subtraction check (" + filename + ")", sqrt(m_s));
  out_finite.close();
  out_sub.close();
  out_real.close();
  exit(0);
}

void NLO_Base::CheckRealRealSub(Vec4D k1, Vec4D k2) {
  // if(k.E() < 20) return;
  // k*=100;
  double real;
  Vec4D _k1 = k1;
  Vec4D _k2 = k2;
  std::string filename1 = "RealReal_k1_subtracted_";
  std::string filename2 = "RealReal_k2_subtracted_";
  std::string filename3 = "RealReal_k1_k2_subtracted_";
  for (auto f : m_flavs) {
    filename1 += f.IDName();
    filename2 += f.IDName();
    filename3 += f.IDName();
  }
  filename1 += ".txt";
  filename2 += ".txt";
  filename3 += ".txt";
  if (ATOOLS::FileExists(filename1))
    ATOOLS::Remove(filename1);
  if (ATOOLS::FileExists(filename2))
    ATOOLS::Remove(filename2);
  if (ATOOLS::FileExists(filename3))
    ATOOLS::Remove(filename3);
  out_sub.open(filename1, std::ios_base::app);
  // Run these soft-limit scans with RR_SOFT_CUT disabled, as for
  // RV_CANCEL_EPS/RV_SOFT_CUT in CheckRealVirtualSub: the whole point here is
  // to walk the photons into the soft region, which is exactly what the guard
  // skips. With it on, CalculateRealReal returns 0 for the soft points, and the
  // third loop below (which breaks on real==0) would terminate immediately.
  SubCheckAccumulator acc1, acc2, acc12;
  for (double i = 1; i < 20; i += 0.02) {
    k1 = k1 / i;
    if(k1.E()< m_isrcut*sqrt(m_s)) break;
    real = CalculateRealReal(k1, k2);
    out_sub << k1.E() << "," << fabs(real) / m_born << std::endl;
    acc1.Add(k1.E(), fabs(real) / m_born);
  }
  acc1.Print("RealReal subtraction check, k1 -> 0 (" + filename1 + ")");
  out_sub.close();
  out_sub.open(filename2, std::ios_base::app);
  k2 = _k2;
  k1 = _k1;
  for (double i = 1; i < 20; i += 0.02) {
    k2 = k2 / i;
    // if(k2.E() <= 1e-16 ) break;
    if(k2.E()< m_isrcut*sqrt(m_s)) break;
    real = CalculateRealReal(k1, k2);
    out_sub << k2.E() << "," << fabs(real) / m_born << std::endl;
    acc2.Add(k2.E(), fabs(real) / m_born);
    // if(IsZero(real)) break;
    // m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
  }
  acc2.Print("RealReal subtraction check, k2 -> 0 (" + filename2 + ")");
  out_sub.close();
  out_sub.open(filename3, std::ios_base::app);
  k2 = _k2;
  k1 = _k1;
  for (double i = 1; i < 20; i += 0.02) {
    k2 = k2 / i;
    k1 = k1 / i;
    if(k1.E()< m_isrcut*sqrt(m_s)) break;
    real = CalculateRealReal(k1, k2);
    out_sub << k1.E() << "," << fabs(real) / m_born << std::endl;
    acc12.Add(k1.E(), fabs(real) / m_born);
    if (k1.E() <= 1e-16 || real == 0 && !m_failcut)
      break;
    // if(IsZero(real)) break;
    // m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
  }
  acc12.Print("RealReal subtraction check, k1&k2 -> 0 (" + filename3 + ")");
  out_sub.close();
  exit(0);
}

void NLO_Base::CheckMassReg() {
  double virt;
  if (m_check_mass_reg == 1 && !m_realvirt) {
    out_sub.open("yfs-sub.txt", std::ios_base::app);
    out_recola.open("virtual-res.txt",
                    std::ios_base::app); // append instead of overwrite
    out_finite.open("yfs-finite.txt", std::ios_base::app);
    if (!HasISR())
      virt = p_virt->Calc(m_bornMomenta, m_born);
    else
      virt = p_virt->Calc(m_plab, m_born);
    if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born() * m_rescale_alpha,
                 1e-6)) {
      msg_Error() << METHOD
                  << "\n Warning! Loop provider's born is different! YFS "
                     "Subtraction likely fails\n"
                  << "Loop Provider " << ":  " << p_virt->p_loop_me->ME_Born()
                  << "Sherpa" << ":  " << m_born;
    }
    double sub = p_dipoles->CalculateVirtualSub();
    std::cout << setprecision(15);
    out_sub << setprecision(15) << m_photonMass << ","
            << -sub * m_born / m_rescale_alpha << std::endl;
    out_recola << setprecision(15) << m_photonMass << "," << virt << std::endl;
    out_finite << setprecision(15) << m_photonMass << ","
               << virt - sub * m_born / m_rescale_alpha << std::endl;
    out_sub.close();
    out_recola.close();
    exit(0);
  }
}

void NLO_Base::RecordSubScatter(const Vec4D &k, double residual,
                                const std::string &tag, double eik) {
  if (!m_subscatter.is_open()) {
    // per-rank filename: every MPI rank accumulates its own file (concatenate
    // for plotting) so ranks do not clobber a shared file.
    int rank = 0;
#ifdef USING__MPI
    if (mpi->Size() > 1) rank = mpi->Rank();
#endif
    std::string fn = "sub_angle_energy";
    for (auto f : m_flavs) { fn += "_"; fn += f.IDName(); }
    fn += "_rank" + std::to_string(rank) + ".txt";
    if (ATOOLS::FileExists(fn)) ATOOLS::Remove(fn);
    m_subscatter.open(fn.c_str(), std::ios_base::app);
    m_subscatter << "# tag  E_gamma  pT_gamma  theta_charged[rad]  nearest_kf  "
                 << "residual  eik(S~product)\n";
  }
  // angle of the photon k to the nearest charged particle
  const Vec4D_Vector &mom =
      (m_reallab.size() >= m_flavs.size() ? m_reallab : m_plab);
  double th(-1.);
  int ni(-1);
  for (size_t j(0); j < m_flavs.size() && j < mom.size(); ++j) {
    if (m_flavs[j].IntCharge() == 0) continue;
    double t = k.Theta(mom[j]);
    if (th < 0. || t < th) { th = t; ni = j; }
  }
  m_subscatter << std::setprecision(10) << tag << " " << k.E() << " "
               << k.PPerp() << " " << th << " "
               << (ni >= 0 ? (m_flavs[ni].IsAnti() ? -(long)m_flavs[ni].Kfcode()
                                                   : (long)m_flavs[ni].Kfcode())
                           : 0)
               << " " << residual << " " << eik << "\n";
}

void NLO_Base::CEEXComparePoint() {
  if (!m_ceex_compare || m_ceex_done) return;
  if (!m_realtool || !m_looptool) {
    msg_Error() << METHOD << ": CEEX_Compare needs both a real and a loop "
                << "provider (NLO_Part: BVR); got real=" << m_realtool
                << " loop=" << m_looptool << "\n";
    m_ceex_done = true;
    return;
  }
  m_ceex_done = true;

  Vec4D k = FixedTestPhoton();
  msg_Out() << std::setprecision(15)
            << "\n=== Sherpa YFS NLO point for the KKMC CEEX comparison ===\n";
  for (size_t i(0); i < m_bornMomenta.size(); ++i)
    msg_Out() << "  born[" << i << "] (" << m_flavs[i] << ") = "
              << m_bornMomenta[i] << "\n";

  /*
    The REAL-EMISSION configuration is what has to be handed to KKMC, not the
    Born momenta plus a photon: born+k does not conserve momentum.

    Note what this point IS. MapMomenta rebuilds the beams from
    Q = sum(p_out) + k. In production k is one of the event's own photons and
    m_plab already carries its recoil, so Q reconstructs the collider beams
    (the "@@@ PHOT" check asserts photon energy + outgoing energy = sqrt(s)).
    FixedTestPhoton() is SYNTHETIC and is not in m_plab, so here the beams come
    out ABOVE the nominal collider energy - 45.6 GeV becomes 51.99 at x=0.3.
    That is still a perfectly physical configuration: a collider at
    sqrt(s)=103.98 radiating a photon down to a 91.2 hard system. It is simply
    not the collider KKMC is initialised for, so the KKMC side must be told
    CMSene = the sqrt(s) printed below, not the nominal beam energy.
  */
  Vec4D_Vector pmap(m_plab);
  MapMomenta(pmap, k);
  pmap.push_back(k);
  msg_Out() << "  --- real-emission momenta actually used (feed THESE to KKMC) ---\n";
  for (size_t i(0); i < pmap.size(); ++i)
    msg_Out() << "  p[" << i << "] = " << pmap[i] << "\n";
  Vec4D bal(pmap[0]+pmap[1]);
  for (size_t i(2); i < pmap.size(); ++i) bal -= pmap[i];
  msg_Out() << "  sqrt(s) of THIS point = " << (pmap[0]+pmap[1]).Mass()
            << "   <-- set KKMC CMSene to this, not the beam energy\n";
  msg_Out() << "  balance (in - out) = " << bal
            << "   max|component| = "
            << Max(Max(dabs(bal[0]),dabs(bal[1])),Max(dabs(bal[2]),dabs(bal[3])))
            << "\n";
  /*
    Export EVERY input KKMC needs, in KKMC's own ReaData format, so that
    nothing has to be matched by hand and the two generators cannot silently
    disagree on a parameter.

    KKee2f::Initialize reads KKMCee_defaults, then ./pro.input with a NEGATIVE
    imax, which means "overwrite, do not zero". A third file read the same way
    therefore overrides both, and this is that file: the KKMC driver applies it
    last. No transcription, no remembering which knob is which index.

    This exists because chasing a disagreement down to the matrix element cost
    a day of comparing things that were never comparable: Sherpa ran ISR+FSR
    against KeyFSR=0, Sherpa's Gamma_Z was the PDG 2.4952 against KKMC's Dizet
    2.5007203 (26% of the residual soft-limit difference), and SIN2THETAW
    differs by 3.9% (0.23155 vs 0.22276773). Sherpa now states all of it.

    Note the sqrt(s) written here is the sqrt(s) of THIS POINT, not the beam
    energy: FixedTestPhoton() is synthetic and not in m_plab, so the mapped
    configuration sits above the nominal collider energy. It is a physical
    point, just at a different sqrt(s) - and KKMC has to be told which.
  */
  {
    const double sqrtsPoint((pmap[0] + pmap[1]).Mass());
    double sin2tw(0.23155);
    if (MODEL::s_model)
      sin2tw = MODEL::s_model->ComplexConstant("csin2_thetaW").real();
    const long kffin(m_flavs[2].Kfcode());
    std::ofstream kp("ceex_kkmc.input");
    kp << std::setprecision(17) << std::scientific;
    kp << "BeginX  <- written by NLO_Base::CEEXComparePoint, do not edit\n";
    kp << "*  Sherpa's inputs, overriding KKMCee_defaults and pro.input.\n";
    kp << "    1  " << sqrtsPoint  << "  CMSene: sqrt(s) OF THIS POINT\n";
    kp << "   30  " << 1./m_alpha  << "  Alfinv0 = 1/alpha(0)\n";
    kp << "  502  " << Flavour(kf_Z).Mass()  << "  MZ\n";
    kp << "  503  " << sin2tw                << "  swsq = sin^2(theta_W)\n";
    kp << "  504  " << Flavour(kf_Z).Width() << "  GamZ\n";
    kp << "   12  " << 0.0 << "  KeyELW: QED only, matching Sherpa's YFS\n";
    kp << "   20  " << (HasISR() ? 1.0 : 0.0) << "  KeyISR\n";
    kp << "   21  " << (HasFSR() ? 1.0 : 0.0) << "  KeyFSR\n";
    kp << "   27  " << ((HasISR() && HasFSR()) ? 2.0 : 0.0) << "  KeyINT (IFI)\n";
    // Final state: index 400+KF, one flag per channel, all others off.
    for (long kf(1); kf <= 16; ++kf) {
      if (kf > 6 && kf < 11) continue;             // no such channel
      kp << "  " << (400 + kf) << "  " << (kf == kffin ? 1.0 : 0.0)
         << "  KFfin " << kf << (kf == kffin ? "  <== selected" : "") << "\n";
    }
    kp << "EndX\n";
    msg_Out() << "  wrote ceex_kkmc.input (KKMC ReaData format): sqrt(s)="
              << sqrtsPoint << ", MZ=" << Flavour(kf_Z).Mass()
              << ", GamZ=" << Flavour(kf_Z).Width() << ", swsq=" << sin2tw
              << ", 1/alpha=" << 1./m_alpha << ", KFfin=" << kffin << "\n";
  }

  // Machine-readable copy so the KKMC driver can consume the point directly
  // instead of it being transcribed by hand.
  std::ofstream pt("ceex_point.dat");
  pt << std::setprecision(17);
  for (size_t i(0); i < pmap.size(); ++i)
    pt << pmap[i][0] << " " << pmap[i][1] << " "
       << pmap[i][2] << " " << pmap[i][3] << "\n";

  // beta_0 is the Born; beta_1 the O(alpha) real + virtual on top of it. These
  // are the same calls the nominal weight uses, so nothing here is a
  // re-derivation of the physics.
  const double b0    = m_born;
  const double real  = CalculateReal(k);
  const double virt  = CalculateVirtual();
  const double beta1 = real + virt;

  msg_Out() << "  beta_0  (Born)                = " << b0    << "\n"
            << "  real    (CalculateReal)       = " << real  << "\n"
            << "  virtual (CalculateVirtual)    = " << virt  << "\n"
            << "  beta_1  = real + virtual      = " << beta1 << "\n"
            << "  beta_0 + beta_1               = " << b0+beta1 << "\n"
            << "\n  KKMC counterpart (kkmc_ceex_crosscheck):\n"
            << "    beta_0        <-> RhoExp0\n"
            << "    beta_0+beta_1 <-> RhoExp1\n"
            << "    beta_1        <-> RhoExp1 - RhoExp0\n"
            << "    real  only    <-> the 'Born+real' row minus its RhoExp0\n"
            << "    virt  only    <-> the 'Born+virtual' row minus its RhoExp0\n"
            << "  NB both sides carry their own overall normalisation (Sherpa's\n"
            << "  beta_0 is the Born ME, KKMC's RhoExp0 the CEEX distribution),\n"
            << "  so compare the RATIOS beta_1/beta_0 vs (RhoExp1-RhoExp0)/RhoExp0,\n"
            << "  not the absolute numbers.\n";
  msg_Out() << "  beta_1/beta_0 = " << (IsZero(b0) ? 0. : beta1/b0) << "\n\n";
  pt << "# beta0 real virtual\n"
     << b0 << " " << real << " " << virt << "\n";
  pt.close();
}
