// Exercises the SHIPPING real-emission subtraction code (Define_Dipoles::
// CalculateRealSub / CalculateRealSubEEX / CalculateRealSubIF) on a fixed
// phase-space point, and checks it against the eikonal algebra independently.
//
// Why this exists
// ---------------
// The unit tests in Tests/Unit/test_YFS_IFI_Eikonal.C fix the reference algebra
// (coherent = incoherent + interference, the +/- charge pattern, oddness of the
// interference), but they compute it themselves - they do not call
// Define_Dipoles, so they cannot catch a regression in the code that actually
// runs. This driver closes that gap in the same standalone style as
// FSR_KKMC_CrossCheck.C / IFI_KKMC_CrossCheck.C.
//
// What it checks
// --------------
//   1. CalculateRealSub(k) is the COHERENT soft factor: it accumulates one
//      current vector across the II and FF dipoles (+Q initial, -Q final) and
//      squares it at the end (Define_Dipoles.C:465), so the square of the sum
//      carries the ISR x FSR cross term. Compared here against an independently
//      built -|J_ini + J_fin|^2.
//   2. CalculateRealSubEEX(k) is the INCOHERENT one: a per-dipole sum of
//      Dipole::Eikonal over II and FF, with the IF loop commented out
//      (Define_Dipoles.C:872), so it carries no interference.
//   3. Their difference is therefore the IFI real emission. It must be ODD
//      under reflecting the photon, since the interference is what generates
//      the forward-backward asymmetry - the same invariant the unit tests fix.
//   4. CalculateRealSubIF(k) sums the IF-dipole eikonals directly. Its only
//      caller in the shipping code is the IFI_EIKONAL debug histogram
//      (NLO_Base.C:599); printed here so its normalisation relative to the
//      difference in (3) is visible rather than assumed.
//
// The point of (1)+(2)+(3) together: the real subtraction ALREADY contains the
// interference. Any fix that adds it again double-counts.
//
// Build/run: ./build_IFI_RealSub_Test.sh

#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "PDF/Main/PDF_Base.H"
#include "YFS/Main/Define_Dipoles.H"
#include "MODEL/Main/Model_Base.H"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace ATOOLS;
using namespace YFS;

class FakeModel : public MODEL::Model_Base {
public:
  FakeModel(double alpha) : MODEL::Model_Base(true) {
    (*p_constants)["alpha_QED"] = alpha;
  }
  void InitVertices() override {}
  void ParticleInit() override {}
  bool ModelInit() override { return true; }
};

namespace {
  // Reference eikonal current for one leg: eta = +1 incoming, -1 outgoing.
  Vec4D J(double q, double eta, const Vec4D &p, const Vec4D &k) {
    return (q * eta / (p * k)) * p;
  }
  double MinusSq(const Vec4D &v) { return -v.Abs2(); }
}

int main(int argc, char* argv[]) {
#ifdef USING__MPI
  MPI_Init(&argc, &argv);
#endif
  ATOOLS::mpi = new My_MPI();
  ATOOLS::exh = new Terminator_Object_Handler();
  ATOOLS::msg = new Message();
  ATOOLS::rpa = new Run_Parameter();

  // MODE must be ISRFSR: Define_Dipoles::MakeDipolesIF returns immediately for
  // pure-FSR mode and HasFSR() gates the IF dipoles entirely.
  char override_arg[] = "YFS: {MODE: ISRFSR, USE_MODEL_ALPHA: 0, IR_CUTOFF: 1e-9}";
  char prog_name[] = "sherpa_ifi_realsub";
  char* fake_argv[] = {prog_name, override_arg};
  Settings::InitializeMainSettings(2, fake_argv);
  ATOOLS::ran = new Random(1234);
  ATOOLS::s_loader = new Library_Loader();
  PDF::pdfdefs = new PDF::PDF_Defaults();

  try {

  const double me    = 0.000511;
  const double mmu   = 0.105658375;
  const double sqrts = 0.7;
  const double pin   = sqrt(sqr(sqrts/2.) - sqr(me));
  const double pout  = sqrt(sqr(sqrts/2.) - sqr(mmu));
  const double ct = cos(0.7), st = sin(0.7);

  rpa->gen.SetEcms(sqrts);

  s_kftable[11] = new Particle_Info(11, me,  0., -3, 0, 1, "e-",  "e^-");
  s_kftable[13] = new Particle_Info(13, mmu, 0., -3, 0, 1, "mu-", "\\mu^-");
  const double alpha = 1./137.035999084;
  MODEL::s_model = new FakeModel(alpha);

  Flavour flEm((kf_code)11, false), flEp((kf_code)11, true);
  Flavour flMm((kf_code)13, false), flMp((kf_code)13, true);

  Vec4D pEm(sqrts/2., 0., 0.,  pin);
  Vec4D pEp(sqrts/2., 0., 0., -pin);
  Vec4D pMm(sqrts/2.,  pout*st, 0.,  pout*ct);
  Vec4D pMp(sqrts/2., -pout*st, 0., -pout*ct);

  Flavour_Vector fl{flEm, flEp, flMm, flMp};
  Vec4D_Vector   mom{pEm, pEp, pMm, pMp};

  Define_Dipoles dip;

  std::cout << std::setprecision(12);
  std::cout << "=== Define_Dipoles real-subtraction test, e+e- -> mu+mu- ===\n";
  std::cout << "sqrt(s) = " << sqrts << "\n\n";

  // Replicate NLO_Base::CalculateReal's exact call sequence (NLO_Base.C:462)
  // and report the dipole counts after each call. MakeDipolesFF() clears
  // m_dipolesFF and then calls Dipole_FF(), which only sorts particles into
  // charged/neutral lists - it never constructs a Dipole. The function that
  // builds FF dipoles is MakeDipoles().
  std::cout << "--- dipole counts through NLO_Base::CalculateReal's sequence ---\n";
  // m_flav_label is reported too: the FSR recoil loop indexes it at
  // NLO_Base.C:482 (p[m_flav_label[f]] = ...), so an empty map silently writes
  // every recoiled momentum into slot 0, i.e. into a beam.
  auto counts = [&](const char *what) {
    std::cout << std::setw(34) << what
              << "   II=" << dip.m_dipolesII.size()
              << "  FF=" << dip.m_dipolesFF.size()
              << "  IF=" << dip.m_dipolesIF.size()
              << "  flav_label=" << dip.m_flav_label.size() << "\n";
  };
  std::cout << "  [old sequence: MakeDipoles, MakeDipolesIF, MakeDipolesFF]\n";
  dip.MakeDipoles(fl, mom, mom);    counts("after MakeDipoles");
  dip.MakeDipolesIF(fl, mom, mom);  counts("after MakeDipolesIF");
  dip.MakeDipolesFF(fl, mom, mom);  counts("after MakeDipolesFF");
  std::cout << "\n  [MakeDipoles last: MakeDipolesII, MakeDipolesIF, MakeDipoles]\n";
  dip.MakeDipolesII(fl, mom, mom);  counts("after MakeDipolesII");
  dip.MakeDipolesIF(fl, mom, mom);  counts("after MakeDipolesIF");
  dip.MakeDipoles(fl, mom, mom);    counts("after MakeDipoles");
  std::cout << "  (NLO_Base.C:471 then iterates GetDipoleFF() for the FSR recoil,\n"
            << "   and line 482 indexes m_flav_label)\n\n";

  // Rebuild in the order used just before CalculateRealSub (NLO_Base.C:497),
  // which ends with MakeDipoles and so does leave the FF dipole in place.
  dip.MakeDipolesII(fl, mom, mom);
  dip.MakeDipolesIF(fl, mom, mom);
  dip.MakeDipoles(fl, mom, mom);
  std::cout << "--- rebuilt as before CalculateRealSub (NLO_Base.C:497-499) ---\n";
  std::cout << "dipoles built: II=" << dip.m_dipolesII.size()
            << "  FF=" << dip.m_dipolesFF.size()
            << "  IF=" << dip.m_dipolesIF.size() << "\n";
  if (dip.m_dipolesII.empty() || dip.m_dipolesFF.empty())
    std::cout << "  !! expected at least one II and one FF dipole\n";
  // The FF dipoles are only summed by CalculateRealSub when flagged resonant.
  for (auto &D : dip.m_dipolesFF)
    std::cout << "  FF dipole IsResonance = " << D.IsResonance() << "\n";
  std::cout << "\n";

  auto line = [&](const Vec4D &k) {
    const double coh_code   = dip.CalculateRealSub(k);
    const double incoh_code = dip.CalculateRealSubEEX(k);
    const double ifi_code   = dip.CalculateRealSubIF(k);
    // independent reference
    const Vec4D J1 = J(-1., +1., pEm, k), J2 = J(-1., -1., pEp, k);
    const Vec4D J3 = J(-1., -1., pMm, k), J4 = J(-1., +1., pMp, k);
    const double pref  = alpha/(4.*M_PI*M_PI);
    const double coh_ref   = pref*MinusSq(J1+J2+J3+J4);
    const double incoh_ref = pref*(MinusSq(J1+J2) + MinusSq(J3+J4));
    return std::vector<double>{coh_code, incoh_code, ifi_code,
                               coh_ref, incoh_ref, coh_ref-incoh_ref};
  };

  // ---- single point ----
  {
    Vec4D k(0.01, 0.01*sin(1.1), 0., 0.01*cos(1.1));
    auto v = line(k);
    std::cout << "--- single photon k = " << k << " ---\n";
    std::cout << "  CalculateRealSub    (coherent, code) = " << v[0] << "\n";
    std::cout << "  reference -|J_i+J_f|^2 * a/4pi^2     = " << v[3]
              << "    ratio = " << (v[3]!=0?v[0]/v[3]:0) << "\n";
    std::cout << "  CalculateRealSubEEX (incoherent,code)= " << v[1] << "\n";
    std::cout << "  reference incoherent                 = " << v[4]
              << "    ratio = " << (v[4]!=0?v[1]/v[4]:0) << "\n";
    std::cout << "  CalculateRealSubIF  (IF eikonals)    = " << v[2] << "\n";
    std::cout << "  reference IFI = coh-incoh            = " << v[5] << "\n\n";
  }

  // ---- angular scan: is the interference in the SHIPPING code odd? ----
  // Reflect the photon through the origin. With the final-state pair back to
  // back this maps the configuration onto itself with ISR/FSR fixed, so the
  // interference must change sign. A failure here means the shipping real
  // subtraction does not carry the asymmetry-generating term correctly.
  std::cout << "--- photon angular scan (energy fixed at 10 MeV) ---\n";
  std::cout << std::setw(8) << "cos(th)"
            << std::setw(20) << "coherent(code)"
            << std::setw(20) << "incoherent(code)"
            << std::setw(20) << "IFI=coh-incoh"
            << std::setw(20) << "IFI(code,SubIF)" << "\n";
  for (int i = -4; i <= 4; ++i) {
    if (i == 0) continue;
    const double c = 0.2*i, s = sqrt(1.-c*c), w = 0.01;
    Vec4D k(w, w*s, 0., w*c);
    auto v = line(k);
    std::cout << std::setw(8) << c
              << std::setw(20) << v[0]
              << std::setw(20) << v[1]
              << std::setw(20) << v[0]-v[1]
              << std::setw(20) << v[2] << "\n";
  }

  // Oddness of the interference. Reflecting the photon ALONE is not a symmetry
  // here - the final-state pair sits at a fixed angle - so the reflection has
  // to act on the pair too: k -> parity(k) together with p3 <-> p4, which for a
  // back-to-back pair is its own parity image. The dipoles must be rebuilt for
  // the swapped final state, which is why this block re-runs MakeDipoles.
  std::cout << "\n--- oddness: photon parity + final-state swap ---\n";
  std::cout << std::setw(8) << "cos(th)"
            << std::setw(22) << "IFI(direct)"
            << std::setw(22) << "IFI(reflected)"
            << std::setw(16) << "sum (->0?)" << "\n";
  {
    // flavours stay put; the mu- takes the mu+'s momentum and vice versa.
    // (Permuting flavours AND momenta together would be the identity.)
    Vec4D_Vector   momSwap{pEm, pEp, pMp, pMm};
    for (int i = 1; i <= 4; ++i) {
      const double c = 0.2*i, s = sqrt(1.-c*c), w = 0.01;
      Vec4D kp(w,  w*s, 0.,  w*c);
      Vec4D km(w, -w*s, 0., -w*c); // parity image of kp

      dip.MakeDipolesII(fl, mom, mom);
      dip.MakeDipolesIF(fl, mom, mom);
      dip.MakeDipoles(fl, mom, mom);
      const double ip = dip.CalculateRealSub(kp) - dip.CalculateRealSubEEX(kp);

      dip.MakeDipolesII(fl, momSwap, momSwap);
      dip.MakeDipolesIF(fl, momSwap, momSwap);
      dip.MakeDipoles(fl, momSwap, momSwap);
      const double im = dip.CalculateRealSub(km) - dip.CalculateRealSubEEX(km);

      std::cout << std::setw(8) << c
                << std::setw(22) << ip
                << std::setw(22) << im
                << std::setw(16) << ip+im << "\n";
    }
    // restore
    dip.MakeDipolesII(fl, mom, mom);
    dip.MakeDipolesIF(fl, mom, mom);
    dip.MakeDipoles(fl, mom, mom);
  }

  } catch (const ATOOLS::Exception &e) {
    std::cerr << "Sherpa exception: " << e << std::endl;
    return 1;
  }
  return 0;
}
