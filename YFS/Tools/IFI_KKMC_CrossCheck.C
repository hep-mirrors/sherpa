// Standalone comparison harness for the INITIAL-FINAL INTERFERENCE form
// factor, in the same spirit as FSR_KKMC_CrossCheck.C: drives the real
// YFS::Dipole / YFS_Form_Factor code with a fixed, hand-picked phase-space
// point and prints every intermediate so it can be compared number-for-number
// against KKMC on identical input.
//
// What is being compared
// ----------------------
// KKMC assembles the interference form factor as a product over the four
// initial-final fermion pairs (SRCee/KKceex.cxx:334):
//
//   Yint = TForFac( +alfpmix, p1,p3, Emin) * TForFac( +alfpmix, p2,p4, Emin)
//        * TForFac( -alfpmix, p1,p4, Emin) * TForFac( -alfpmix, p2,p3, Emin)
//
// with (SRCee/KKceex.cxx:3944)
//
//   TForFac(a,p1,p2,Emin,MasPhot) = exp( Btilda(a,...,Emin,MasPhot)
//                                      + TBvirt(a,...,MasPhot) )
//   alfpmix = alfpi * ChaIni * ChaFin
//
// Sherpa builds the same object in Define_Dipoles::FormFactorSum() as
//
//   Y_IF = sum over m_dipolesIF of ChargeNorm(D) * IFForFac(D, omega)
//
// with IFForFac = BVR_full(...) + BVirtT(...) (YFS_Form_Factor.C:995), i.e.
// the same Btilda + t-channel-virtual split, and ChargeNorm() = -QiQj*thetaij
// carrying KKMC's +/- pattern. So the two should satisfy
//
//   exp( Y_IF )  ==  Yint      pair by pair, and in the product.
//
// Dipole_IF (Define_Dipoles.C:438) builds the four pairs in the order
// (e-,mu-), (e-,mu+), (e+,mu-), (e+,mu+), i.e. KKMC's (p1,p3), (p1,p4),
// (p2,p3), (p2,p4).
//
// The omega scan at the end is the point of the exercise: Y_IF depends on the
// soft cutoff, and that dependence is only physical if it cancels against the
// real interference integrated over the generated photons. The scan quantifies
// how big that uncancelled dependence is.
//
// Build/run: ./build_IFI_KKMC_CrossCheck.sh
// The matching KKMC-side driver is Test/SherpaCompare/kkmc_ifi_crosscheck.cxx.

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
#include "YFS/Main/Dipole.H"
#include "YFS/Main/YFS_Form_Factor.H"
#include "MODEL/Main/Model_Base.H"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace ATOOLS;
using namespace YFS;

// Same minimal stand-in as FSR_KKMC_CrossCheck.C: YFS_Base::RegisterSettings()
// queries s_model->ScalarConstant("alpha_QED") on one branch regardless.
class FakeModel : public MODEL::Model_Base {
public:
  FakeModel(double alpha) : MODEL::Model_Base(true) {
    (*p_constants)["alpha_QED"] = alpha;
  }
  void InitVertices() override {}
  void ParticleInit() override {}
  bool ModelInit() override { return true; }
};

int main(int argc, char* argv[]) {
#ifdef USING__MPI
  MPI_Init(&argc, &argv);
#endif
  ATOOLS::mpi = new My_MPI();
  ATOOLS::exh = new Terminator_Object_Handler();
  ATOOLS::msg = new Message();
  ATOOLS::rpa = new Run_Parameter();

  char override_arg[] = "YFS: {USE_MODEL_ALPHA: 0, KKMC_ANG: 1, IR_CUTOFF: 1e-9}";
  char prog_name[] = "sherpa_ifi_driver";
  char* fake_argv[] = {prog_name, override_arg};
  Settings::InitializeMainSettings(2, fake_argv);
  ATOOLS::ran = new Random(1234);
  ATOOLS::s_loader = new Library_Loader();
  PDF::pdfdefs = new PDF::PDF_Defaults();

  try {

  // ---- fixed phase-space point: e+e- -> mu+ mu- at sqrt(s)=0.7 GeV ----
  // Matches the muon runcard (Sherpa.Mu.NLO.yaml): BEAM_ENERGIES 0.7/2.
  const double me    = 0.000511;
  const double mmu   = 0.105658375;
  const double sqrts = 0.7;
  const double pin   = sqrt(sqr(sqrts/2.) - sqr(me));
  const double pout  = sqrt(sqr(sqrts/2.) - sqr(mmu));
  // 40 degrees, so the interference (which is odd in cos(theta)) is not
  // accidentally sitting at a zero.
  const double ct = cos(0.7), st = sin(0.7);

  rpa->gen.SetEcms(sqrts);

  Vec4D pEm(sqrts/2., 0., 0.,  pin);   // p1  incoming e-
  Vec4D pEp(sqrts/2., 0., 0., -pin);   // p2  incoming e+
  Vec4D pMm(sqrts/2.,  pout*st, 0.,  pout*ct);  // p3  outgoing mu-
  Vec4D pMp(sqrts/2., -pout*st, 0., -pout*ct);  // p4  outgoing mu+

  s_kftable[11] = new Particle_Info(11, me,  0., -3, 0, 1, "e-",  "e^-");
  s_kftable[13] = new Particle_Info(13, mmu, 0., -3, 0, 1, "mu-", "\\mu^-");
  const double alpha = 1./137.035999084;
  MODEL::s_model = new FakeModel(alpha);

  Flavour flEm((kf_code)11, false), flEp((kf_code)11, true);
  Flavour flMm((kf_code)13, false), flMp((kf_code)13, true);

  std::cout << std::setprecision(15);
  std::cout << "=== IFI form-factor cross-check, e+e- -> mu+mu- ===\n";
  std::cout << "sqrt(s)   = " << sqrts << "\n";
  std::cout << "alpha     = " << alpha << "\n";
  std::cout << "p1 (e- )  = " << pEm << "\n";
  std::cout << "p2 (e+ )  = " << pEp << "\n";
  std::cout << "p3 (mu-)  = " << pMm << "\n";
  std::cout << "p4 (mu+)  = " << pMp << "\n\n";

  // ---- the four initial-final pairs, in Dipole_IF's order ----
  struct Pair { const char *name; Flavour fi, fj; Vec4D pi, pj; };
  std::vector<Pair> pairs = {
    {"(e- ,mu-) = KKMC (p1,p3)", flEm, flMm, pEm, pMm},
    {"(e- ,mu+) = KKMC (p1,p4)", flEm, flMp, pEm, pMp},
    {"(e+ ,mu-) = KKMC (p2,p3)", flEp, flMm, pEp, pMm},
    {"(e+ ,mu+) = KKMC (p2,p4)", flEp, flMp, pEp, pMp},
  };

  YFS_Form_Factor ff;

  // omega = sqrt(s)/2 is what FormFactorSum() passes (Define_Dipoles.C:652).
  auto report = [&](double omega, bool verbose) {
    double sum(0.);
    if (verbose) {
      std::cout << "--- per-pair breakdown at omega = " << omega << " ---\n";
      std::cout << std::setw(28) << "pair"
                << std::setw(14) << "ChargeNorm"
                << std::setw(20) << "Breal(Btilda)"
                << std::setw(20) << "Bvirt(TBvirt)"
                << std::setw(20) << "IFForFac"
                << std::setw(20) << "ChargeNorm*IFF" << "\n";
    }
    for (auto &pr : pairs) {
      Flavour_Vector fl{pr.fi, pr.fj};
      Vec4D_Vector   mm{pr.pi, pr.pj};
      Vec4D_Vector   bm{pr.pi, pr.pj};
      Dipole D(fl, mm, bm, dipoletype::ifi, alpha);
      D.SetResonance(false);

      // IFForFac's momentum order (GetBornMomenta(1) then (0)) is internal to
      // it; reproduce the pieces here only for printing.
      Vec4D q1 = D.GetBornMomenta(1), q2 = D.GetBornMomenta(0);
      const double breal = ff.BVR_full(q1*q2, q1.E(), q2.E(), q1.Mass(), q2.Mass(),
                                       omega, ff.m_photonMass, 0);
      const double bvirt = ff.BVirtT(q1, q2, omega);
      const double iff   = ff.IFForFac(D, omega);
      const double cn    = D.ChargeNorm();
      sum += cn * iff;
      if (verbose) {
        std::cout << std::setw(28) << pr.name
                  << std::setw(14) << cn
                  << std::setw(20) << breal
                  << std::setw(20) << bvirt
                  << std::setw(20) << iff
                  << std::setw(20) << cn*iff << "\n";
      }
    }
    if (verbose) {
      std::cout << "\n  Y_IF = sum ChargeNorm*IFForFac = " << sum << "\n";
      std::cout << "  exp(Y_IF)  [ == KKMC Yint ]     = " << exp(sum) << "\n\n";
    }
    return sum;
  };

  const double omega_prod = sqrts/2.;
  report(omega_prod, true);

  // ---- omega dependence ----
  // Y_IF is a soft/virtual object: its cutoff dependence is only physical if
  // the real interference integrated up to the same cutoff cancels it. This
  // scan is the size of what has to cancel.
  std::cout << "--- omega dependence of Y_IF and exp(Y_IF) ---\n";
  std::cout << std::setw(18) << "omega"
            << std::setw(22) << "Y_IF"
            << std::setw(22) << "exp(Y_IF)"
            << std::setw(22) << "exp(Y)/exp(Y_prod)" << "\n";
  const double yprod = report(omega_prod, false);
  for (double om : {sqrts/2., sqrts/4., 0.05, 0.01, 1e-3, 1e-4}) {
    const double y = report(om, false);
    std::cout << std::setw(18) << om
              << std::setw(22) << y
              << std::setw(22) << exp(y)
              << std::setw(22) << exp(y - yprod) << "\n";
  }
  std::cout << "\nIf Y_IF is to be cancelled by the real interference, the\n"
            << "log-slope below must match d/dlog(omega) of the integrated\n"
            << "real IF eikonal over the generated photons:\n";
  {
    const double y1 = report(0.01, false), y2 = report(0.02, false);
    std::cout << "  dY_IF/dlog(omega) ~ " << (y2-y1)/log(2.) << "\n";
  }

  // ---- s-channel (II and FF) form factors ----
  // FormFactorSum() gives the II and FF dipoles BVR_full(D, omega), which for
  // the default FULL_FORM=1 is BVR_full(...) + BVV_full(...), i.e. the same
  // Btilda + s-channel-virtual split KKMC uses in SForFac (KKceex.cxx:3927):
  //     SForFac = exp( Btilda(...,Emin,MasPhot) + SBvirt(...,MasPhot) )
  // These dominate the rate, unlike the IF term, and had never been compared.
  std::cout << "\n--- s-channel form factors (compare against KKMC SForFac) ---\n";
  std::cout << std::setw(18) << "dipole"
            << std::setw(22) << "Breal(Btilda)"
            << std::setw(22) << "Bvirt(SBvirt)"
            << std::setw(22) << "sum" << "\n";
  {
    struct SPair { const char *name; Vec4D a, b; };
    // II uses omega = sqrt(s)/2 (Define_Dipoles.C:603); FF uses
    // sqrt(D.Sprime())/2, which for the Born configuration is also sqrt(s)/2.
    std::vector<SPair> sp = {{"II (e-,e+)", pEm, pEp}, {"FF (mu-,mu+)", pMm, pMp}};
    for (auto &pr : sp) {
      const double p1p2 = pr.a * pr.b;
      const double breal = ff.BVR_full(p1p2, pr.a.E(), pr.b.E(),
                                       pr.a.Mass(), pr.b.Mass(),
                                       omega_prod, ff.m_photonMass, 0);
      const double bvirt = ff.BVV_full(pr.a, pr.b, ff.m_photonMass, omega_prod, 0);
      std::cout << std::setw(18) << pr.name
                << std::setw(22) << breal
                << std::setw(22) << bvirt
                << std::setw(22) << breal+bvirt << "\n";
    }
  }

  // ---- angular scan ----
  // Y_IF is the object that generates the forward-backward asymmetry, so what
  // matters for A_FB is not its value at one angle but its dependence on
  // cos(theta) - and in particular that it is ODD, i.e. Y_IF(-c) = -Y_IF(c)
  // once the pair is reflected. A constant offset shifts sigma; only the
  // angular shape moves A_FB. Emitted in a parseable form so it can be diffed
  // directly against the KKMC driver's identical scan.
  // How big is the ODD part of Y_IF as omega is lowered from sqrt(s)/2 (what
  // FormFactorSum passes) towards the actual generation cutoff? That
  // difference is the piece the real interference would have to cancel, and
  // it is the natural scale of any residual A_FB discrepancy.
  std::cout << "\n--- odd part of Y_IF at cos(theta)=0.5 vs omega ---\n";
  std::cout << std::setw(14) << "omega" << std::setw(24) << "Y_IF(c=0.5)"
            << std::setw(24) << "exp(Y_IF)-1 [%]" << "\n";
  for (double om : {sqrts/2., 0.05, 0.01, 1e-3, 1e-5, 1e-7, 1e-9}) {
    const double c = 0.5, s = sqrt(1.-c*c);
    Vec4D q3(sqrts/2.,  pout*s, 0.,  pout*c);
    Vec4D q4(sqrts/2., -pout*s, 0., -pout*c);
    std::vector<Pair> sc = {
      {"13", flEm, flMm, pEm, q3}, {"14", flEm, flMp, pEm, q4},
      {"23", flEp, flMm, pEp, q3}, {"24", flEp, flMp, pEp, q4},
    };
    double sum(0.);
    for (auto &pr : sc) {
      Flavour_Vector fl{pr.fi, pr.fj};
      Vec4D_Vector mm{pr.pi, pr.pj}, bm{pr.pi, pr.pj};
      Dipole D(fl, mm, bm, dipoletype::ifi, alpha);
      D.SetResonance(false);
      sum += D.ChargeNorm()*ff.IFForFac(D, om);
    }
    std::cout << std::setw(14) << om << std::setw(24) << sum
              << std::setw(24) << 100.*(exp(sum)-1.) << "\n";
  }

  std::cout << "\n#SCAN costheta Y_IF exp(Y_IF)\n";
  for (int i = -9; i <= 9; ++i) {
    const double c = 0.1*i;
    const double s = sqrt(1.-c*c);
    Vec4D q3(sqrts/2.,  pout*s, 0.,  pout*c);
    Vec4D q4(sqrts/2., -pout*s, 0., -pout*c);
    std::vector<Pair> sc = {
      {"13", flEm, flMm, pEm, q3}, {"14", flEm, flMp, pEm, q4},
      {"23", flEp, flMm, pEp, q3}, {"24", flEp, flMp, pEp, q4},
    };
    double sum(0.);
    for (auto &pr : sc) {
      Flavour_Vector fl{pr.fi, pr.fj};
      Vec4D_Vector mm{pr.pi, pr.pj}, bm{pr.pi, pr.pj};
      Dipole D(fl, mm, bm, dipoletype::ifi, alpha);
      D.SetResonance(false);
      sum += D.ChargeNorm()*ff.IFForFac(D, omega_prod);
    }
    std::cout << "#SCAN " << std::setw(6) << c
              << std::setw(24) << sum
              << std::setw(24) << exp(sum) << "\n";
  }

  } catch (const ATOOLS::Exception &e) {
    std::cerr << "Sherpa exception: " << e << std::endl;
    return 1;
  }
  return 0;
}
