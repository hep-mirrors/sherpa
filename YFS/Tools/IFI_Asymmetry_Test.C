// Probes the angular (forward-backward) asymmetry induced by the IF-dipole
// YFS form factor that Define_Dipoles::FormFactorSum() adds when IFI_Sub=1:
//
//     Y_IF(costh) = sum_{IF dipoles} D.ChargeNorm() * BVR_full(D, omega)
//
// At sqrt(s)=0.7 GeV, e+e- -> mu+mu- is pure s-channel gamma*, so the Born
// angular distribution is exactly 1+cos^2(theta) and AFB(Born) == 0. Any AFB in
// the generator therefore comes entirely from IFI/box physics, which makes this
// term directly comparable against the measured AFB.
//
// Prints Y_IF and exp(Y_IF) at +/-costh, and the AFB it induces on a
// 1+cos^2 Born inside the analysis acceptance, for the omega actually used in
// FormFactorSum() (sqrt(s)/2) and for KKMC's Yint choice (Emin).
//
// Build/run: see build_IFI_Asymmetry_Test.sh in this directory.
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
#include "YFS/Tools/Dipole.H"
#include "YFS/Main/YFS_Form_Factor.H"
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

int main(int argc, char* argv[]) {
#ifdef USING__MPI
  MPI_Init(&argc, &argv);
#endif
  ATOOLS::mpi = new My_MPI();
  ATOOLS::exh = new Terminator_Object_Handler();
  ATOOLS::msg = new Message();
  ATOOLS::rpa = new Run_Parameter();
  // IR_CUTOFF 1e-6 and FULL_FORM 1 match Sherpa.Mu.NLO.yaml
  char override_arg[] = "YFS: {USE_MODEL_ALPHA: 0, KKMC_ANG: 0, IR_CUTOFF: 1e-6, FULL_FORM: 1}";
  char prog_name[] = "ifi_asym";
  char* fake_argv[] = {prog_name, override_arg};
  Settings::InitializeMainSettings(2, fake_argv);
  ATOOLS::ran = new Random(1234);
  ATOOLS::s_loader = new Library_Loader();
  PDF::pdfdefs = new PDF::PDF_Defaults();

  try {

  double mmu = 0.105658375, mel = 0.000511;
  double sqrts = 0.7, s = sqrts*sqrts;
  rpa->gen.SetEcms(sqrts);
  double pbeam = sqrt(sqr(sqrts/2.)-sqr(mel));
  double pmu   = sqrt(sqr(sqrts/2.)-sqr(mmu));

  // icharge is in units of e/3 (so -3 -> charge -1), spin=1 for a fermion.
  s_kftable[11] = new Particle_Info(11, mel, 0., -3, 1, 0, "e-", "e^-");
  s_kftable[13] = new Particle_Info(13, mmu, 0., -3, 1, 0, "mu-", "\\mu^-");
  double alpha = 1./137.035999084;
  MODEL::s_model = new FakeModel(alpha);

  Flavour flEm((kf_code)11, false);   // e-  (beam 0)
  Flavour flEp((kf_code)11, true);    // e+  (beam 1)
  Flavour flMum((kf_code)13, false);  // mu-
  Flavour flMup((kf_code)13, true);   // mu+

  Vec4D pEm( sqrts/2., 0., 0.,  pbeam);
  Vec4D pEp( sqrts/2., 0., 0., -pbeam);

  YFS_Form_Factor form;
  double Emin = 0.5 * s * (1e-6/sqrts);   // FSR::Initialize()'s m_Emin formula

  // Y_IF as a function of cos(theta) of the mu-, for a given omega.
  auto YIF = [&](double costh, double omega) {
    double sinth = sqrt(1.-costh*costh);
    Vec4D pMum( sqrts/2.,  pmu*sinth, 0.,  pmu*costh);
    Vec4D pMup( sqrts/2., -pmu*sinth, 0., -pmu*costh);
    Flavour_Vector fl{flEm, flEp, flMup, flMum};
    Vec4D_Vector   mom{pEm, pEp, pMup, pMum};
    double sum = 0.;
    for (size_t i = 0; i < 2; ++i)
      for (size_t j = 2; j < 4; ++j) {
        Flavour_Vector ff{fl[i], fl[j]};
        Vec4D_Vector   mm{mom[i], mom[j]};
        Dipole D(ff, mm, mm, dipoletype::ifi, alpha);
        D.SetResonance(false);
        // force-correct in case Flavour::Mass()/Charge() came back 0 in this
        // minimal (no MODEL) setup
        D.m_masses  = {fl[i].Mass(), fl[j].Mass()};
        D.m_charges = {fl[i].Charge(), fl[j].Charge()};
        sum += D.ChargeNorm() * form.BVR_full(D, omega);
      }
    return sum;
  };

  std::cout << std::setprecision(6);
  struct Choice { const char* name; double omega; };
  Choice choices[2] = {
    {"omega = sqrt(s)/2  (what FormFactorSum() uses)", sqrts/2.},
    {"omega = Emin       (what KKMC's Yint uses)    ", Emin},
  };
  std::cout << "Emin = " << Emin << " GeV,  sqrt(s)/2 = " << sqrts/2. << " GeV\n";

  for (auto &ch : choices) {
    std::cout << "\n=== " << ch.name << " ===\n";
    std::cout << "  costh      Y_IF(+c)     Y_IF(-c)    exp(Y(+c))  exp(Y(-c))   ratio\n";
    // AFB induced on a 1+cos^2 Born over the analysis acceptance |costh|<0.55
    double F = 0., B = 0.;
    const int N = 220;
    for (int i = 0; i < N; ++i) {
      double c = -0.55 + (i+0.5)*1.10/N;
      double born = 1. + c*c;
      double w = exp(YIF(c, ch.omega));
      if (c > 0) F += born*w; else B += born*w;
    }
    for (double c : {0.1, 0.2, 0.3, 0.4, 0.5}) {
      double yp = YIF(c, ch.omega), ym = YIF(-c, ch.omega);
      std::cout << "  " << std::setw(5) << c
                << std::setw(13) << yp << std::setw(13) << ym
                << std::setw(13) << exp(yp) << std::setw(12) << exp(ym)
                << std::setw(10) << exp(yp)/exp(ym) << "\n";
    }
    std::cout << "  -> AFB induced on (1+cos^2) Born, |costh|<0.55 : "
              << std::showpos << (F-B)/(F+B) << std::noshowpos << "\n";
  }
  std::cout << "\n(measured: KKMC -0.01015 / BabaYaga -0.01104 / Sherpa -0.01605 for AFB(mu+),"
            << "\n so AFB(mu-) = +those; Sherpa's excess is ~ +0.005 in |AFB|)\n";

  // ---- NLO-only piece: the exponentiated IF virtual (BVV_full, photon-mass
  // regularised, what BVR_full(D,omega) uses at FULL_FORM=1) vs the IF virtual
  // that CalculateVirtualSubEps() subtracts from the OpenLoops one-loop
  // (BVV_full_eps, dim-reg). Any t/u-dependent difference between the two is a
  // spurious asymmetry that appears ONLY once NLO is switched on, because the
  // subtraction is only called then. ----
  auto IFvirt = [&](double costh, int which) {
    double sinth = sqrt(1.-costh*costh);
    Vec4D pMum( sqrts/2.,  pmu*sinth, 0.,  pmu*costh);
    Vec4D pMup( sqrts/2., -pmu*sinth, 0., -pmu*costh);
    Flavour_Vector fl{flEm, flEp, flMup, flMum};
    Vec4D_Vector   mom{pEm, pEp, pMup, pMum};
    double sum = 0.;
    for (size_t i = 0; i < 2; ++i)
      for (size_t j = 2; j < 4; ++j) {
        Flavour_Vector ff{fl[i], fl[j]};
        Vec4D_Vector   mm{mom[i], mom[j]};
        Dipole D(ff, mm, mm, dipoletype::ifi, alpha);
        D.SetResonance(false);
        D.m_masses  = {fl[i].Mass(), fl[j].Mass()};
        D.m_charges = {fl[i].Charge(), fl[j].Charge()};
        double v = (which == 0)
          ? form.BVV_full(D, form.m_photonMass, sqrts/2., 0)   // in exp(Y)
          : form.BVV_full_eps(D, sqrts/2., 3).Finite();        // subtracted at NLO
        sum += D.ChargeNorm() * v;
      }
    return sum;
  };

  std::cout << "\n=== IF virtual: exponentiated (BVV_full) vs subtracted (BVV_full_eps) ===\n";
  std::cout << "  costh    exp'd(+c)    sub'd(+c)     diff(+c)     diff(-c)\n";
  for (double c : {0.1, 0.2, 0.3, 0.4, 0.5}) {
    double a = IFvirt(c,0), b = IFvirt(c,1);
    double am = IFvirt(-c,0), bm = IFvirt(-c,1);
    std::cout << "  " << std::setw(5) << c << std::setw(13) << a << std::setw(13) << b
              << std::setw(13) << (a-b) << std::setw(13) << (am-bm) << "\n";
  }
  {
    double F=0., B=0.; const int N=220;
    for (int i=0;i<N;++i) {
      double c = -0.55 + (i+0.5)*1.10/N;
      double born = 1.+c*c;
      double w = 1. + (IFvirt(c,0) - IFvirt(c,1));  // O(alpha) residual
      if (c>0) F += born*w; else B += born*w;
    }
    std::cout << "  -> AFB induced by the exp'd-minus-subtracted IF virtual mismatch: "
              << std::showpos << (F-B)/(F+B) << std::noshowpos << "\n";
  }

  } catch (const ATOOLS::Exception& e) {
    std::cerr << "ATOOLS::Exception: " << e << std::endl;
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "std::exception: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
