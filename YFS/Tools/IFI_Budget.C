// A_FB budget for the initial-final interference, e+e- -> mu+mu- at 0.7 GeV.
//
// IFI_KKMC_CrossCheck.C established that Sherpa's per-pair IF form factor
// (Btilda + TBvirt) agrees with KKMC's TForFac number-for-number. This tool
// answers the next question: given that the ingredients agree, what does the
// ASSEMBLY in Define_Dipoles::FormFactorSum() actually produce for A_FB, and
// how far is each candidate assembly from KKMC's measured -0.0103?
//
// Every variant below is evaluated by driving the real YFS::Dipole /
// YFS_Form_Factor code, so any convention bug in them shows up here too.
//
//   WAS       exp( sum -ChargeNorm * IFForFac(D, sqrt(D.Sprime())/2) )
//             what FormFactorSum() ran before this was diagnosed
//   WAS,+sign exp( sum +ChargeNorm * IFForFac(D, sqrt(D.Sprime())/2) )
//   NOW       exp( sum +ChargeNorm * IFForFac(D, IFIOmega()) ), IFIOmega()
//             being sqrt(s)/2 with IFI_Real off. This is the coefficient
//             CalculateVirtualSub()/TFormFactor() use on the same dipoles.
//   YINT      exp( sum +ChargeNorm * IFForFac(D, Emin) )  == KKMC's Yint,
//             which in KKMC is compensated by the real interference summed
//             over generated photons; uncompensated here
//   BOX       the finite (IR-subtracted) gamma-gamma box KKMC adds on top of
//             Yint (KKbvir::CBoxGG - KKbvir::IntIR, KKceex.cxx:800-812).
//             Sherpa's YFS has no equivalent term at all.
//
// Two things about the WAS line: the sign was opposite to every other
// IF-dipole consumer in Define_Dipoles.C, and omega = sqrt(Sprime)/2 is an
// ANGLE-DEPENDENT cutoff, since for an initial-final pair
// Sprime = (p_in + p_out)^2 ~ (s/2)(1 -+ beta cos(theta)). A soft cutoff that
// is itself odd in cos(theta) injects asymmetry by construction, and that is
// where essentially all of Sherpa's A_FB was coming from.
//
// Build/run: ./build_IFI_Budget.sh

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
#include <functional>

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
  const double me    = 0.000511;
  const double mmu   = 0.105658375;
  const double sqrts = 0.7;
  const double svar  = sqrts*sqrts;
  const double alpha = 1./137.035999084;
  double pin, pout, beta;
}

// KKMC KKbvir::CBoxGG and KKbvir::IntIR, real parts, IR-subtracted exactly as
// KKceex::BornFoam2 does. The photon-mass logs cancel between them, so what is
// left is finite and t<->u antisymmetric, i.e. pure A_FB.
// Returns {2*Re BoxGGtu, 2*Re BoxGGut} / (alpha/pi) already multiplied by
// Coef = (alpha/pi)*Qe*Qf = +alpha/pi for e -> mu.
static void ReBox(double s, double t, double u, double &tu, double &ut) {
  const double alpi = alpha/M_PI;
  const double L  = log(t/u);        // t,u < 0 so the ratio is positive
  const double lt = log(-t/s), lu = log(-u/s);
  const double ir = L*log(sqrt(t*u)/s) - 0.5*L;   // CBoxGG_IR - IntIR
  tu = alpi*( ir + 0.25*s*(u-t)/(u*u)*lt*lt - 0.5*(s/u)*lt );
  ut = alpi*( ir - 0.25*s*(t-u)/(t*t)*lu*lu + 0.5*(s/t)*lu );
}

int main(int argc, char* argv[]) {
#ifdef USING__MPI
  MPI_Init(&argc, &argv);
#endif
  ATOOLS::mpi = new My_MPI();
  ATOOLS::exh = new Terminator_Object_Handler();
  ATOOLS::msg = new Message();
  ATOOLS::rpa = new Run_Parameter();

  // Same YFS settings as the production muon card (FULL_FORM 1, IR_CUTOFF 1e-6).
  char override_arg[] = "YFS: {USE_MODEL_ALPHA: 0, KKMC_ANG: 0, IR_CUTOFF: 1e-6, FULL_FORM: 1}";
  char prog_name[] = "ifi_budget";
  char* fake_argv[] = {prog_name, override_arg};
  Settings::InitializeMainSettings(2, fake_argv);
  ATOOLS::ran = new Random(1234);
  ATOOLS::s_loader = new Library_Loader();
  PDF::pdfdefs = new PDF::PDF_Defaults();

  try {

  pin  = sqrt(sqr(sqrts/2.) - sqr(me));
  pout = sqrt(sqr(sqrts/2.) - sqr(mmu));
  beta = pout/(sqrts/2.);
  rpa->gen.SetEcms(sqrts);

  s_kftable[11] = new Particle_Info(11, me,  0., -3, 1, 0, "e-",  "e^-");
  s_kftable[13] = new Particle_Info(13, mmu, 0., -3, 1, 0, "mu-", "\\mu^-");
  MODEL::s_model = new FakeModel(alpha);

  Flavour flEm((kf_code)11, false), flEp((kf_code)11, true);
  Flavour flMm((kf_code)13, false), flMp((kf_code)13, true);

  const Vec4D pEm(sqrts/2., 0., 0.,  pin);
  const Vec4D pEp(sqrts/2., 0., 0., -pin);

  YFS_Form_Factor ff;
  const double Emin = 0.5*svar*(1e-6/sqrts);   // FSR::Initialize()'s m_Emin

  // Build the four IF dipoles for a given cos(theta) of the mu-, in exactly
  // the order Define_Dipoles::Dipole_IF() does: (e-,mu-) (e-,mu+) (e+,mu-) (e+,mu+).
  auto dipoles = [&](double c) {
    const double s_ = sqrt(1.-c*c);
    const Vec4D pMm(sqrts/2.,  pout*s_, 0.,  pout*c);
    const Vec4D pMp(sqrts/2., -pout*s_, 0., -pout*c);
    std::vector<Dipole> D;
    const Flavour fi[2] = {flEm, flEp};
    const Vec4D   pi[2] = {pEm, pEp};
    const Flavour fj[2] = {flMm, flMp};
    const Vec4D   pj[2] = {pMm, pMp};
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j) {
        Flavour_Vector fl{fi[i], fj[j]};
        Vec4D_Vector   mm{pi[i], pj[j]};
        Dipole d(fl, mm, mm, dipoletype::ifi, alpha);
        d.SetResonance(false);
        d.m_masses  = {fi[i].Mass(), fj[j].Mass()};
        d.m_charges = {fi[i].Charge(), fj[j].Charge()};
        D.push_back(d);
      }
    return D;
  };

  // --- the exponents ---------------------------------------------------
  auto Y_old = [&](double c) {             // what Define_Dipoles.C:652 used to run
    double sum(0.); for (auto &d : dipoles(c))
      sum += -d.ChargeNorm()*ff.IFForFac(d, sqrt(d.Sprime())/2.);
    return sum;
  };
  auto Y_oldsign = [&](double c) {         // same omega, sign as documented
    double sum(0.); for (auto &d : dipoles(c))
      sum +=  d.ChargeNorm()*ff.IFForFac(d, sqrt(d.Sprime())/2.);
    return sum;
  };
  // +ChargeNorm and one shared cutoff = what FormFactorSum() runs now, with
  // IFIOmega() returning sqrt(s)/2 whenever IFI_Real is off.
  auto Y_kkmcsign = [&](double c) {
    double sum(0.); for (auto &d : dipoles(c))
      sum +=  d.ChargeNorm()*ff.IFForFac(d, sqrts/2.);
    return sum;
  };
  auto Y_yint = [&](double c) {            // KKMC's Yint, uncompensated
    double sum(0.); for (auto &d : dipoles(c))
      sum +=  d.ChargeNorm()*ff.IFForFac(d, Emin);
    return sum;
  };
  auto Y_orig = [&](double c) {            // the formula before the IFI commits
    double sum(0.); for (auto &d : dipoles(c))
      sum += 0.5*d.ChargeNorm()*ff.BVR_full(d, sqrt(d.Sprime())/2.);
    return sum;
  };
  // +ChargeNorm and a COMMON omega, i.e. KKMC's assembly with a free cutoff.
  auto Y_omega = [&](double c, double om) {
    double sum(0.); for (auto &d : dipoles(c))
      sum +=  d.ChargeNorm()*ff.IFForFac(d, om);
    return sum;
  };
  auto Y_box = [&](double c) {             // 2*Re(box), weighted by helicity classes
    const double t = sqr(mmu) - 0.5*svar*(1.-beta*c);
    const double u = sqr(mmu) - 0.5*svar*(1.+beta*c);
    double btu, but; ReBox(svar, t, u, btu, but);
    return (u*u*2.*btu + t*t*2.*but)/(t*t+u*u);
  };

  // A_FB(mu+) induced on the Born (2-beta^2) + beta^2 c^2 inside |c| < cmax.
  // cos(theta) throughout is the mu- direction, so A_FB(mu+) = -A_FB(mu-).
  auto AFB = [&](const std::function<double(double)> &Y, double cmax, bool expo) {
    double F(0.), B(0.);
    const int N = 2000;
    for (int i = 0; i < N; ++i) {
      const double c = -cmax + (i+0.5)*2.*cmax/N;
      const double born = (2.-beta*beta) + beta*beta*c*c;
      const double w = expo ? exp(Y(c)) : 1.+Y(c);
      if (c > 0) F += born*w; else B += born*w;
    }
    return -(F-B)/(F+B);          // minus: quoted for the mu+
  };

  std::cout << std::setprecision(6);
  std::cout << "=== IFI A_FB budget, e+e- -> mu+mu-, sqrt(s) = " << sqrts << " GeV ===\n";
  std::cout << "beta_mu = " << beta << ",  Emin = " << Emin << " GeV\n\n";

  // --- the angle-dependent cutoff the live formula uses -----------------
  std::cout << "--- omega = sqrt(D.Sprime())/2 per IF dipole (the LIVE choice) ---\n";
  std::cout << std::setw(8) << "cos" << std::setw(14) << "(e-,mu-)"
            << std::setw(14) << "(e-,mu+)" << std::setw(14) << "(e+,mu-)"
            << std::setw(14) << "(e+,mu+)" << "     [sqrt(s)/2 = " << sqrts/2. << "]\n";
  for (double c : {-0.5, -0.25, 0., 0.25, 0.5}) {
    std::cout << std::setw(8) << c;
    for (auto &d : dipoles(c)) std::cout << std::setw(14) << sqrt(d.Sprime())/2.;
    std::cout << "\n";
  }

  // --- how the two halves of IFForFac balance ---------------------------
  std::cout << "\n--- Breal vs Bvirt in sum ChargeNorm*IFForFac (omega = sqrt(s)/2) ---\n";
  std::cout << std::setw(8) << "cos" << std::setw(16) << "sum CN*Breal"
            << std::setw(16) << "sum CN*Bvirt" << std::setw(16) << "sum (odd part)" << "\n";
  for (double c : {0.1, 0.3, 0.5, 0.7, 0.9}) {
    double br(0.), bv(0.);
    for (auto &d : dipoles(c)) {
      const Vec4D q1 = d.GetBornMomenta(1), q2 = d.GetBornMomenta(0);
      br += d.ChargeNorm()*ff.BVR_full(q1*q2, q1.E(), q2.E(), q1.Mass(), q2.Mass(),
                                       sqrts/2., ff.m_photonMass, 0);
      bv += d.ChargeNorm()*ff.BVirtT(q1, q2, sqrts/2.);
    }
    std::cout << std::setw(8) << c << std::setw(16) << br
              << std::setw(16) << bv << std::setw(16) << br+bv << "\n";
  }
  std::cout << "\n  The two are each ~20x the sum: A_FB lives on a cancellation between\n"
            << "  Breal and Bvirt, so any mismatch in sign or in omega between them is\n"
            << "  amplified by that factor.\n";

  // --- the budget -------------------------------------------------------
  struct Row { const char *name; std::function<double(double)> Y; bool expo; };
  std::vector<Row> rows = {
    {"ORIG     +0.5CN * BVR_full(sqrt(Sprime)/2)", Y_orig,     true},
    {"WAS      -CN * IFF(sqrt(Sprime)/2)        ", Y_old,      true},
    {"WAS,+sign+CN * IFF(sqrt(Sprime)/2)        ", Y_oldsign,  true},
    {"NOW      +CN * IFF(IFIOmega() = sqrt(s)/2)", Y_kkmcsign, true},
    {"YINT     +CN * IFF(Emin)                  ", Y_yint,     true},
    {"BOX      finite gamma-gamma box           ", Y_box,      false},
  };
  for (double cmax : {0.55, 0.9}) {
    std::cout << "\n--- A_FB(mu+) induced inside |cos(theta)| < " << cmax << " ---\n";
    for (auto &r : rows)
      std::cout << "  " << r.name << " : " << std::showpos
                << std::setw(12) << AFB(r.Y, cmax, r.expo) << std::noshowpos << "\n";
    std::cout << "  " << "KKMCSIGN + BOX                            " << " : " << std::showpos
              << std::setw(12)
              << AFB([&](double c){ return Y_kkmcsign(c)+Y_box(c); }, cmax, true)
              << std::noshowpos << "\n";
    std::cout << "  " << "------------------------------------------  ------------\n";
    std::cout << "  " << "KKMC CEEX2 (measured)                      " << " :    -0.010300\n";
    std::cout << "  " << "BabaYaga   (measured)                      " << " :    -0.011040\n";
    std::cout << "  " << "Sherpa     (measured, IFI_Asymmetry_Test.C)" << " :    -0.016050\n";
  }

  // --- is there a common omega that reproduces KKMC? --------------------
  // With +ChargeNorm and one cutoff shared by all four pairs -- KKMC's actual
  // assembly -- A_FB is a straight line in log(omega), because Btilda's only
  // cutoff dependence is (p1p2*A - 1)*log(4 omega^2 / MasPhot^2). Reproducing
  // -0.0103 then just picks an omega off that line; it is a fit, not a
  // prediction, and the number it picks is not a scale in the problem.
  std::cout << "\n--- A_FB(mu+) vs a COMMON cutoff omega, |cos| < 0.55 ---\n";
  std::cout << std::setw(14) << "omega [GeV]" << std::setw(16) << "A_FB(mu+)" << "\n";
  for (double om : {sqrts/2., 0.2, 0.1, 0.056, 0.05, 0.01, 1e-3, 1e-5, Emin})
    std::cout << std::setw(14) << om << std::setw(16) << std::showpos
              << AFB([&](double c){ return Y_omega(c, om); }, 0.55, true)
              << std::noshowpos << "\n";

  // --- where the live formula's asymmetry comes from --------------------
  // --- would the real interference actually cancel the cutoff? ----------
  // Define_Dipoles::CalculateRealSubIF(k) already computes sum_IF D.Eikonal(k),
  // the initial-final part of the radiation function, but nothing consumes it:
  // its only caller is a histogram fill (NLO_Base.C:609). It never reaches the
  // event weight, which is why omega has nothing to cancel against.
  //
  // Integrating it over photon phase space between two cutoffs must reproduce
  // the shift in Btilda between the same two cutoffs. If it does, reweighting
  // each generated photon by the IF eikonal makes the omega dependence drop
  // out and the cutoff stops being a tuning knob.
  {
    const double c = 0.5, om1 = 0.01, om2 = 0.05;
    auto D = dipoles(c);
    // int d^3k/k^0 = int k dk dOmega, at Born kinematics (soft approximation).
    // "prod" is what CalculateRealSubIF() would hand over, i.e. through
    // Dipole::Eikonal(k,p1,p2). "text" drops that overload's norm = -1 for
    // like-charge pairs, which is the only difference between it and the
    // Dipole::Eikonal(k) overload two functions below it in Dipole.C.
    double real_prod(0.), real_text(0.);
    const int NE = 400, NC = 200, NP = 64;
    const double dlk = (log(om2)-log(om1))/NE;
    for (int ie = 0; ie < NE; ++ie) {
      const double k = exp(log(om1) + (ie+0.5)*dlk);
      for (int ic = 0; ic < NC; ++ic) {
        const double ck = -1. + (ic+0.5)*2./NC, sk = sqrt(1.-ck*ck);
        for (int ip = 0; ip < NP; ++ip) {
          const double ph = (ip+0.5)*2.*M_PI/NP;
          const Vec4D kk(k, k*sk*cos(ph), k*sk*sin(ph), k*ck);
          const double meas = k*k*dlk * (2./NC) * (2.*M_PI/NP);
          for (auto &d : D) {
            const Vec4D p1 = d.GetMomenta(0), p2 = d.GetMomenta(1);
            real_prod += d.Eikonal(kk, p1, p2) * meas;
            real_text += -d.ChargeNorm()*alpha/(4.*M_PI*M_PI)
                         * (p1/(p1*kk) - p2/(p2*kk)).Abs2() * meas;
          }
        }
      }
    }
    double dB(0.);
    for (auto &d : D) {
      const Vec4D q1 = d.GetBornMomenta(1), q2 = d.GetBornMomenta(0);
      dB += d.ChargeNorm()*(
              ff.BVR_full(q1*q2, q1.E(), q2.E(), q1.Mass(), q2.Mass(),
                          om2, ff.m_photonMass, 0)
            - ff.BVR_full(q1*q2, q1.E(), q2.E(), q1.Mass(), q2.Mass(),
                          om1, ff.m_photonMass, 0));
    }
    std::cout << "\n--- does the real IF eikonal cancel the cutoff? (cos = " << c
              << ", " << om1 << " < E_gamma < " << om2 << ") ---\n"
              << "  sum_IF ChargeNorm*[Btilda(w2)-Btilda(w1)] = " << dB << "\n"
              << "  int, via Dipole::Eikonal(k,p1,p2)         = " << real_prod
              << "   ratio " << real_prod/dB << "\n"
              << "  int, same but without the norm=-1 flag    = " << real_text
              << "   ratio " << real_text/dB << "\n"
              << "  (ratio +1 => the real emission between the two cutoffs restores\n"
              << "   exactly what raising omega took out of Btilda, so omega cancels)\n";
  }

  std::cout << "\n--- exponent vs cos(theta) ---\n";
  std::cout << std::setw(8) << "cos" << std::setw(16) << "Y_was"
            << std::setw(16) << "Y_was,+sign" << std::setw(16) << "Y_now"
            << std::setw(16) << "Y_box" << "\n";
  for (int i = -9; i <= 9; i += 2) {
    const double c = 0.1*i;
    std::cout << std::setw(8) << c << std::setw(16) << Y_old(c)
              << std::setw(16) << Y_oldsign(c) << std::setw(16) << Y_kkmcsign(c)
              << std::setw(16) << Y_box(c) << "\n";
  }

  // Even part: an exactly odd exponent leaves sigma alone. Anything even in
  // here is a rate shift masquerading as an interference term.
  std::cout << "\n--- even part 0.5*(Y(+c)+Y(-c)), which shifts sigma rather than A_FB ---\n";
  std::cout << std::setw(8) << "cos" << std::setw(16) << "Y_was"
            << std::setw(16) << "Y_now" << "\n";
  for (double c : {0.1, 0.3, 0.5, 0.7, 0.9})
    std::cout << std::setw(8) << c
              << std::setw(16) << 0.5*(Y_old(c)+Y_old(-c))
              << std::setw(16) << 0.5*(Y_kkmcsign(c)+Y_kkmcsign(-c)) << "\n";

  } catch (const ATOOLS::Exception &e) {
    std::cerr << "Sherpa exception: " << e << std::endl;
    return 1;
  }
  return 0;
}
