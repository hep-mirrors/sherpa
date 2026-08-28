// Computes the magnitude of Define_Dipoles::FormFactorSum()'s FF-dipole term
// (D.ChargeNorm()*BVR_full(D, sqrt(m_s)/2)) for the e+e- -> pi+pi- FSR dipole
// at Born kinematics - the piece currently skipped by FormFactorSum() whenever
// HIDE_PHOTONS=1 (the default). Evaluated at fixed Born kinematics/omega, so
// this is a genuine constant for FSR-only mode.
//
// Build/run: see build_FF_FormFactor_Magnitude.sh in this directory.
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
  char override_arg[] = "YFS: {USE_MODEL_ALPHA: 0, KKMC_ANG: 1, IR_CUTOFF: 1e-9, FULL_FORM: 1}";
  char prog_name[] = "ff_formfactor_mag";
  char* fake_argv[] = {prog_name, override_arg};
  Settings::InitializeMainSettings(2, fake_argv);
  ATOOLS::ran = new Random(1234);
  ATOOLS::s_loader = new Library_Loader();
  PDF::pdfdefs = new PDF::PDF_Defaults();

  try {

  double mpi    = 0.13957039;
  double sqrts  = 0.7;
  double pmag   = sqrt(sqr(sqrts/2.)-sqr(mpi));
  rpa->gen.SetEcms(sqrts);

  Vec4D pPip(sqrts/2., 0., 0.,  pmag);
  Vec4D pPim(sqrts/2., 0., 0., -pmag);
  s_kftable[211] = new Particle_Info(211, mpi, 0., 3, 0, 0, "pi+", "\\pi^+");
  double alpha = 1./137.035999084;
  MODEL::s_model = new FakeModel(alpha);

  Flavour flPip((kf_code)211, false);
  Flavour flPim((kf_code)211, true);
  Flavour_Vector fl{flPip, flPim};
  Vec4D_Vector   mom{pPip, pPim};
  Vec4D_Vector   born{pPip, pPim};

  Dipole dip(fl, mom, born, dipoletype::final, alpha);
  dip.m_masses  = {mpi, mpi};
  dip.m_charges = {1., -1.};
  dip.SetFlavLab(0, 1);
  dip.SetResonance(true);
  dip.SetBorn(1.0);
  dip.BoostToQFM(false);

  YFS_Form_Factor form;
  form.SetCharge(1.0);

  double omega = sqrt(sqr(sqrts))/2.; // sqrt(m_s)/2
  double bvr = form.BVR_full(dip, omega);
  double chargeNorm = dip.ChargeNorm();
  double term = chargeNorm * bvr;

  std::cout << std::setprecision(10);
  std::cout << "omega = sqrt(s)/2      = " << omega << "\n";
  std::cout << "dip.ChargeNorm()       = " << chargeNorm << "\n";
  std::cout << "BVR_full(dip, omega)   = " << bvr << "\n";
  std::cout << "missing FormFactorSum term (ChargeNorm*BVR_full) = " << term << "\n";
  std::cout << "exp(term) [multiplicative constant if this term were included] = " << exp(term) << "\n";
  std::cout << "current m_formfactor for FSR-only (=exp(0), since FF skipped)  = " << exp(0.0) << "\n";
  std::cout << "relative shift if included = " << (exp(term)-1.0)*100. << " %\n";

  } catch (const ATOOLS::Exception& e) {
    std::cerr << "ATOOLS::Exception: " << e << std::endl;
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "std::exception: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
