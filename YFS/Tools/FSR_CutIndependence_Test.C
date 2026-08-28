// Standalone MC-average test: does FSR's total per-event weight, averaged
// over many random photon configurations at FIXED Born kinematics, actually
// come out independent of FSR_CUT (the technical crude/hard-photon
// generation cutoff)? Physically it must - FSR_CUT is a purely technical
// parameter separating "hidden"/resummed soft photons from explicit ones,
// and the whole point of the crude/exact mass-weight + hideW volume-factor
// construction is that it cancels out of the total cross section.
//
// Usage: FSR_CutIndependence_Test <FSR_CUT value> <n events> [seed]
// Run twice with different FSR_CUT values (same seed) and compare the
// printed average +/- MC error.
//
// Build/run: see build_FSR_CutIndependence_Test.sh in this directory.
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
#include "YFS/Main/FSR.H"
#include "MODEL/Main/Model_Base.H"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <string>

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
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <FSR_CUT value> <n events> [seed] [IR_CUTOFF]\n";
    return 1;
  }
  std::string fsrcut_str = argv[1];
  long nevents = atol(argv[2]);
  int seed = argc > 3 ? atoi(argv[3]) : 1234;
  std::string ircutoff_str = argc > 4 ? argv[4] : "7e-11";

  ATOOLS::mpi = new My_MPI();
  ATOOLS::exh = new Terminator_Object_Handler();
  ATOOLS::msg = new Message();
  ATOOLS::rpa = new Run_Parameter();
  std::string override_arg_s =
    "YFS: {USE_MODEL_ALPHA: 0, KKMC_ANG: 1, IR_CUTOFF: " + ircutoff_str + ", FSR_CUT: " + fsrcut_str + "}";
  std::vector<char> override_arg(override_arg_s.begin(), override_arg_s.end());
  override_arg.push_back('\0');
  char prog_name[] = "sherpa_fsr_cut_test";
  char* fake_argv[] = {prog_name, override_arg.data()};
  Settings::InitializeMainSettings(2, fake_argv);
  ATOOLS::ran = new Random(seed);
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

  FSR fsr;
  fsr.SetV(0.0);
  std::cout << std::setprecision(10);

  double sumw = 0.0, sumw2 = 0.0;
  double sum_massW = 0.0, sum_hideW = 0.0, sum_wt2 = 0.0, sum_n = 0.0;
  long naccepted_finite = 0;
  long nrejected = 0;
  // mxx histogram, matching the rivet-plots binning: 300-710 MeV, 80 bins.
  const int NBINS = 80;
  const double MXX_LO = 300.0, MXX_HI = 710.0;
  double mxx_hist[NBINS] = {0};
  double mxx_hist_cut[NBINS] = {0};
  const double PMIN = 0.45*0.7/2.; // strong.cc's cmd analysis cut: pm/pp < 0.45*sqrts/2 -> veto
  for (long iev = 0; iev < nevents; ++iev) {
    // A fresh Dipole every event, matching production (YFS_Handler builds a
    // new Dipole per event) - Dipole's constructor is what clears m_ghost for
    // dipoletype::final; reusing one Dipole object across events left stale
    // ghosts from a previous event sitting at GetGhost(0)/(1) once
    // AddToGhosts() had been called more than twice total.
    Dipole dip(fl, mom, born, dipoletype::final, alpha);
    dip.m_masses  = {mpi, mpi};
    dip.m_charges = {1., -1.};
    dip.SetFlavLab(0, 1);
    dip.SetResonance(true);
    dip.SetBorn(1.0);
    dip.BoostToQFM(false);

    fsr.Reset();
    if (!fsr.Initialize(dip)) { std::cerr << "FSR::Initialize failed\n"; return 1; }
    fsr.m_mass = {mpi, mpi};

    bool ok = fsr.MakeFSR();
    double w = 0.0;
    if (ok) {
      dip.SetMomentum(0, fsr.m_dipole[0]);
      dip.SetMomentum(1, fsr.m_dipole[1]);
      dip.AddToGhosts(fsr.m_r1);
      dip.AddToGhosts(fsr.m_r2);
      if (!fsr.F()) {
        w = 0.0;
      } else {
        dip.SetNPhoton((int)fsr.m_photons.size());
        Vec4D_Vector photons = fsr.m_photons;
        dip.AddPhotonsToDipole(photons);
        dip.Boost();
        if (!fsr.YFS_FORM()) {
          w = 0.0;
        } else {
          fsr.HidePhotons();
          fsr.Weight();
          w = fsr.GetWeight();
          if (!IsBad(w)) {
            double mxx_mev = sqrt(fsr.m_sQ) * 1000.0;
            int bin = (int)((mxx_mev - MXX_LO) / (MXX_HI - MXX_LO) * NBINS);
            if (bin >= 0 && bin < NBINS) {
              mxx_hist[bin] += w;
              Vec4D pPlus  = dip.GetNewMomenta(0);
              Vec4D pMinus = dip.GetNewMomenta(1);
              double pp = Vec3D(pPlus).Abs();
              double pm = Vec3D(pMinus).Abs();
              if (pp >= PMIN && pm >= PMIN) mxx_hist_cut[bin] += w;
            }
          }
          if (!IsBad(w) && !IsBad(fsr.m_massW) && !IsBad(fsr.m_hideW) && !IsBad(fsr.m_wt2)) {
            sum_massW += fsr.m_massW;
            sum_hideW += fsr.m_hideW;
            sum_wt2   += fsr.m_wt2;
            sum_n     += fsr.m_photons.size();
            naccepted_finite++;
          }
        }
      }
    } else {
      nrejected++;
    }
    if (IsBad(w)) w = 0.0;
    sumw  += w;
    sumw2 += w*w;
  }

  double avg = sumw / nevents;
  double var = sumw2/nevents - avg*avg;
  double err = sqrt(std::max(var,0.0) / nevents);
  std::cout << "IR_CUTOFF = " << ircutoff_str << "  FSR_CUT (raw, pre-sqrt(s) scaling) = " << fsrcut_str << "\n";
  std::cout << "actual fsr.m_fsrcut (dimensionless) = " << fsr.m_fsrcut << "  fsr.m_Emin (GeV) = " << fsr.m_Emin << "\n";
  std::cout << "n events = " << nevents << "  n rejected (below-threshold) = " << nrejected << "\n";
  std::cout << "average total FSR weight = " << avg << " +/- " << err << "\n";
  std::cout << "--- per-factor averages over accepted+finite events (n=" << naccepted_finite << ") ---\n";
  std::cout << "<massW> = " << sum_massW/naccepted_finite << "\n";
  std::cout << "<hideW> = " << sum_hideW/naccepted_finite << "\n";
  std::cout << "<wt2>   = " << sum_wt2/naccepted_finite << "\n";
  std::cout << "<n_photons> = " << sum_n/naccepted_finite << "\n";
  std::cout << "--- mxx histogram (MeV), 300-710, 80 bins: uncut, then with cmd's pp/pm>" << PMIN << " cut ---\n";
  for (int i = 0; i < NBINS; ++i) {
    double lo = MXX_LO + i*(MXX_HI-MXX_LO)/NBINS;
    std::cout << lo << "  " << mxx_hist[i] << "  " << mxx_hist_cut[i] << "\n";
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
