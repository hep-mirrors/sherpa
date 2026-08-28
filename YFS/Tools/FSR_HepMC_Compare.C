// Reads real Sherpa-generated (post-recoil dipole, final photons) tuples from
// a parsed HepMC sample (see /tmp .../parse_hepmc.py) and, for each event,
// computes the FSR mass-weight (F()) and form-factor (YFS_FORM()) using the
// REAL FSR class - exactly as production does - to check internal
// self-consistency of the formula across many real phase-space points
// (rather than the hand-picked points used in FSR_KKMC_CrossCheck.C).
//
// Input file format (one event per block):
//   EVENT weight=... nphot=N
//   BORN_PIP px py pz E      <- pi+ momentum AFTER FSR recoil (Sherpa's m_dipole[0])
//   BORN_PIM px py pz E      <- pi- momentum AFTER FSR recoil (Sherpa's m_dipole[1])
//   PHOTON px py pz E        <- one line per photon (N lines)
//
// Build/run: see build_FSR_HepMC_Compare.sh in this directory.
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "PDF/Main/PDF_Base.H"
#include "YFS/Main/Dipole.H"
#include "YFS/Main/FSR.H"
#include "MODEL/Main/Model_Base.H"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
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

struct Event {
  double weight;
  int nphot;
  Vec4D pip, pim;
  Vec4D_Vector photons;
};

static std::vector<Event> ReadEvents(const std::string &path) {
  std::vector<Event> events;
  std::ifstream fin(path);
  std::string line;
  Event cur;
  bool have = false;
  while (std::getline(fin, line)) {
    std::istringstream iss(line);
    std::string tag;
    iss >> tag;
    if (tag == "EVENT") {
      if (have) events.push_back(cur);
      cur = Event();
      have = true;
      std::string wtok, ntok;
      iss >> wtok >> ntok;
      cur.weight = atof(wtok.substr(wtok.find('=')+1).c_str());
      cur.nphot = atoi(ntok.substr(ntok.find('=')+1).c_str());
    } else if (tag == "BORN_PIP") {
      double px,py,pz,E; iss >> px >> py >> pz >> E;
      cur.pip = Vec4D(E,px,py,pz);
    } else if (tag == "BORN_PIM") {
      double px,py,pz,E; iss >> px >> py >> pz >> E;
      cur.pim = Vec4D(E,px,py,pz);
    } else if (tag == "PHOTON") {
      double px,py,pz,E; iss >> px >> py >> pz >> E;
      cur.photons.push_back(Vec4D(E,px,py,pz));
    }
  }
  if (have) events.push_back(cur);
  return events;
}

int main(int argc, char* argv[]) {
#ifdef USING__MPI
  MPI_Init(&argc, &argv);
#endif
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <events file> [IR_CUTOFF]\n";
    return 1;
  }
  std::string eventfile = argv[1];
  std::string ircutoff_str = argc > 2 ? argv[2] : "1e-9"; // matches Sherpa.Pion.FSR.yaml

  ATOOLS::mpi = new My_MPI();
  ATOOLS::exh = new Terminator_Object_Handler();
  ATOOLS::msg = new Message();
  ATOOLS::rpa = new Run_Parameter();
  std::string override_arg_s =
    "YFS: {USE_MODEL_ALPHA: 0, KKMC_ANG: 1, IR_CUTOFF: " + ircutoff_str + ", FSR_CRU: 1}";
  std::vector<char> override_arg(override_arg_s.begin(), override_arg_s.end());
  override_arg.push_back('\0');
  char prog_name[] = "sherpa_fsr_hepmc_compare";
  char* fake_argv[] = {prog_name, override_arg.data()};
  Settings::InitializeMainSettings(2, fake_argv);
  ATOOLS::ran = new Random(1234);
  ATOOLS::s_loader = new Library_Loader();
  PDF::pdfdefs = new PDF::PDF_Defaults();

  try {

  double mpi_mass = 0.13957039;
  double sqrts = 0.7;
  double pmag = sqrt(sqr(sqrts/2.)-sqr(mpi_mass));
  rpa->gen.SetEcms(sqrts);

  Vec4D bornPip(sqrts/2., 0., 0.,  pmag);
  Vec4D bornPim(sqrts/2., 0., 0., -pmag);
  s_kftable[211] = new Particle_Info(211, mpi_mass, 0., 3, 0, 0, "pi+", "\\pi^+");
  double alpha = 1./137.035999084;
  MODEL::s_model = new FakeModel(alpha);

  Flavour flPip((kf_code)211, false);
  Flavour flPim((kf_code)211, true);
  Flavour_Vector fl{flPip, flPim};
  Vec4D_Vector   mom{bornPip, bornPim};
  Vec4D_Vector   born{bornPip, bornPim};

  auto events = ReadEvents(eventfile);
  std::cerr << "Read " << events.size() << " events\n";

  FSR fsr;
  fsr.SetV(0.0);
  std::cout << std::setprecision(10);
  std::cout << "# idx nphot massW hideW wt2 total m_Emin m_fsrcut\n";

  for (size_t idx = 0; idx < events.size(); ++idx) {
    Event &ev = events[idx];

    Dipole dip(fl, mom, born, dipoletype::final, alpha);
    dip.m_masses  = {mpi_mass, mpi_mass};
    dip.m_charges = {1., -1.};
    dip.SetFlavLab(0, 1);
    dip.SetResonance(true);
    dip.SetBorn(1.0);
    dip.BoostToQFM(false);

    fsr.Reset();
    if (!fsr.Initialize(dip)) { std::cerr << "Initialize failed at " << idx << "\n"; continue; }
    fsr.m_mass = {mpi_mass, mpi_mass};

    // ---- inject the real (dipole, photons) tuple, bypassing NPhotons()/
    // GeneratePhotonMomentum()/RescalePhotons() entirely, since these ARE
    // already the post-rescale physical momenta. ----
    fsr.m_dist1.clear(); fsr.m_dist2.clear(); fsr.m_del1.clear(); fsr.m_del2.clear();
    fsr.m_cos.clear(); fsr.m_sin.clear(); fsr.m_fbarvec.clear(); fsr.m_MassWls.clear();
    fsr.m_photons = ev.photons;
    fsr.m_photonspreboost = ev.photons;
    fsr.m_photonSum = Vec4D();
    for (auto &k : ev.photons) fsr.m_photonSum += k;
    fsr.m_dipole = {ev.pip, ev.pim};
    fsr.m_Q = ev.pip + ev.pim;
    fsr.m_sQ = fsr.m_Q.Abs2();
    fsr.m_sprim = fsr.m_sQ;
    fsr.m_n = (int)ev.photons.size();
    // Reconstruct wt2 = yy*(1+preE) from the given (post-rescale) photon sum:
    // yy = sQ/dip_sp by definition (m_sQ = m_dip_sp*m_yy always), and the
    // pre-rescale photon energy fraction is preE = photonSum.E()/ener with
    // ener = sqrt(sQ)/2 (RescalePhotons()'s own "ener" scale).
    {
      double yy = fsr.m_sQ / fsr.m_dip_sp;
      double ener = sqrt(fsr.m_sQ)/2.;
      double preE = ev.photons.empty() ? 0.0 : fsr.m_photonSum.E()/ener;
      fsr.m_wt2 = ev.photons.empty() ? 1.0 : yy*(1.+preE);
    }

    // reconstruct m_cos[i]: the boost chain from the XFM (angle-generation)
    // frame to this (lab==CM, symmetric beams) frame is a pure rotation for
    // this fixed-sqrt(s), no-ISR setup (all Poincare boosts involved have
    // their reference 4-vector already at rest, i.e. (sqrts,0,0,0)), so
    // relative angles between photons and dipole[0] are preserved and can be
    // read off directly from the given (already-final) momenta.
    Vec3D d0hat = Vec3D(ev.pip) * (1./Vec3D(ev.pip).Abs());
    for (auto &k : ev.photons) {
      Vec3D khat = Vec3D(k) * (1./Vec3D(k).Abs());
      double c = d0hat * khat;
      double s = sqrt(std::max(0., 1.-c*c));
      fsr.m_cos.push_back(c);
      fsr.m_sin.push_back(s);
      double del1_s = 1. - fsr.m_beta1*c;
      double del2_s = 1. + fsr.m_beta2*c;
      fsr.m_fbarvec.push_back(1./(del1_s*del2_s) * (1.+fsr.m_beta1*fsr.m_beta2)/2.);
      fsr.m_MassWls.push_back(1.0);
    }

    bool fok = fsr.F();
    if (!fok) { std::cerr << "F() failed at " << idx << "\n"; continue; }

    // ghost momenta: same direction as dipole[0]/[1] (MakePair always builds
    // along its own z axis, but since dipole[0]/[1] are themselves related to
    // that same z axis by a pure rotation here, and BVR_cru/BVR_full only
    // depend on r1.r2 (invariant), r1[0]/r2[0] (energies, frame-dependent but
    // rotation-invariant since it's a 0-component... wait energies ARE
    // rotation-invariant) - only the ENERGIES and the INVARIANT dot product
    // r1*r2 matter for BVR_cru/BVR_full, not the explicit direction, so no
    // rotation of r1/r2 into the dipole's actual direction is needed at all.
    double masc1 = fsr.m_mass[0] * sqrt(fsr.m_sQ / fsr.m_dip_sp);
    double masc2 = fsr.m_mass[1] * sqrt(fsr.m_sQ / fsr.m_dip_sp);
    double eta1_ghost, eta2_ghost;
    fsr.MakePair(sqrt(fsr.m_sQ), fsr.m_r1, fsr.m_r2, masc1, masc2, eta1_ghost, eta2_ghost);
    fsr.CalculateBetaBar();
    dip.SetMomentum(0, fsr.m_dipole[0]);
    dip.SetMomentum(1, fsr.m_dipole[1]);
    dip.AddToGhosts(fsr.m_r1);
    dip.AddToGhosts(fsr.m_r2);
    dip.SetNPhoton((int)ev.photons.size());
    Vec4D_Vector photons_for_dip = ev.photons;
    dip.AddPhotonsToDipole(photons_for_dip);
    // Boost() is what actually populates m_newmomenta (GetNewMomenta(), read
    // by YFS_FORM()) from m_momenta - without it m_newmomenta stays at its
    // constructor-time (fixed Born) value. For this fixed-sqrt(s), no-ISR,
    // symmetric-beam setup the boost/rotate built at BoostToXFM() time are
    // all near-identity (born pair already at rest, already along z), so this
    // just correctly copies momenta through / sets m_photonSum.
    dip.Boost();

    bool yok = fsr.YFS_FORM();
    if (!yok) { std::cerr << "YFS_FORM() failed at " << idx << "\n"; continue; }
    fsr.HidePhotons();
    fsr.Weight();

    double total = fsr.GetWeight();

    std::cout << idx << " " << ev.nphot << " " << fsr.m_massW << " " << fsr.m_hideW
               << " " << fsr.m_wt2 << " " << total << " " << fsr.m_Emin << " " << fsr.m_fsrcut << "\n";
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
