#ifndef TESTS_UNIT_PHASIC_TEST_UTILS_H
#define TESTS_UNIT_PHASIC_TEST_UTILS_H

#include <catch2/catch_all.hpp>

#include "PHASIC++/Process/Process_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/KF_Table.H"

#include <cstdlib>
#include <memory>

namespace PHASIC_Test {

// Seed the globals a bare sherpa_tests run leaves uninitialised. Seed the
// settings only here: InitializeMainSettings replaces the singleton, which
// would discard every default and override registered by earlier tests.
inline void SeedGlobals()
{
  if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
  if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
  static bool seeded = false;
  if (!seeded) {
    ATOOLS::Settings::InitializeMainSettings(std::string{"{}"});
    seeded = true;
  }
  if (ATOOLS::s_kftable.find(kf_none) == ATOOLS::s_kftable.end())
    ATOOLS::AddParticle(kf_none, 0., 0., 0., 0, 0, true, 1, "none", "none");
}

// Temporarily replace one KF-table entry, configured through Info(), and
// restore the previous entry (or its absence) on destruction. Unlike
// AddParticle, this neither throws on nor alters an entry other tests share.
class ScopedParticle {
  ATOOLS::Particle_Info m_info;
  ATOOLS::Particle_Info* p_previous;
public:
  explicit ScopedParticle(kf_code kf)
  {
    const auto it = ATOOLS::s_kftable.find(kf);
    p_previous = it == ATOOLS::s_kftable.end() ? nullptr : it->second;
    m_info.m_kfc = kf;
    ATOOLS::s_kftable[kf] = &m_info;
  }
  ~ScopedParticle()
  {
    if (p_previous) ATOOLS::s_kftable[m_info.m_kfc] = p_previous;
    else ATOOLS::s_kftable.erase(m_info.m_kfc);
  }
  ScopedParticle(const ScopedParticle&) = delete;
  ScopedParticle& operator=(const ScopedParticle&) = delete;
  ATOOLS::Particle_Info& Info() { return m_info; }
  void SetMass(double mass)
  {
    m_info.m_mass = m_info.m_hmass = mass;
    m_info.m_massive = true;
  }
};

class Jet_Group {
  ScopedParticle m_jet{kf_jet};
public:
  explicit Jet_Group(const ATOOLS::Flavour& constituent)
  {
    m_jet.Info().SetIsGroup(true);
    m_jet.Info().Add(constituent);
  }
};

class Bunches {
  ATOOLS::Flavour m_first, m_second;
public:
  Bunches(const ATOOLS::Flavour& first, const ATOOLS::Flavour& second)
    : m_first(ATOOLS::rpa->gen.Bunch(0)),
      m_second(ATOOLS::rpa->gen.Bunch(1))
  {
    ATOOLS::rpa->gen.SetBunch(first, 0);
    ATOOLS::rpa->gen.SetBunch(second, 1);
  }
  ~Bunches()
  {
    ATOOLS::rpa->gen.SetBunch(m_first, 0);
    ATOOLS::rpa->gen.SetBunch(m_second, 1);
  }
  Bunches(const Bunches&) = delete;
  Bunches& operator=(const Bunches&) = delete;
};

// File and directory helpers consult mpi->Rank(), even in a serial run.
class ScopedMPI {
  std::unique_ptr<ATOOLS::My_MPI> m_mpi;
public:
  ScopedMPI()
  {
#ifdef USING__MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
      REQUIRE(MPI_Init(nullptr, nullptr) == MPI_SUCCESS);
      std::atexit([] {
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized) MPI_Finalize();
      });
    }
#endif
    if (!ATOOLS::mpi) {
      m_mpi.reset(new ATOOLS::My_MPI());
      ATOOLS::mpi = m_mpi.get();
    }
  }
  ~ScopedMPI() { if (m_mpi) ATOOLS::mpi = nullptr; }
  ScopedMPI(const ScopedMPI&) = delete;
  ScopedMPI& operator=(const ScopedMPI&) = delete;
};

// Install a seeded ATOOLS::ran for one test and restore the previous generator
// afterwards, so the draws do not depend on which tests ran before. Random
// registers itself with the terminator handler, so deregister it before it is
// destroyed.
class ScopedRandom {
  ATOOLS::Random* p_previous;
  std::unique_ptr<ATOOLS::Random> m_ran;
public:
  explicit ScopedRandom(long seed = 1234) : p_previous(ATOOLS::ran)
  {
    if (!ATOOLS::exh) ATOOLS::exh = new ATOOLS::Terminator_Object_Handler();
    m_ran.reset(new ATOOLS::Random(seed));
    ATOOLS::ran = m_ran.get();
  }
  ~ScopedRandom()
  {
    ATOOLS::exh->RemoveTerminatorObject(m_ran.get());
    ATOOLS::ran = p_previous;
  }
  ScopedRandom(const ScopedRandom&) = delete;
  ScopedRandom& operator=(const ScopedRandom&) = delete;
};

// Install a model that supplies only alpha_QED, as a constant and as a
// coupling function, without initializing vertices or particles.
class QEDModel : public MODEL::Model_Base {
  class AlphaQED : public ATOOLS::Function_Base {
  public:
    explicit AlphaQED(double alpha) { m_name = "Alpha_QED"; m_defval = alpha; }
  };
  MODEL::Model_Base* p_previous;
public:
  explicit QEDModel(double alpha = 1. / 137.)
    : Model_Base(false), p_previous(MODEL::s_model)
  {
    Constants()->emplace("alpha_QED", alpha);
    p_functions->emplace("alpha_QED", new AlphaQED(alpha));
    MODEL::s_model = this;
  }
  ~QEDModel() override { MODEL::s_model = p_previous; }
private:
  void InitVertices() override {}
  void ParticleInit() override {}
  bool ModelInit() override { return true; }
};

// Supply external-leg metadata without constructing an ME generator. The
// event-generation members are inert stubs.
class Process : public PHASIC::Process_Base {
public:
  Process(const ATOOLS::Flavour_Vector& flavours, size_t nin)
  {
    m_flavs = flavours;
    m_nin = nin;
    m_nout = flavours.size() - nin;
  }
  void SetScale(const PHASIC::Scale_Setter_Arguments&) override {}
  void SetKFactor(const PHASIC::KFactor_Setter_Arguments&) override {}
  void InitializeTheReweighting(ATOOLS::Variations_Mode) override {}
  size_t Size() const override { return 1; }
  PHASIC::Process_Base* operator[](const size_t&) override { return this; }
  ATOOLS::Weight_Info* OneEvent(int, ATOOLS::Variations_Mode, int) override
  { return nullptr; }
  ATOOLS::Weights_Map Differential(const ATOOLS::Vec4D_Vector&,
                                   ATOOLS::Variations_Mode) override
  { return ATOOLS::Weights_Map{0.0}; }
  bool CalculateTotalXSec(const std::string&, bool) override { return false; }
  void SetLookUp(bool lookup) override { m_lookup = lookup; }
};

} // namespace PHASIC_Test

#endif
