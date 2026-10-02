#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/YFS_Process.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/My_MPI.H"

#include <cstdlib>
#include <memory>
#include <unistd.h>

namespace {
class GeneratorSettings {
  std::string m_old_path, m_path;
  std::unique_ptr<ATOOLS::My_MPI> m_mpi;

public:
  GeneratorSettings()
  {
    PHASIC_Test::SeedGlobals();
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
    // Seed only the generator setting this standalone lifecycle test needs.
    // Settings overrides are immutable; tests needing real ME backends must
    // initialize their own main settings, as other configuration tests do.
    auto setting = ATOOLS::Settings::GetMainSettings()["ME_GENERATORS"];
    m_old_path = ATOOLS::rpa->gen.Variable("SHERPA_CPP_PATH");
    char path[] = "/tmp/phasic-yfs-process-XXXXXX";
    const auto dir = mkdtemp(path);
    REQUIRE(dir != nullptr);
    m_path = dir;
    ATOOLS::rpa->gen.SetVariable("SHERPA_CPP_PATH", m_path);
    setting.OverrideVector(std::vector<std::string>{"None"});
  }

  ~GeneratorSettings()
  {
    ATOOLS::rpa->gen.SetVariable("SHERPA_CPP_PATH", m_old_path);
    rmdir((m_path + "/Process").c_str());
    rmdir(m_path.c_str());
    if (m_mpi) ATOOLS::mpi = nullptr;
  }
};

class Generator : public PHASIC::ME_Generator_Base {
public:
  Generator() : ME_Generator_Base("YFSTest") {}
  bool Initialize(MODEL::Model_Base*, BEAM::Beam_Spectra_Handler*,
                  PDF::ISR_Handler*, YFS::YFS_Handler*) override { return false; }
  PHASIC::Process_Base* InitializeProcess(const PHASIC::Process_Info&,
                                         bool) override { return nullptr; }
  int PerformTests() override { return 1; }
  bool NewLibraries() override { return false; }
};
}

TEST_CASE("YFS_Process can be destroyed before process initialization",
          "[PHASIC::YFS_Process]")
{
  GeneratorSettings settings;
  PHASIC::ME_Generators generators;
  {
    PHASIC::YFS_Process process(generators, nullptr);
    CHECK(process.p_virtproc == nullptr);
    CHECK(process.p_yfs == nullptr);
  }
  SUCCEED("The base-owned integrator is destroyed once");
}

TEST_CASE("YFS_Process exposes its generator through Process_Base",
          "[PHASIC::YFS_Process]")
{
  GeneratorSettings settings;
  PHASIC::ME_Generators generators;
  PHASIC_Test::Jet_Group jets(ATOOLS::Flavour(kf_none));
  Generator generator;
  PHASIC::YFS_Process process(generators, nullptr);
  PHASIC::Process_Base& base = process;

  base.SetGenerator(&generator);
  CHECK(base.Generator() == &generator);
  CHECK(process.Generator() == &generator);
  base.SetGenerator(nullptr);
  CHECK(base.Generator() == nullptr);
  CHECK(process.Generator() == nullptr);
}
