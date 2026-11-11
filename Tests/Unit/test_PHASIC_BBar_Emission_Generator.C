#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"
#include "PHASIC++/Channels/BBar_Emission_Generator.H"
#include "ATOOLS/Org/My_MPI.H"

#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <new>

namespace {
void SeedSettings()
{
  PHASIC_Test::SeedGlobals();
  ATOOLS::Settings::GetMainSettings()["DIPOLES"]["AMIN"].SetDefault(
    ATOOLS::Max(ATOOLS::rpa->gen.Accu(), 1.e-8));
}

class CacheEnvironment {
  std::unique_ptr<ATOOLS::My_MPI> m_mpi;
  std::string m_directory;
public:
  CacheEnvironment()
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
    char directory[] = "/tmp/phasic-eeg-cache-XXXXXX";
    const auto path = mkdtemp(directory);
    REQUIRE(path != nullptr);
    m_directory = path;
  }
  ~CacheEnvironment()
  {
    std::filesystem::remove_all(m_directory);
    if (m_mpi) ATOOLS::mpi = nullptr;
  }
  std::string Prefix() const { return m_directory + "/cache"; }
  void Write(const std::string& contents) const
  {
    // Legacy File_IO_Base concatenates path and file verbatim; this is a
    // filename prefix, not a directory (matching the production writer).
    std::ofstream file(Prefix() + "_CS_EEG_PV");
    REQUIRE(file.is_open());
    file << contents;
    REQUIRE(file.good());
  }
};
}

TEST_CASE("BBar_Emission_Generator starts without an active emission or weight",
          "[PHASIC::BBar_Emission_Generator]")
{
  SeedSettings();
  using Generator = PHASIC::BBar_Emission_Generator;
  // Pattern ordinary placement-new storage so allocator zeroing cannot hide
  // missing constructor initialization. Do not inspect any private members.
  alignas(Generator) std::array<unsigned char, sizeof(Generator)> storage;
  storage.fill(0x3f);
  auto destroy = [](Generator* generator) { generator->~Generator(); };
  std::unique_ptr<Generator, decltype(destroy)> generator(
    new (storage.data()) Generator(), destroy);
  REQUIRE_THAT(generator->Weight(), Catch::Matchers::WithinAbs(0., 1.e-14));
  const auto active = generator->Active(nullptr);
  CHECK(active.p_dip == nullptr);
  CHECK(active.p_p == nullptr);
  CHECK(active.m_procs.empty());
  CHECK_THAT(active.m_weight, Catch::Matchers::WithinAbs(0., 1.e-14));
  CHECK_FALSE(generator->GenerateWeight(nullptr));
}

TEST_CASE("BBar_Emission_Generaistinguishes unavailable and malformed optimization data",
          "[PHASIC::BBar_Emission_Generator]")
{
  SeedSettings();
  CacheEnvironment cache;
  PHASIC::BBar_Emission_Generator generator;
  SECTION("a missing metadata file returns false") {
    CHECK_FALSE(generator.ReadIn(cache.Prefix()));
  }
  SECTION("an empty metadata file returns false") {
    cache.Write("");
    CHECK_FALSE(generator.ReadIn(cache.Prefix()));
  }
  SECTION("an optimization row with extra fields still throws") {
    cache.Write("5 6\n");
    CHECK_THROWS_AS(generator.ReadIn(cache.Prefix()), ATOOLS::fatal_error);
  }
  SECTION("a scalar optimization row is accepted without dipoles") {
    cache.Write("5\n");
    CHECK(generator.ReadIn(cache.Prefix()));
  }
}
