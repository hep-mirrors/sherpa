#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Coupling_Test_Utils.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/My_MPI.H"

#include <cstdlib>
#include <memory>

namespace {
class Environment {
  PHASIC_Test::MasslessQCDParticles m_particles;
  MODEL::Running_AlphaS* p_previous_as;
  std::unique_ptr<MODEL::Running_AlphaS> m_as;
  std::unique_ptr<ATOOLS::My_MPI> m_mpi;
public:
  Environment()
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
    p_previous_as = MODEL::as;
    m_as.reset(new MODEL::Running_AlphaS(.118, 91.1876 * 91.1876, 0, 1));
    MODEL::as = m_as.get();
  }
  ~Environment()
  {
    MODEL::as = p_previous_as;
    m_as.reset();
    if (m_mpi) ATOOLS::mpi = nullptr;
  }
};
}

TEST_CASE("Fastjet scale evaluation accommodates more clustering nodes than expressions",
          "[PHASIC::Fastjet_Scale_Setter]")
{
  Environment environment;
  const ATOOLS::Flavour gluon(kf_gluon);
  PHASIC_Test::Process process(ATOOLS::Flavour_Vector(7, gluon), 2);
  PHASIC::Scale_Setter_Arguments args(nullptr, "FASTJET[M=0,R=0.1]{H_T2}");
  args.p_proc = &process;
  std::unique_ptr<PHASIC::Scale_Setter_Base> setter(
    PHASIC::Scale_Setter_Base::Scale_Getter_Function::GetObject("FASTJET", args));
  REQUIRE(setter != nullptr);
  ATOOLS::Vec4D_Vector momenta = {ATOOLS::Vec4D(100., 0., 0., 100.),
                                 ATOOLS::Vec4D(100., 0., 0., -100.)};
  double ht = 0.;
  for (int i = 0; i < 5; ++i) {
    const double pt = 10. + i;
    const double phi = 2. * M_PI * i / 5.;
    momenta.emplace_back(pt, pt * std::cos(phi), pt * std::sin(phi), 0.);
    ht += pt;
  }
  // the couplings keep pointers into the scale vector from the constructor,
  // so storing the clustering nodes must not move it
  const double* storage = setter->Scales().data();
  for (int i = 0; i < 2; ++i) {
    CAPTURE(i);
    CHECK_THAT(setter->CalculateScale(momenta), Catch::Matchers::WithinRel(ht * ht, 1.e-12));
    CHECK(setter->Scales().size() == ATOOLS::stp::size);
    CHECK(setter->Scales().data() == storage);
    for (const double scale : setter->Scales()) CHECK(std::isfinite(scale));
  }
}
