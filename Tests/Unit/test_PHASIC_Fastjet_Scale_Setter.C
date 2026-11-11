#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Scale_Test_Utils.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include <memory>

TEST_CASE("Fastjet scale evaluation accommodates more clustering nodes than expressions",
          "[PHASIC::Fastjet_Scale_Setter]")
{
  PHASIC_Test::MasslessAlphaS coupling;
  PHASIC_Test::ScopedMPI mpi;
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
