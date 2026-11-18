#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Scale_Test_Utils.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include <cmath>
#include <memory>
#include <string>

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

TEST_CASE("Fastjet scale tags pair the jets with the flavours of the clustered partons",
          "[PHASIC::Fastjet_Scale_Setter]")
{
  PHASIC_Test::MasslessAlphaS coupling;
  PHASIC_Test::ScopedMPI mpi;
  PHASIC_Test::ScopedParticle bottom(kf_b), photon(kf_photon);
  bottom.SetMass(4.75);
  bottom.Info().m_strong = 3;
  const ATOOLS::Flavour g(kf_gluon), b(kf_b), a(kf_photon);
  // g g -> b bbar gamma with massive b quarks in Sherpa's order, the massive
  // quarks before the photon; the setter lists the unclustered photon before
  // the jets, each b quark is its own jet
  PHASIC_Test::Process process({g, g, b, b.Bar(), a}, 2);
  const auto massive = [](double m, double pt, double y, double phi) {
    const double mt = std::sqrt(m * m + pt * pt);
    return ATOOLS::Vec4D(mt * std::cosh(y), pt * std::cos(phi), pt * std::sin(phi),
                         mt * std::sinh(y));
  };
  const ATOOLS::Vec4D_Vector p = {ATOOLS::Vec4D(100., 0., 0., 100.),
                                 ATOOLS::Vec4D(100., 0., 0., -100.),
                                 massive(4.75, 20., 0.5, 0.), massive(4.75, 30., -0.3, 2.5),
                                 PHASIC_Test::Massless(10., 1., -2.)};
  std::string scale;
  double expected = 0.;
  SECTION("P_TM2 multiplies the transverse masses of the b jets") {
    scale = "FASTJET[M=0,R=0.1]{P_TM2}";
    expected = p[2].MPerp() * p[3].MPerp();
  }
  SECTION("hH_T2 sums the b jets, not the photon") {
    scale = "FASTJET[M=0,R=0.1]{hH_T2}";
    expected = std::pow(p[2].PPerp() + p[3].PPerp(), 2);
  }

  PHASIC::Scale_Setter_Arguments args(nullptr, scale);
  args.p_proc = &process;
  std::unique_ptr<PHASIC::Scale_Setter_Base> setter(
    PHASIC::Scale_Setter_Base::Scale_Getter_Function::GetObject("FASTJET", args));
  REQUIRE(setter != nullptr);
  CHECK_THAT(setter->CalculateScale(p), Catch::Matchers::WithinRel(expected, 1.e-10));
}
