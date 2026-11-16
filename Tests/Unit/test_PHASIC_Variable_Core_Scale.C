#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Scale_Test_Utils.H"
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

#include <cmath>
#include <memory>
#include <string>
#include <utility>

TEST_CASE("Core-scale hadronic H_T tags sum p_T or m_T of coloured final-state legs",
          "[PHASIC::Variable_Core_Scale]")
{
  PHASIC_Test::MasslessAlphaS coupling;
  PHASIC_Test::ScopedParticle jet(kf_jet), electron(kf_e), charm(kf_c);
  electron.Info().m_icharge = -3;
  charm.SetMass(1.5);
  charm.Info().m_strong = 3;
  // e- g -> e- c cbar; the tags only read transverse quantities, so momentum
  // balance is irrelevant here
  const ATOOLS::Flavour_Vector flavours = {
      ATOOLS::Flavour(kf_e), ATOOLS::Flavour(kf_gluon), ATOOLS::Flavour(kf_e),
      ATOOLS::Flavour(kf_c), ATOOLS::Flavour(kf_c, true)};
  PHASIC_Test::Process process(flavours, 2);
  const ATOOLS::Vec4D_Vector momenta = {
      ATOOLS::Vec4D(-10., 0., 0., -10.), ATOOLS::Vec4D(-10., 0., 0., 10.),
      ATOOLS::Vec4D(10., -1., 0., std::sqrt(99.)),
      ATOOLS::Vec4D(std::sqrt(7.25), 2., 0., 1.), ATOOLS::Vec4D(3.5, -1., 0., -3.)};
  auto destroy = [](ATOOLS::Cluster_Amplitude* ampl) { ampl->Delete(); };
  std::unique_ptr<ATOOLS::Cluster_Amplitude, decltype(destroy)> amplitude(
      ATOOLS::Cluster_Amplitude::New(), destroy);
  amplitude->SetNIn(2);
  for (size_t i = 0; i < momenta.size(); ++i)
    amplitude->CreateLeg(momenta[i], flavours[i]);
  const double htm = momenta[3].MPerp() + momenta[4].MPerp();
  for (const auto& [tag, expected] : {std::pair<std::string, double>{"hH_T2", 9.},
                                      std::pair<std::string, double>{"hH_TM2", htm * htm}}) {
    CAPTURE(tag);
    const std::string scale = "VAR{" + tag + "}";
    std::unique_ptr<PHASIC::Core_Scale_Setter> setter(
        PHASIC::Core_Scale_Getter::GetObject(
            scale, PHASIC::Core_Scale_Arguments(&process, scale)));
    REQUIRE(setter != nullptr);
    const auto scales = setter->Calculate(amplitude.get());
    CHECK_THAT(scales.m_mu2, Catch::Matchers::WithinRel(expected, 1.e-12));
  }
}
