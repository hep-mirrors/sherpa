#ifndef TESTS_UNIT_PHASIC_WJET_CORE_SCALE_TEST_UTILS_H
#define TESTS_UNIT_PHASIC_WJET_CORE_SCALE_TEST_UTILS_H

#include <catch2/catch_all.hpp>
#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include <memory>

namespace PHASIC_Test {
inline void CheckWJetCoreScale(const std::string& name)
{
  SeedGlobals();
  ScopedParticle w(kf_Wplus), quark(kf_u);
  w.SetMass(80.);
  quark.Info().m_strong = 3;
  std::unique_ptr<PHASIC::Core_Scale_Setter> setter(
    PHASIC::Core_Scale_Getter::GetObject(name, PHASIC::Core_Scale_Arguments(nullptr)));
  REQUIRE(setter != nullptr);
  for (bool wfirst : {true, false}) {
    CAPTURE(name, wfirst);
    auto destroy = [](ATOOLS::Cluster_Amplitude* ampl) { ampl->Delete(); };
    std::unique_ptr<ATOOLS::Cluster_Amplitude, decltype(destroy)> amplitude(
      ATOOLS::Cluster_Amplitude::New(), destroy);
    amplitude->SetNIn(2);
    const ATOOLS::Vec4D w(std::sqrt(80. * 80. + 30. * 30.), 30., 0., 0.);
    const ATOOLS::Vec4D jet(30., -30., 0., 0.);
    const double energy = (w[0] + jet[0]) / 2.;
    amplitude->CreateLeg(ATOOLS::Vec4D(-energy, 0., 0., -energy), ATOOLS::Flavour(kf_u));
    amplitude->CreateLeg(ATOOLS::Vec4D(-energy, 0., 0., energy), ATOOLS::Flavour(kf_u).Bar());
    amplitude->CreateLeg(wfirst ? w : jet, ATOOLS::Flavour(wfirst ? kf_Wplus : kf_u));
    amplitude->CreateLeg(wfirst ? jet : w, ATOOLS::Flavour(wfirst ? kf_u : kf_Wplus));
    const auto scales = setter->Calculate(amplitude.get());
    CHECK_THAT(scales.m_op, Catch::Matchers::WithinRel(w.MPerp2(), 1.e-12));
    CHECK_THAT(scales.m_kt2, Catch::Matchers::WithinRel(w.MPerp2(), 1.e-12));
    CHECK_THAT(scales.m_mu2, Catch::Matchers::WithinRel(jet.PPerp2(), 1.e-12));
  }
}
}
#endif
