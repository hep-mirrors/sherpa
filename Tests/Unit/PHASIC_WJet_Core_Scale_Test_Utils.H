#ifndef TESTS_UNIT_PHASIC_WJET_CORE_SCALE_TEST_UTILS_H
#define TESTS_UNIT_PHASIC_WJET_CORE_SCALE_TEST_UTILS_H

#include <catch2/catch_all.hpp>
#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include <array>
#include <memory>

namespace PHASIC_Test {
class WJetParticles {
  std::array<ATOOLS::Particle_Info, 2> m_particles;
  std::array<ATOOLS::Particle_Info*, 2> m_previous;
  const std::array<kf_code, 2> m_codes = {kf_Wplus, kf_u};
public:
  WJetParticles()
  {
    for (size_t i = 0; i < m_codes.size(); ++i) {
      const auto it = ATOOLS::s_kftable.find(m_codes[i]);
      m_previous[i] = it == ATOOLS::s_kftable.end() ? nullptr : it->second;
      m_particles[i].m_kfc = m_codes[i];
      m_particles[i].m_mass = i == 0 ? 80. : 0.;
      m_particles[i].m_strong = i == 0 ? 0 : 3;
      ATOOLS::s_kftable[m_codes[i]] = &m_particles[i];
    }
  }
  ~WJetParticles()
  {
    for (size_t i = 0; i < m_codes.size(); ++i) {
      if (m_previous[i]) ATOOLS::s_kftable[m_codes[i]] = m_previous[i];
      else ATOOLS::s_kftable.erase(m_codes[i]);
    }
  }
};

inline void CheckWJetCoreScale(const std::string& name)
{
  SeedGlobals();
  WJetParticles metadata;
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
