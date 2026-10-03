#ifndef TESTS_UNIT_PHASIC_COUPLING_TEST_UTILS_H
#define TESTS_UNIT_PHASIC_COUPLING_TEST_UTILS_H

#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"

namespace PHASIC_Test {

// Running_AlphaS needs the light-quark thresholds as well as the gluon.
// A gluon-only table reaches ContinueAlphaS(0), whose preceding threshold
// access is outside the supported layout. Restore every temporary KF entry.
class MasslessQCDParticles {
  ATOOLS::Particle_Info m_particles[4];
  ATOOLS::Particle_Info* p_previous[4];
public:
  MasslessQCDParticles()
  {
    SeedGlobals();
    for (int i = 0; i < 4; ++i) {
      const auto code = static_cast<kf_code>(i < 3 ? i + 1 : kf_gluon);
      const auto it = ATOOLS::s_kftable.find(code);
      p_previous[i] = it == ATOOLS::s_kftable.end() ? nullptr : it->second;
      m_particles[i].m_kfc = code;
      m_particles[i].m_strong = i < 3 ? 3 : 8;
      m_particles[i].m_mass = m_particles[i].m_hmass = 0.;
      ATOOLS::s_kftable[code] = &m_particles[i];
    }
  }
  ~MasslessQCDParticles()
  {
    for (int i = 0; i < 4; ++i) {
      const auto code = static_cast<kf_code>(i < 3 ? i + 1 : kf_gluon);
      if (p_previous[i]) ATOOLS::s_kftable[code] = p_previous[i];
      else ATOOLS::s_kftable.erase(code);
    }
  }
  MasslessQCDParticles(const MasslessQCDParticles&) = delete;
  MasslessQCDParticles& operator=(const MasslessQCDParticles&) = delete;
};

} // namespace PHASIC_Test

#endif
