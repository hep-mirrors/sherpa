#ifndef TESTS_UNIT_PHASIC_SCALE_TEST_UTILS_H
#define TESTS_UNIT_PHASIC_SCALE_TEST_UTILS_H

#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "MODEL/Main/Running_AlphaS.H"

#include <deque>
#include <memory>
#include <string>

namespace PHASIC_Test {

// Install a Running_AlphaS over massless d, u, s quarks and the gluon as
// MODEL::as. The light-quark thresholds are needed as well as the gluon: a
// gluon-only table reaches ContinueAlphaS(0), whose preceding threshold
// access is outside the supported layout.
class MasslessAlphaS {
  std::deque<ScopedParticle> m_particles;
  MODEL::Running_AlphaS* p_previous;
  std::unique_ptr<MODEL::Running_AlphaS> m_as;
public:
  MasslessAlphaS() : p_previous(MODEL::as)
  {
    SeedGlobals();
    for (const kf_code kf : {kf_d, kf_u, kf_s, kf_gluon}) {
      m_particles.emplace_back(kf);
      m_particles.back().Info().m_strong = kf == kf_gluon ? 8 : 3;
    }
    m_as.reset(new MODEL::Running_AlphaS(.118, 91.1876 * 91.1876, 0, 1));
    MODEL::as = m_as.get();
  }
  ~MasslessAlphaS() { MODEL::as = p_previous; }
  MasslessAlphaS(const MasslessAlphaS&) = delete;
  MasslessAlphaS& operator=(const MasslessAlphaS&) = delete;
};

// Global inputs of the clustering scale setters (METS, E_Gamma, VBF, MINLO):
// alphaS, their settings defaults, a 1 TeV collider and unit scale factors.
class ScaleEnvironment {
  MasslessAlphaS m_as;
  double m_ecms;
  std::string m_factors[3];
  inline static const std::string factor_names[3] = {
    "RENORMALIZATION_SCALE_FACTOR", "FACTORIZATION_SCALE_FACTOR", "RESUMMATION_SCALE_FACTOR"};
public:
  ScaleEnvironment()
  {
    auto& settings = ATOOLS::Settings::GetMainSettings();
    settings["ASSOCIATED_CONTRIBUTIONS_VARIATIONS"].SetDefault(
      std::vector<std::vector<std::string>>{});
    settings["DIPOLES"]["NF_GSPLIT"].SetDefault(0);
    settings["SHOWER"]["KFACTOR_SCHEME"].SetDefault(1);
    m_ecms = ATOOLS::rpa->gen.Ecms();
    ATOOLS::rpa->gen.SetEcms(1000.);
    for (int i = 0; i < 3; ++i) {
      try { m_factors[i] = ATOOLS::rpa->gen.Variable(factor_names[i]); }
      catch (const ATOOLS::fatal_error&) { m_factors[i] = "1"; }
      ATOOLS::rpa->gen.SetVariable(factor_names[i], "1");
    }
  }
  ~ScaleEnvironment()
  {
    for (int i = 0; i < 3; ++i)
      ATOOLS::rpa->gen.SetVariable(factor_names[i], m_factors[i]);
    ATOOLS::rpa->gen.SetEcms(m_ecms);
  }
  ScaleEnvironment(const ScaleEnvironment&) = delete;
  ScaleEnvironment& operator=(const ScaleEnvironment&) = delete;
};

} // namespace PHASIC_Test

#endif
