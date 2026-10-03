#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"
#include "PHASIC++/Scales/Tag_Setter.H"

namespace {
class Scale : public PHASIC::Scale_Setter_Base {
public:
  explicit Scale(const PHASIC::Scale_Setter_Arguments& args) : Scale_Setter_Base(args) {}
private:
  double Calculate(const ATOOLS::Vec4D_Vector&, const size_t&) override
  {
    m_scale = {100., 200., 300.};
    return m_scale[0];
  }
};

class MassiveParticle {
  ATOOLS::Particle_Info m_particle;
  ATOOLS::Particle_Info* p_previous;
public:
  MassiveParticle()
  {
    auto it = ATOOLS::s_kftable.find(kf_h0);
    p_previous = it == ATOOLS::s_kftable.end() ? nullptr : it->second;
    m_particle.m_kfc = kf_h0;
    m_particle.m_mass = 4.;
    m_particle.m_massive = true;
    ATOOLS::s_kftable[kf_h0] = &m_particle;
  }
  ~MassiveParticle()
  {
    if (p_previous) ATOOLS::s_kftable[kf_h0] = p_previous;
    else ATOOLS::s_kftable.erase(kf_h0);
  }
};
}

TEST_CASE("Scale tags distinguish the massive transverse product from indexed scales",
          "[PHASIC::Tag_Setter]")
{
  PHASIC_Test::SeedGlobals();
  MassiveParticle metadata;
  PHASIC_Test::Process process({ATOOLS::Flavour(kf_none), ATOOLS::Flavour(kf_none),
                                ATOOLS::Flavour(kf_h0), ATOOLS::Flavour(kf_h0),
                                ATOOLS::Flavour(kf_none)}, 2);
  PHASIC::Scale_Setter_Arguments args(nullptr);
  args.p_proc = &process;
  Scale setter(args);
  const ATOOLS::Vec4D_Vector momenta = {ATOOLS::Vec4D(10., 0., 0., 10.),
                                       ATOOLS::Vec4D(10., 0., 0., -10.),
                                       ATOOLS::Vec4D(5., 3., 0., 0.),
                                       ATOOLS::Vec4D(4., 0., 0., 0.),
                                       ATOOLS::Vec4D(2., 0., 2., 0.)};
  setter.CalculateScale(momenta);
  PHASIC::Tag_Setter tags(&setter);
  ATOOLS::Algebra_Interpreter calculator;
  calculator.SetTagReplacer(&tags);
  tags.SetCalculator(&calculator);
  tags.SetTags(&calculator);
  SECTION("P_TM2 is the squared geometric mean over massive final legs only") {
    calculator.Interprete("P_TM2");
    CHECK_THAT(calculator.Calculate()->Get<double>(),
               Catch::Matchers::WithinRel(momenta[2].MPerp() * momenta[3].MPerp(), 1.e-12));
  }
  SECTION("indexed scales and momentum tags retain their mapping") {
    for (size_t i = 0; i < setter.Scales().size(); ++i) {
      calculator.Interprete("MU_" + std::to_string(i) + "2");
      CHECK_THAT(calculator.Calculate()->Get<double>(),
                 Catch::Matchers::WithinRel(setter.Scales()[i], 1.e-12));
    }
    calculator.Interprete("PPerp2(p[4])");
    CHECK_THAT(calculator.Calculate()->Get<double>(),
               Catch::Matchers::WithinRel(momenta[4].PPerp2(), 1.e-12));
  }
}
