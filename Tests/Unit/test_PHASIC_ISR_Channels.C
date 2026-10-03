#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"
#include "PHASIC++/Channels/ISR_Channels.H"
#include "PHASIC++/Channels/Threshold_Channels.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PDF/Main/Intact.H"
#include "MODEL/Main/Model_Base.H"

namespace {
// Only the alpha_QED constant is consumed while constructing an off-mode YFS
// handler. Do not initialize a physical model or alter a shared model's table.
class Model : public MODEL::Model_Base {
  MODEL::Model_Base* p_previous;
public:
  Model() : Model_Base(false), p_previous(MODEL::s_model)
  {
    Constants()->emplace("alpha_QED", 1. / 137.03599976);
    MODEL::s_model = this;
  }
  ~Model() override { MODEL::s_model = p_previous; }
private:
  void InitVertices() override {}
  void ParticleInit() override {}
  bool ModelInit() override { return true; }
};

// Coulomb reads W metadata even when YFS is off; restore the shared KF slot.
class WMetadata {
  ATOOLS::Particle_Info m_w;
  ATOOLS::Particle_Info* p_previous;
public:
  WMetadata()
  {
    const auto it = ATOOLS::s_kftable.find(kf_Wplus);
    p_previous = it == ATOOLS::s_kftable.end() ? nullptr : it->second;
    m_w.m_kfc = kf_Wplus;
    ATOOLS::s_kftable[kf_Wplus] = &m_w;
  }
  ~WMetadata()
  {
    if (p_previous) ATOOLS::s_kftable[kf_Wplus] = p_previous;
    else ATOOLS::s_kftable.erase(kf_Wplus);
  }
};

class Hadron : public PDF::ISR_Base {
public:
  Hadron() : ISR_Base(nullptr)
  {
    m_type = PDF::isrtype::hadron;
    m_bunch = ATOOLS::Flavour(kf_none);
    m_on = true;
  }
  bool CalculateWeight(double, double, double, double, int) override
  { return true; }
  double Weight(ATOOLS::Flavour) override { return 1.; }
};
}

TEST_CASE("ISR_Channels constructs threshold samplers for threshold metadata",
          "[PHASIC::ISR_Channels]")
{
  PHASIC_Test::SeedGlobals();
  Model model;
  WMetadata w;
  YFS::YFS_Handler yfs;
  REQUIRE_FALSE(yfs.HasISR());
  const ATOOLS::Flavour none(kf_none);
  PDF::ISR_Handler isr({new Hadron(), new PDF::Intact(none)});
  PHASIC_Test::Process process({none, none, none, none}, 2);
  process.Integrator()->Initialize(nullptr, &isr, &yfs);
  PHASIC::Phase_Space_Handler phase_space(process.Integrator(), .01, "", "");
  PHASIC::ISR_Channels channels(&phase_space, "threshold-test");
  channels.Parameters().emplace_back(PHASIC::channel_type::threshold, 10., 2., 0.);
  REQUIRE(channels.Initialize());
  REQUIRE(channels.Number() == 1);
  CHECK(dynamic_cast<PHASIC::Threshold_Central*>(channels.Channel(0)) != nullptr);
}
