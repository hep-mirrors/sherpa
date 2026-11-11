#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "PHASIC++/Channels/ISR_Channels.H"
#include "PHASIC++/Channels/Threshold_Channels.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PDF/Main/Intact.H"

namespace {
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
  // An off-mode YFS handler still constructs Coulomb, which reads alpha_QED
  // and W metadata.
  PHASIC_Test::QEDModel model;
  PHASIC_Test::ScopedParticle w(kf_Wplus);
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
