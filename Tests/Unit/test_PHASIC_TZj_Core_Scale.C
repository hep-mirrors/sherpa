#include "Tests/Unit/PHASIC_WJet_Core_Scale_Test_Utils.H"

TEST_CASE("TZj W-plus-jet core uses the W transverse mass in either leg order",
          "[PHASIC::TZj_Core_Scale]")
{
  PHASIC_Test::CheckWJetCoreScale("TZj");
}
