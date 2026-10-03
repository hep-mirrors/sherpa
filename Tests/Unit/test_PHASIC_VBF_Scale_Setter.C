#include "Tests/Unit/PHASIC_Merging_Scale_Test_Utils.H"

TEST_CASE("VBF rejects an unterminated inline option header",
          "[PHASIC::VBF_Scale_Setter]")
{
  PHASIC_Test::CheckMalformedScaleHeader("VBF");
}
