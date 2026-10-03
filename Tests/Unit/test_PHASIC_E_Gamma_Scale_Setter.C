#include "Tests/Unit/PHASIC_Merging_Scale_Test_Utils.H"

TEST_CASE("E_Gamma rejects an unterminated inline option header",
          "[PHASIC::E_Gamma_Scale_Setter]")
{
  PHASIC_Test::CheckMalformedScaleHeader("E_Gamma");
}

TEST_CASE("E_Gamma owns the unordered core requested by its inline option",
          "[PHASIC::E_Gamma_Scale_Setter]")
{
  PHASIC_Test::CheckInlineUnorderedScale("E_Gamma");
}
