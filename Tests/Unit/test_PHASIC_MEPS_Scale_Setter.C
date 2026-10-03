#include "Tests/Unit/PHASIC_Merging_Scale_Test_Utils.H"

TEST_CASE("METS owns the unordered core requested by its inline option",
          "[PHASIC::MEPS_Scale_Setter]")
{
  PHASIC_Test::CheckInlineUnorderedScale("METS");
}
