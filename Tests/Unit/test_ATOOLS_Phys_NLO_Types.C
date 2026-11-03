#include <catch2/catch_all.hpp>

#include "Test_Setup.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace ATOOLS;

// KP:FACTORISATION_SCHEME is read as a facscheme enum, so operator<< and
// operator>> must round-trip exactly (SetDefault stringifies via ToString)
TEST_CASE("facscheme round-trips through ToString and operator>>",
          "[ATOOLS::facscheme]") {
  for (auto fs : {facscheme::MSbar, facscheme::DIS, facscheme::DISgamma}) {
    facscheme::code back;
    MyStrStream str{ ToString(fs) };
    str >> back;
    CHECK(back == fs);
  }
}

TEST_CASE("facscheme accepts the integer codes", "[ATOOLS::facscheme]") {
  facscheme::code fs;
  MyStrStream("0") >> fs;  CHECK(fs == facscheme::MSbar);
  MyStrStream("1") >> fs;  CHECK(fs == facscheme::DIS);
  MyStrStream("2") >> fs;  CHECK(fs == facscheme::DISgamma);
}

TEST_CASE("KP:FACTORISATION_SCHEME reads back through Settings",
          "[ATOOLS::facscheme]") {
  Test_Setup::InitGlobals();
  Settings s{ "KP: {FACTORISATION_SCHEME: DISgamma}" };
  CHECK(s["KP"]["FACTORISATION_SCHEME"]
          .SetDefault(facscheme::MSbar).Get<facscheme::code>()
        == facscheme::DISgamma);
}
