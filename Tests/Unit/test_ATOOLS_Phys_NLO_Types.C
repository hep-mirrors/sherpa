#include <catch2/catch_all.hpp>

#include "ATOOLS/Phys/NLO_Types.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

// The KP:FACTORISATION_SCHEME setting selects the collinear scheme; DISgamma
// gates the resolved-photon pointlike-KP mask in KP_Terms::Calculate. That gate
// is only correct if the facscheme enum round-trips through YAML: SetDefault
// stringifies via ToString (operator<<) and the value is parsed back through
// operator>>. "DIS" is a substring of "DISgamma", so a regression in the
// operator>> ordering would silently collapse DISgamma to DIS (or throw).

TEST_CASE("facscheme stringifies to a single whitespace-free token",
          "[ATOOLS::facscheme]") {
  CHECK(ToString(facscheme::MSbar)    == "MSbar");
  CHECK(ToString(facscheme::DIS)      == "DIS");
  CHECK(ToString(facscheme::DISgamma) == "DISgamma");
}

TEST_CASE("facscheme round-trips ToString -> operator>>",
          "[ATOOLS::facscheme]") {
  for (auto fs : {facscheme::MSbar, facscheme::DIS, facscheme::DISgamma}) {
    facscheme::code back;
    MyStrStream str{ ToString(fs) };
    str >> back;
    CHECK(back == fs);
  }
}

TEST_CASE("facscheme operator>> keeps DISgamma distinct from DIS",
          "[ATOOLS::facscheme]") {
  // The explicit regression guard for the substring hazard.
  facscheme::code fs;
  MyStrStream("DISgamma") >> fs;
  CHECK(fs == facscheme::DISgamma);
  MyStrStream("DIS") >> fs;
  CHECK(fs == facscheme::DIS);
}

TEST_CASE("facscheme accepts the legacy integer codes",
          "[ATOOLS::facscheme]") {
  facscheme::code fs;
  MyStrStream("0") >> fs;  CHECK(fs == facscheme::MSbar);
  MyStrStream("1") >> fs;  CHECK(fs == facscheme::DIS);
  MyStrStream("2") >> fs;  CHECK(fs == facscheme::DISgamma);
}

// The read path FactorisationScheme() depends on: SetDefault(MSbar).Get<> over
// a KP:FACTORISATION_SCHEME entry. A local Settings exercises the exact
// stringify/parse round-trip without touching the main-settings singleton.
TEST_CASE("KP:FACTORISATION_SCHEME reads back through Settings",
          "[ATOOLS::facscheme]") {
  if (!msg) msg = new Message();

  SECTION("explicit DISgamma") {
    Settings s{ "KP: {FACTORISATION_SCHEME: DISgamma}" };
    CHECK(s["KP"]["FACTORISATION_SCHEME"]
            .SetDefault(facscheme::MSbar).Get<facscheme::code>()
          == facscheme::DISgamma);
  }
  SECTION("explicit MSbar") {
    Settings s{ "KP: {FACTORISATION_SCHEME: MSbar}" };
    CHECK(s["KP"]["FACTORISATION_SCHEME"]
            .SetDefault(facscheme::MSbar).Get<facscheme::code>()
          == facscheme::MSbar);
  }
  SECTION("default when absent") {
    Settings s{ "" };
    CHECK(s["KP"]["FACTORISATION_SCHEME"]
            .SetDefault(facscheme::MSbar).Get<facscheme::code>()
          == facscheme::MSbar);
  }
}
