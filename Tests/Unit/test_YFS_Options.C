#include <catch2/catch_test_macros.hpp>

#include "ATOOLS/Org/Message.H"
#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/YFS_Options.H"

#include <sstream>
#include <string>
#include <vector>

namespace {

  void BootMessages()
  {
    if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
  }

  template <typename Code>
  Code ParseTag(const std::string &tag)
  {
    BootMessages();
    std::istringstream str(tag);
    Code code;
    str >> code;
    return code;
  }

  template <typename Code>
  std::string PrintCode(Code code)
  {
    std::ostringstream str;
    str << code;
    return str.str();
  }

  template <typename Code>
  void CheckRoundTrip(const std::vector<Code> &values)
  {
    for (Code v : values) {
      INFO("value " << static_cast<int>(v) << " printed as " << PrintCode(v));
      CHECK(ParseTag<Code>(PrintCode(v)) == v);
      CHECK(ParseTag<Code>(std::to_string(static_cast<int>(v))) == v);
    }
  }

}

TEST_CASE("YFS options: every value survives print and parse", "[YFS][Options]") {
  using namespace YFS;
  CheckRoundTrip<tristate::code>({tristate::automatic, tristate::off, tristate::on});
  CheckRoundTrip<realmap::code>({realmap::rest_frame, realmap::beam_axis,
                                 realmap::scaled, realmap::invariant});
  CheckRoundTrip<realfsrmap::code>({realfsrmap::rest_frame, realfsrmap::rescale_all,
                                    realfsrmap::dipole, realfsrmap::pre_emission});
  CheckRoundTrip<realsubeik::code>(
    {realsubeik::coherent_point, realsubeik::crude, realsubeik::coherent_born,
     realsubeik::crude_fsr, realsubeik::post_emission_fsr,
     realsubeik::coherent_born_fsr, realsubeik::assignment_born_fsr,
     realsubeik::multichannel, realsubeik::multichannel_prefsr});
  // a table with a gap: 0, 1, 4
  CheckRoundTrip<sub8fsrsub::code>({sub8fsrsub::generation, sub8fsrsub::point,
                                    sub8fsrsub::blend});
}

TEST_CASE("YFS options: names are case-insensitive", "[YFS][Options]") {
  using namespace YFS;
  CHECK(ParseTag<realsubeik::code>("MultiChannel_PreFSR") == realsubeik::multichannel_prefsr);
  CHECK(ParseTag<realmap::code>("SCALED") == realmap::scaled);
  CHECK(ParseTag<tristate::code>("Auto") == tristate::automatic);
}

TEST_CASE("YFS options: aliases parse but are not printed", "[YFS][Options]") {
  using namespace YFS;
  CHECK(ParseTag<tristate::code>("true") == tristate::on);
  CHECK(ParseTag<tristate::code>("False") == tristate::off);
  CHECK(PrintCode(tristate::on) == "on");
  // REAL_FSR_MAP: 4 is an old integer that behaves as rescale_all
  CHECK(ParseTag<realfsrmap::code>("4") == realfsrmap::rescale_all);
  CHECK(PrintCode(realfsrmap::rescale_all) == "rescale_all");
}

TEST_CASE("YFS options: unknown values are fatal", "[YFS][Options]") {
  using namespace YFS;
  CHECK_THROWS(ParseTag<realmap::code>("bogus"));
  CHECK_THROWS(ParseTag<realmap::code>("7"));
  CHECK_THROWS(ParseTag<realsubeik::code>(""));
  // inside the table's range but not one of its values
  CHECK_THROWS(ParseTag<sub8fsrsub::code>("2"));
  // a name must match whole
  CHECK_THROWS(ParseTag<realmap::code>("scale"));
}

TEST_CASE("YFS: MODE parses its own names, case-insensitively", "[YFS][Options]") {
  using namespace YFS;
  CHECK(ParseTag<yfsmode::code>("Off") == yfsmode::off);
  CHECK(ParseTag<yfsmode::code>("ISR") == yfsmode::isr);
  CHECK(ParseTag<yfsmode::code>("isrfsr") == yfsmode::isrfsr);
  CHECK(ParseTag<yfsmode::code>("ISR+FSR") == yfsmode::isrfsr);
  CHECK(ParseTag<yfsmode::code>("fsr") == yfsmode::fsr);
  for (yfsmode::code m : {yfsmode::off, yfsmode::isr, yfsmode::isrfsr, yfsmode::fsr})
    CHECK(ParseTag<yfsmode::code>(PrintCode(m)) == m);
}

TEST_CASE("YFS: MODE leaves YFS off for the values Photons++ reads", "[YFS][Options]") {
  using namespace YFS;
  // Photons++ shares the key: its values must neither throw nor switch YFS on
  for (const char *tag : {"None", "Soft", "Full", "0", "1", "2", "3", ""}) {
    INFO("YFS: MODE " << tag);
    yfsmode::code m(yfsmode::isr);
    CHECK_NOTHROW(m = ParseTag<yfsmode::code>(tag));
    CHECK(m == yfsmode::off);
  }
}
