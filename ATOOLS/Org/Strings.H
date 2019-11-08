#ifndef ATOOLS_ORG_STRINGS_H
#define ATOOLS_ORG_STRINGS_H

#include <string>

namespace ATOOLS {

namespace Strings {

constexpr char ProcessesSyntaxExamples[] = R"|(A valid PROCESSES example is:

PROCESSES:
- 93 93 -> 11 -11:
    Order: {QCD: 0, EW: 2}
    Print_Graphs: Graphs
- 93 93 -> 13 -13:
    Order: {QCD: 0, EW: 2}

This defines two processes, each with some process-specific settings.
Note that the process-specific settings must be indented with respect
to the first character of the process specification. See the "Processes"
section in the manual for a more detailed description.
)|";

} // namespace Strings

} // namespace ATOOLS

#endif
