#ifndef ATOOLS_ORG_STRINGS_H
#define ATOOLS_ORG_STRINGS_H

#include <string>

namespace ATOOLS {

namespace Strings {

constexpr char ProcessesSyntaxExamples[] =
R"|(A valid PROCESSES example is:

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

constexpr char NoConfigFilesWarning[] =
R"|(WARNING: No config file will be loaded. While it is possible to run Sherpa
without config files, it is likely that this is not intended. Make sure that
there is either a config file named "Sherpa.yaml" in the current working
directory, or provide a config file on the command line like this:
"Sherpa MySherpaFile.yaml".
)|";

constexpr char OverrideAlphaSSplittingScalesWarning[] =
R"|(WARNING: Shower variations are disabled (CSS_REWEIGHT: false), but the
reweighting of the AlphaS input scales in splittings is enabled
(REWEIGHT_SPLITTING_ALPHAS_SCALES: true). This would lead to an inconsistent
reweighting of the splitting scales in the shower and in the reconstructed CKKW
shower history. Therefore, the value of REWEIGHT_SPLITTING_ALPHAS_SCALES will
be overriden to be false. Set it to false explicitly or remove your custom
REWEIGHT_SPLITTING_ALPHAS_SCALES setting to silence this warning.
)|";

constexpr char OverridePDFSplittingScalesWarning[] =
R"|(WARNING: Shower variations are disabled (CSS_REWEIGHT: false), but the
reweighting of the PDF input scales in splittings is enabled
(REWEIGHT_SPLITTING_PDF_SCALES: true). This would lead to an inconsistent
reweighting of the splitting scales in the shower and in the reconstructed CKKW
shower history. Therefore, the value of REWEIGHT_SPLITTING_PDF_SCALES will
be overriden to be false. Set it to false explicitly or remove your custom
REWEIGHT_SPLITTING_PDF_SCALES setting to silence this warning.
)|";

} // namespace Strings

} // namespace ATOOLS

#endif
