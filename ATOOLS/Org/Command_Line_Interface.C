#include "ATOOLS/Org/Command_Line_Interface.H"

#include "ATOOLS/Org/Option_Parser.H"
#include "ATOOLS/Org/Exception.H"

// include the actual definition of command line options
// NOTE: go to this header if you want to modify or add options
#include "ATOOLS/Org/Command_Line_Options.H"

#include <numeric>
#include <cstring>

using namespace ATOOLS;

Command_Line_Interface::Command_Line_Interface(int argc, char* argv[])
{
  // skip program name argv[0] if present
  if (argc > 0) {
    --argc;
    ++argv;
  }

  Parse(argc, argv);
}

void Command_Line_Interface::Parse(int argc, char* argv[])
{
  // prepare static memory for parser operations
  const auto allow_mixing_options_and_nonoptions = true;
  Option_Parser::Stats  stats(allow_mixing_options_and_nonoptions,
                              usage,
                              argc, argv);
  std::vector<Option_Parser::Option> options(stats.options_max);
  std::vector<Option_Parser::Option> buffer(stats.buffer_max);

  // parse arguments
  Option_Parser::Parser parser(allow_mixing_options_and_nonoptions,
                               usage,
                               argc, argv,
                               &options.front(), &buffer.front());

  if (parser.error()) {
    msg_Error() << "Command line syntax error.\n";
    PrintUsageAndExit();
  }

  if (options[HELP]) {
    PrintUsageAndExit();
  }

  if (options[VERSION]) {
    msg_Out() << "Sherpa version "
              << SHERPA_VERSION << "." << SHERPA_SUBVERSION
              << " (" << SHERPA_NAME << ")" << std::endl;
    exit(0);
  }

  // parse parameters; order matters here, we want to give options precedence
  // over non-options
  bool success(true);
  success = (ParseOptions(options)    && success);
  success = (ParseNoneOptions(parser) && success);

  if (!success) {
    PrintUsageAndExit();
  }

  Yaml_Reader::Parse(m_yamlstream);
}

bool Command_Line_Interface::ParseOptions(
    std::vector<Option_Parser::Option> &options)
{
  bool success(true);

  // check unknown
  if (options[UNKNOWN]) {
    for (Option_Parser::Option* opt = options[UNKNOWN];
         opt;
         opt = opt->next()) {
      msg_Error() << "ERROR: Unknown option: '" << opt->name << "'"
                  << std::endl;
    }
    success = false;
  }

  // fill parameter->value map
  typedef String_Option_Map::const_iterator It;
  for (It it(parameter_index_map.begin());
      it != parameter_index_map.end();
      ++it) {
    const char* clivalue(options[it->second].last()->arg);
    std::string finalvalue;
    if (clivalue) {
      finalvalue = std::string(clivalue);
    } else if (options[it->second].last()->type() == DISABLE) {
      finalvalue = "0";
    } else if (options[it->second].last()->type() == ENABLE) {
      finalvalue = "1";
    }
    if (finalvalue != "")
      m_yamlstream << it->first << ": " << finalvalue << "\n";
  }

  return success;
}

bool Command_Line_Interface::ParseNoneOptions(Option_Parser::Parser& parser)
{
  String_Map legacysyntaxtags;
  auto didfinddyamltags = false;
  for (int i {parser.nonOptionsCount()-1}; i >= 0; --i) {
    auto nonOption = StringTrim(parser.nonOption(i));

    // find legacy-syntax tag specifications
    const auto pos = nonOption.find(":=");
    if (pos != std::string::npos) {
      const auto tagname = nonOption.substr(0, pos);
      const auto tagvalue = nonOption.substr(pos + 2);
      legacysyntaxtags[tagname] = tagvalue;
      continue;
    } else if (nonOption.substr(0, 4) == "TAGS") {
      didfinddyamltags = true;
    }

    // convert legacy-delimiter '=' into yaml-delimiter ':'
    const auto equalpos = nonOption.find('=');
    const auto colonpos = nonOption.find(':');
    if (equalpos != std::string::npos && colonpos == std::string::npos) {
      nonOption[equalpos] = ':';
    }

    // convert ":" into ": ", as this allows to specify tags on the command
    // line without spaces (and hence potentially without quotation marks)
    nonOption = StringReplace(nonOption, ":", ": ");

    m_yamlstream << nonOption << '\n';
  }

  if (!legacysyntaxtags.empty()) {
    // first guard against mixed yaml/legacy tag specifications
    if (didfinddyamltags) {
      msg_Error()
          << "You can not specify tags on the command line"
          <<  " using both the yaml-style \"TAGS: {X: Y}\" and the legacy-style"
          <<  " \"X:=Y\" syntaxes. Please only use one kind of syntax.\n";
      exit(1);
    }
    m_yamlstream << "TAGS: {";
    bool first = true;
    for (const auto tag : legacysyntaxtags) {
      if (first)
        first = false;
      else
        m_yamlstream << ", ";
      m_yamlstream << tag.first << ": " << tag.second;
    }
    m_yamlstream << "}\n";
  }

  return true;
}

void Command_Line_Interface::PrintUsageAndExit()
{
  msg_Out() << std::endl;
  int columns = getenv("COLUMNS") ? atoi(getenv("COLUMNS")) : 80;
  Option_Parser::printUsage(msg_Out(), usage, columns);
  exit(0);
}
