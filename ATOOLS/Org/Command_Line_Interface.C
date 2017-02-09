#include "ATOOLS/Org/Command_Line_Interface.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Option_Parser.H"

// include the actual definition of command line options
// NOTE: go to this header if you want to modify or add options
#include "ATOOLS/Org/Command_Line_Options.H"

using namespace ATOOLS;

Command_Line_Interface::Command_Line_Interface(int argc, char* argv[])
{
  // skip program name argv[0] if present
  if (argc > 0) {
    --argc;
    ++argv;
  }

  // prepare static memory for parser operations
  const bool allow_mixing_options_and_nonoptions(true);
  Option_Parser::Stats  stats(allow_mixing_options_and_nonoptions, usage,
                              argc, argv);
  SP(Option_Parser::Option) options(new Option_Parser::Option[stats.options_max]);
  SP(Option_Parser::Option) buffer(new Option_Parser::Option[stats.buffer_max]);

  // parse arguments
  Option_Parser::Parser parser(allow_mixing_options_and_nonoptions, usage,
                               argc, argv, --options, --buffer);

  if (parser.error()) {
    THROW(fatal_error, "Command line syntax error");
  }

  if ((--options)[HELP]) {
    PrintUsageAndExit();
  }

  if ((--options)[VERSION]) {
    msg_Out() << "Sherpa version "
              << SHERPA_VERSION << "." << SHERPA_SUBVERSION
              << " (" << SHERPA_NAME << ")" << std::endl;
    exit(0);
  }

  // parse parameters; order matters here, we want to give options precedence
  // over non-options
  bool success(true);
  success = (ParseNoneOptions(parser) && success);
  success = (ParseOptions(options)    && success);

  if (!success) {
    PrintUsageAndExit();
  }

  /* exit(0); */
}

bool Command_Line_Interface::ParseOptions(SP(Option_Parser::Option)& options)
{
  bool success(true);

  // check unknown
  if ((--options)[UNKNOWN]) {
    for (Option_Parser::Option* opt = (--options)[UNKNOWN]; opt; opt = opt->next()) {
      msg_Error() << "ERROR: Unknown option: '" << opt->name << "'" << std::endl;
    }
    success = false;
  }

  // fill parameter->value map
  typedef String_Option_Map::const_iterator It;
  for (It it(parameter_index_map.begin());
      it != parameter_index_map.end();
      ++it) {
    const char* arg((--options)[it->second].last()->arg);
    if (arg) {
      parameter_value_map[it->first] = std::string(arg);
    } else if ((--options)[it->second].last()->type() == DISABLE) {
      parameter_value_map[it->first] = "0";
    } else if ((--options)[it->second].last()->type() == ENABLE) {
      parameter_value_map[it->first] = "1";
    }
  }

  return success;
}

bool Command_Line_Interface::ParseNoneOptions(Option_Parser::Parser& parser)
{
  bool success(true);

  for (int i = 0; i < parser.nonOptionsCount(); ++i) {

    std::string nonOption(parser.nonOption(i));

    // process tag definition
    size_t tagdelimiterpos = nonOption.find(":=");
    if (tagdelimiterpos != std::string::npos) {
      std::string tagname(nonOption.substr(0, tagdelimiterpos));
      std::string tagvalue(nonOption.substr(tagdelimiterpos + 2));
      tag_value_map[tagname] = tagvalue;
      continue;
    }

    // process parameter definition
    size_t paramdelimiterpos = nonOption.find("=");
    if (paramdelimiterpos != std::string::npos) {
      std::string paramname(nonOption.substr(0, paramdelimiterpos));
      std::string paramvalue(nonOption.substr(paramdelimiterpos + 1));
      parameter_value_map[paramname] = paramvalue;
      continue;
    }

    // malformed parameter
    msg_Error() << "ERROR: Can not parse argument: '" << nonOption << "'" << std::endl;
    success = false;
  }

  return success;
}

void Command_Line_Interface::PrintUsageAndExit()
{
  msg_Out() << std::endl;
  int columns = getenv("COLUMNS") ? atoi(getenv("COLUMNS")) : 80;
  Option_Parser::printUsage(msg_Out(), usage, columns);
  exit(0);
}

