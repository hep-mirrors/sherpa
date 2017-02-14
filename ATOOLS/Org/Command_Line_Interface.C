#include "ATOOLS/Org/Command_Line_Interface.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Option_Parser.H"
#include "ATOOLS/Org/Data_Reader.H"

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

void Command_Line_Interface::Parse(int argc, char* argv[], bool allowoverwrite)
{
  // prepare static memory for parser operations
  const bool allow_mixing_options_and_nonoptions(true);
  Option_Parser::Stats  stats(allow_mixing_options_and_nonoptions, usage,
                              argc, argv);
  SP(Option_Parser::Option) options(
      new Option_Parser::Option[stats.options_max]);
  SP(Option_Parser::Option) buffer(
      new Option_Parser::Option[stats.buffer_max]);

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
  success = (ParseNoneOptions(parser, allowoverwrite) && success);
  success = (ParseOptions(options,    allowoverwrite) && success);

  if (!success) {
    PrintUsageAndExit();
  }
}

bool Command_Line_Interface::ParseOptions(SP(Option_Parser::Option)& options,
                                          bool allowoverwrite)
{
  bool success(true);

  // check unknown
  if ((--options)[UNKNOWN]) {
    for (Option_Parser::Option* opt = (--options)[UNKNOWN];
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
    const char* arg((--options)[it->second].last()->arg);
    if (arg) {
      SetParameterValue(it->first, std::string(arg), allowoverwrite);
    } else if ((--options)[it->second].last()->type() == DISABLE) {
      SetParameterValue(it->first, "0", allowoverwrite);
    } else if ((--options)[it->second].last()->type() == ENABLE) {
      SetParameterValue(it->first, "1", allowoverwrite);
    }
  }

  return success;
}

bool Command_Line_Interface::ParseNoneOptions(Option_Parser::Parser& parser,
                                              bool allowoverwrite)
{
  bool success(true);

  for (int i = 0; i < parser.nonOptionsCount(); ++i) {

    std::string nonOption(parser.nonOption(i));

    // process tag definition
    size_t tagdelimiterpos = nonOption.find(":=");
    if (tagdelimiterpos != std::string::npos) {
      std::string tagname(nonOption.substr(0, tagdelimiterpos));
      std::string tagvalue(nonOption.substr(tagdelimiterpos + 2));
      SetTagValue(tagname, tagvalue, allowoverwrite);
      continue;
    }

    // process parameter definition
    size_t paramdelimiterpos = nonOption.find("=");
    if (paramdelimiterpos != std::string::npos) {
      std::string paramname(nonOption.substr(0, paramdelimiterpos));
      std::string paramvalue(nonOption.substr(paramdelimiterpos + 1));
      SetParameterValue(paramname, paramvalue, allowoverwrite);
      continue;
    }

    // malformed parameter
    msg_Error() << "ERROR: Can not parse argument: '" << nonOption << "'"
                << std::endl;
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

void Command_Line_Interface::SetParameterValue(const std::string& paramname,
                                               const std::string& value,
                                               const bool& allowoverwrite)
{
  if (parameter_value_map[paramname] == "" || allowoverwrite) {
    parameter_value_map[paramname] = value;
  }
}

void Command_Line_Interface::SetTagValue(const std::string& tagname,
                                         const std::string& value,
                                         const bool& allowoverwrite)
{
  if (parameter_value_map[tagname] == "" || allowoverwrite) {
    parameter_value_map[tagname] = value;
  }
}

void Command_Line_Interface::AddArgumentsFromCommandFile(
    const std::string& file)
{
  // read command file
  Data_Reader reader;
  reader.SetInputFile(file);
  String_Matrix args;
  reader.MatrixFromFile(args);

  int argc(args.size());
  char** argv = new char*[argc];

  // copy rows into argv
  size_t j(0);
  typedef String_Matrix::const_iterator It;
  for (It it(args.begin()); it != args.end(); ++it) {
    std::string row = std::accumulate(it->begin(), it->end(), std::string(""));
    argv[j] = new char[row.size() + 1];
    strcpy(argv[j], row.c_str());
    ++j;
  }

  // parse, but do not overwrite existing values, because the command line
  // should have precedence over the command file
  Parse(argc, argv, false);

  for (int i(0); i < argc; ++i) delete [] argv[i];
  delete [] argv;
}

