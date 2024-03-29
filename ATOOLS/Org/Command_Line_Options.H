#ifndef ATOOLS_Org_Command_Line_Options_H
#define ATOOLS_Org_Command_Line_Options_H

#include "ATOOLS/Org/Option_Parser.H"

// NOTE: to add an option, declare it in the `optionIndex' enum, define it in
// the `usage' array below, and define its mapping to a Sherpa parameter in the
// `create_parameter_index_map' function below

/// declaration of command line options
enum optionIndex {
  // special options
  UNKNOWN,
  HELP,
  VERSION,
  CMDLINE_DEBUG,

  // options that will be mapped to SHERPA parameters
  RUNDATA,
  PATH,
  EVENTS,
  EVENT_TYPE,
  RESULT_DIRECTORY,
  SHERPA_CPP_PATH,
  RANDOM_SEED,
  ME_GENERATORS,
  MI_HANDLER,
  EVENT_GENERATION_MODE,
  SHOWER_GENERATOR,
  FRAGMENTATION,
  ANALYSIS,
  ANALYSIS_OUTPUT,
  OUTPUT,
  EVT_OUTPUT,
  MSG_LIMIT,
  LOG_FILE,

  // options used to enable/disable features
  BATCH_MODE,
  INIT_ONLY,
  GENERATE_RESULT_DIRECTORY,
  PRINT_VERSION_INFO
};

/// declaration of command line option types
enum optionType  {REGULAR, DISABLE, ENABLE};

/// actual definition of command line arguments and usage output
const ATOOLS::Option_Parser::Descriptor usage[] = {
  {
    UNKNOWN,
    REGULAR,
    "",
    "",
    ATOOLS::Option_Parser::Arg::None,
    "USAGE: Sherpa [options] ['<parameter>: <value>'] ['TAGS: {<tag>: <value>}'] [Sherpa.yaml]\n\nOptions:"
  },
  {
    HELP,
    REGULAR,
    "h",
    "help",
    ATOOLS::Option_Parser::Arg::None,
    "  --help, -h  \tPrint usage and exit."
  },
  {
    VERSION,
    REGULAR,
    "v",
    "version",
    ATOOLS::Option_Parser::Arg::None,
    "  --version, -v  \tPrint version and exit."
  },
  {
    CMDLINE_DEBUG,
    REGULAR,
    "c",
    "cmdline-debug",
    ATOOLS::Option_Parser::Arg::None,
    "  --cmdline-debug, -c  \tPrint command line after parsing to YAML and exit."
  },
  {
    RUNDATA,
    REGULAR,
    "f",
    "run-data",
    ATOOLS::Option_Parser::Arg::Required,
    "  --run-data, -f <file>  \tRead run data from <file>. [DEPRECATED, use positional arguments instead]"
  },
  {
    PATH,
    REGULAR,
    "p",
    "path",
    ATOOLS::Option_Parser::Arg::Required,
    "  --path, -p <path>  \tRead run data in <path>."
  },
  {
    EVENTS,
    REGULAR,
    "e",
    "events",
    ATOOLS::Option_Parser::Arg::Required,
    "  --events, -e <N_events>  \tSet number of events <N_events>."
  },
  {
    EVENT_TYPE,
    REGULAR,
    "e",
    "event-type",
    ATOOLS::Option_Parser::Arg::Required,
    "  --event-type, -t <event_type>  \tSet event type <event_type>."
  },
  {
    RESULT_DIRECTORY,
    REGULAR,
    "r",
    "result-directory",
    ATOOLS::Option_Parser::Arg::Required,
    "  --result-directory, -r <path>  \tSet result directory <path>."
  },
  {
    SHERPA_CPP_PATH,
    REGULAR,
    "L",
    "sherpa-lib-path",
    ATOOLS::Option_Parser::Arg::Required,
    "  --sherpa-lib-path, -L <path>  \tSet <path> for dynamically created C++ libs."
  },
  {
    RANDOM_SEED,
    REGULAR,
    "R",
    "random-seed",
    ATOOLS::Option_Parser::Arg::Required,
    "  --random-seed, -R <seed(s)>  \tSet random <seed(s)>."
  },
  {
    ME_GENERATORS,
    REGULAR,
    "m",
    "me-generators",
    ATOOLS::Option_Parser::Arg::Required,
    "  --me-generators, -m <generator>  \tSet matrix element <generator>."
  },
  {
    MI_HANDLER,
    REGULAR,
    "M",
    "mi-handler",
    ATOOLS::Option_Parser::Arg::Required,
    "  --mi-handler, -M <handler>  \tSet multiple interaction <handler>."
  },
  {
    EVENT_GENERATION_MODE,
    REGULAR,
    "w",
    "event-generation-mode",
    ATOOLS::Option_Parser::Arg::Required,
    "  --event-generation-mode, -w <mode>  \tSet event generation <mode>."
  },
  {
    SHOWER_GENERATOR,
    REGULAR,
    "s",
    "shower-generator",
    ATOOLS::Option_Parser::Arg::Required,
    "  --shower-generator, -s <generator>  \tSet parton-shower <generator>."
  },
  {
    FRAGMENTATION,
    REGULAR,
    "F",
    "fragmentation",
    ATOOLS::Option_Parser::Arg::Required,
    "  --fragmentation, -F <module>  \tSet fragmentation <module>."
  },
  {
    ANALYSIS,
    REGULAR,
    "a",
    "analyis",
    ATOOLS::Option_Parser::Arg::Required,
    "  --analyis, -a <handler>  \tSet the analysis <handler>."
  },
  {
    ANALYSIS_OUTPUT,
    REGULAR,
    "A",
    "analyis-output",
    ATOOLS::Option_Parser::Arg::Required,
    "  --analyis-output, -A <path>  \tSet the analysis output <path>."
  },
  {
    OUTPUT,
    REGULAR,
    "O",
    "output",
    ATOOLS::Option_Parser::Arg::Required,
    "  --output, -O <level>  \tSet the output <level>."
  },
  {
    EVT_OUTPUT,
    REGULAR,
    "o",
    "event-output",
    ATOOLS::Option_Parser::Arg::Required,
    "  --event-output, -o <level>  \tSet the event output <level>."
  },
  {
    MSG_LIMIT,
    REGULAR,
    "W",
    "msg-limit",
    ATOOLS::Option_Parser::Arg::Required,
    "  --msg-limit, -W <n>  \tSet number of warnings emitted before silencing them."
  },
  {
    LOG_FILE,
    REGULAR,
    "l",
    "log-file",
    ATOOLS::Option_Parser::Arg::Required,
    "  --log-file, -l <file>  \tSet the log <file>."
  },
  {
    BATCH_MODE,
    DISABLE,
    "b",
    "disable-batch-mode",
    ATOOLS::Option_Parser::Arg::None,
    "  --disable-batch-mode, -b  \tDisable batch mode."
  },
  {
    INIT_ONLY,
    ENABLE,
    "I",
    "enable-init-only",
    ATOOLS::Option_Parser::Arg::None,
    "  --enable-init-only, -I  \tExit after completing initialisation."
  },
  {
    GENERATE_RESULT_DIRECTORY,
    DISABLE,
    "g",
    "disable-result-directory-generation",
    ATOOLS::Option_Parser::Arg::None,
    "  --disable-result-directory-generation, -g  \tDisable result directory generation."
  },
  {
    PRINT_VERSION_INFO,
    ENABLE,
    "V",
    "print-version-info",
    ATOOLS::Option_Parser::Arg::None,
    "  --print-version-info, -V  \tPrint version info during runtime."
  },
  {
    UNKNOWN,
    REGULAR,
    "",
    "",
    ATOOLS::Option_Parser::Arg::None,
    "\nExamples:\n"
    "  Sherpa -e 1k\n"
    "  Sherpa -w Unweighted MySherpaConfig.yaml\n"
  },
  {0, 0, 0, 0, 0, 0}
};

// mapping from parameter names to option indizes
typedef std::map<std::string, optionIndex> String_Option_Map;
String_Option_Map create_parameter_index_map() {
  static String_Option_Map map;
  if (map.size() == 0) {
    map["RUNDATA"]                   = RUNDATA;
    map["PATH"]                      = PATH;
    map["EVENTS"]                    = EVENTS;
    map["EVENT_TYPE"]                = EVENT_TYPE;
    map["RESULT_DIRECTORY"]          = RESULT_DIRECTORY;
    map["SHERPA_CPP_PATH"]           = SHERPA_CPP_PATH;
    map["RANDOM_SEED"]               = RANDOM_SEED;
    map["ME_GENERATORS"]             = ME_GENERATORS;
    map["MI_HANDLER"]                = MI_HANDLER;
    map["EVENT_GENERATION_MODE"]     = EVENT_GENERATION_MODE;
    map["SHOWER_GENERATOR"]          = SHOWER_GENERATOR;
    map["FRAGMENTATION"]             = FRAGMENTATION;
    map["ANALYSIS"]                  = ANALYSIS;
    map["ANALYSIS_OUTPUT"]           = ANALYSIS_OUTPUT;
    map["OUTPUT"]                    = OUTPUT;
    map["EVT_OUTPUT"]                = EVT_OUTPUT;
    map["MSG_LIMIT"]                 = MSG_LIMIT;
    map["LOG_FILE"]                  = LOG_FILE;
    map["BATCH_MODE"]                = BATCH_MODE;
    map["INIT_ONLY"]                 = INIT_ONLY;
    map["PRINT_VERSION_INFO"]        = PRINT_VERSION_INFO;
    map["GENERATE_RESULT_DIRECTORY"] = GENERATE_RESULT_DIRECTORY;
  }
  return map;
}
const String_Option_Map parameter_index_map =
create_parameter_index_map();

#endif
