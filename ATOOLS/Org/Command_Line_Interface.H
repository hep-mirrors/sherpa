#ifndef ATOOLS_Org_Command_Line_Interface_H
#define ATOOLS_Org_Command_Line_Interface_H

#include "ATOOLS/Org/Smart_Pointer.H"

#include <string>
#include <map>

namespace ATOOLS {

  namespace Option_Parser {
    class Option;
    class Parser;
  }

  class Command_Line_Interface {

    typedef std::map<std::string, std::string> String_Map;

    public:

      Command_Line_Interface(int argc, char* argv[]);

      // get parameter/tag maps
      const String_Map& GetParameters() const { return parameter_value_map; }
      const String_Map& GetTags()       const { return tag_value_map; }

      // get parameter/tag values
      std::string& GetParameterValue(const std::string& paramname)
      { return parameter_value_map[paramname]; }
      void SetParameterValue(const std::string& paramname,
                             const std::string& value,
                             const bool& allowoverwrite=true);
      std::string& GetTagValue(const std::string& tagname)
      { return tag_value_map[tagname]; }
      void SetTagValue(const std::string& tagname,
                       const std::string& value,
                       const bool& allowoverwrite=true);

      /// add arguments from a command file (where each row is an argument)
      void AddArgumentsFromCommandFile(const std::string&);

    private:

      void Parse(int argc, char* argv[], bool allowoverwrite=true);

      bool ParseOptions(SP(Option_Parser::Option)&, bool allowoverwrite=true);
      bool ParseNoneOptions(Option_Parser::Parser&, bool allowoverwrite=true);

      void PrintUsageAndExit();

      String_Map parameter_value_map, tag_value_map;

  };

}

#endif

