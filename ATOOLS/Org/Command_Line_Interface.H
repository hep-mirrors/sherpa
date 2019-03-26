#ifndef ATOOLS_Org_Command_Line_Interface_H
#define ATOOLS_Org_Command_Line_Interface_H

#include "ATOOLS/Org/Yaml_Reader.H"

#include "ATOOLS/Org/MyStrStream.H"

namespace ATOOLS {

  namespace Option_Parser {
    class Option;
    class Parser;
  }

  class Command_Line_Interface : public Yaml_Reader {

    public:

      Command_Line_Interface(int argc, char* argv[]);

    private:

      void Parse(int argc, char* argv[]);
      bool ParseOptions(std::vector<Option_Parser::Option>&);
      bool ParseNoneOptions(Option_Parser::Parser&);

      void PrintUsageAndExit();

      MyStrStream m_yamlstream;

  };

}

#endif
