#ifndef ATOOLS_Org_Default_Reader_H
#define ATOOLS_Org_Default_Reader_H

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Data_Writer.H"
#include "ATOOLS/Org/Message.H"

namespace ATOOLS {

  /**
   * Use this Data_Reader wrapper for consistent defaults retrieval.
   *
   * If WRITE_DEFAULTS_FILE is set to 1, a requested_defaults.mmd is written
   * after the run with a list of all requested defaults and their default and
   * custom values (if present).
   *
   * Examples for reading defaults:
   *
   *   Default_Reader reader;
   *   Type value = reader.Get("PARAMETER_NAME", default_value);
   *
   * If you provide a (human readable) parameter name and a method name, an info
   * message is printed if a custom value is used instead of the default one:
   *
   *   Type value = reader.Get("PARAMETER_NAME", default_value,
   *                           "parameter name", METHOD);
   *
   * would print "CallingMethod(): Set parameter name = value." using msg_Info().
   *
   * Use Read instead of Get if you need to know if a custom value has been
   * read (as opposed to falling back to the default value):
   *
   *   Type value;
   *   bool didread = reader.Read(&value, "PARAMETER_NAME", default_value);
   *
   * For parameters that accept disabling values like "None", "Off" and "0",
   * a convenience function can (and should) be used:
   *
   *   if (reader.IsDisabled("PARAMETER_NAME")) { ... }
   */
  class Default_Reader: public Data_Reader {

    // facilities for preparing a Markdown output of all requested parameters
    typedef std::map<std::string, std::pair<std::string, std::string> > RequestedDefaults;
    static RequestedDefaults s_requesteddefaults;
    static bool s_shouldwritedefaults;
    std::string EncodeForMarkdown(const std::string &) const;
    void Configure();

  public:
    // constructors
    Default_Reader() : Data_Reader(" ", ";", "#", "=") {
      AddComment("!");
      AddWordSeparator("\t");
      Configure();
    };

    bool IsDisabled(const std::string &parameter, const bool &def=false);

    template <class Read_Type> bool
    Read(Read_Type &val, const std::string &param, const Read_Type def)
    {
      const bool success = Data_Reader::ReadFromFile(val, param);
      std::string val_string;
      if (success) {
        val_string = ToString<Read_Type>(val);
      } else {
        val = def;
        val_string = "";
      }
      if (s_shouldwritedefaults) {
        s_requesteddefaults[param] = std::make_pair(ToString<Read_Type>(def),
                                                    val_string);
      }
      return success;
    }

    template <class Read_Type> bool
    ReadVector(std::vector<Read_Type> &val, const std::string &param)
    {
      VectorFromFile<Read_Type>(val, param);
      if (s_shouldwritedefaults) {
        s_requesteddefaults[param] = std::make_pair(std::string(""),
                                                    VectorToString<Read_Type>(val));
      }
      return (val.size() != 0);
    }

    bool ReadStringVectorNormalisingNoneLikeValues(std::vector<std::string> &val,
                                                   const std::string &param);

    template <class Read_Type> bool
    ReadMatrix(std::vector<std::vector<Read_Type> > &val,
               const std::string &param)
    {
      MatrixFromFile<Read_Type>(val, param);
      if (s_shouldwritedefaults) {
        s_requesteddefaults[param] = std::make_pair(std::string("<empty>"),
                                                    MatrixToString<Read_Type>(val));
      }
      return (val.size() != 0);
    }

    template <class Read_Type> Read_Type
    Get(const std::string &param,
        const Read_Type &def,
        const std::string &outputname="",
        const std::string &methodname="")
    {
      Read_Type val;
      if (Read(val, param, def) && outputname != "") {
        if (methodname != "") {
          msg_Info() << methodname << "(): ";
        }
        msg_Info() << "Set " << outputname << " = " << val << "." << std::endl;
      }
      return val;
    }

    /// call this if you are done using this class
    void Finalize() const;

  }; // end of class Default_Reader

} // end of namespace ATOOLS

#endif // ifndef ATOOLS_Org_Default_Reader_H
