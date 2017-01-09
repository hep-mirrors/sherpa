#include "ATOOLS/Org/Default_Reader.H"

#include "ATOOLS/Org/MyStrStream.H"

using std::string;
using std::vector;
using namespace ATOOLS;

Default_Reader::RequestedDefaults Default_Reader::s_requesteddefaults;

bool Default_Reader::s_shouldwritedefaults = false;

string Default_Reader::EncodeForMarkdown(const string &data) const
{
  string buffer;
  buffer.reserve(data.size());
  for(size_t pos = 0; pos != data.size(); ++pos) {
    switch(data[pos]) {
      case '\n':
        buffer.append("<br />");
        break;
      case '|':
      case '_':
      case '*':
      case '\\':
      case '`':
      case '{':
      case '}':
      case '[':
      case ']':
      case '(':
      case ')':
      case '#':
      case '+':
      case '-':
      case '.':
      case '!':
      case '<':
      case '>':
        buffer.append("\\");
      default:
        buffer.append(&data[pos], 1);
        break;
    }
  }
  return buffer;
}

void Default_Reader::Configure()
{
  static bool was_executed = false;
  if (was_executed) {
    return;
  }
  if (!Data_Reader::ReadFromFile(s_shouldwritedefaults, "WRITE_DEFAULTS_FILE")) {
    s_shouldwritedefaults = false;
  }
  was_executed = true;
}

bool Default_Reader::IsDisabled(const string &parameter, const bool &def)
{
  const string value = GetStringNormalisingNoneLikeValues(parameter, "");
  if (def && value == "") return true;
  if (value == "None")    return true;
  return false;
}

bool Default_Reader::ReadStringVectorNormalisingNoneLikeValues(vector<string> &val,
                                                               const string &param)
{
  StringVectorFromFileNormalisingNoneLikeValues(val, param);
  if (s_shouldwritedefaults) {
    s_requesteddefaults[param] = std::make_pair(string("<empty>"),
                                                VectorToString<string>(val));
  }
  return (val.size() != 0);
}

void Default_Reader::Finalize() const
{
  if (!s_shouldwritedefaults) {
    return;
  }
  Data_Writer writer("","","", "");
  writer.SetOutputFile("requested_defaults.mmd");

  // write table caption and header
  writer.WriteToFile(string("[Defaults used in the last run, "));
  writer.WriteToFile(string("beginning with customised ones]"));
  writer.WriteToFile(string("| parameter | default value | custom value |"));
  writer.WriteToFile(string("| --- | --- | --- |"));
  // iterate twice, first writing out customised defaults, then writing ones
  // that were kept at their default value
  for (size_t i(0); i < 2; ++i) {
    for (RequestedDefaults::const_iterator it = s_requesteddefaults.begin();
         it != s_requesteddefaults.end();
         ++it) {
      if ((i == 0) != (it->second.second == "")) {
        writer.WriteToFile(EncodeForMarkdown(it->first), "| ", false);
        writer.WriteToFile(EncodeForMarkdown(it->second.first), " | ", false);
        if (it->second.second == "") {
          writer.WriteToFile(string(" "), " | ", false);
        } else {
          writer.WriteToFile(EncodeForMarkdown(it->second.second), " | ", false);
        }
        writer.WriteToFile(string(""), " |", true);
      }
    }
  }
}
